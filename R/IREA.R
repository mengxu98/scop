#' Prepare an Immune Dictionary reference for IREA-style analysis
#'
#' Reads the cell-type Seurat object and signature spreadsheets distributed by
#' the Immune Dictionary portal. The experimental reference is mouse lymph-node
#' cells, including when a human orthologue table is selected.
#'
#' @param directory Directory containing the original portal downloads, or the
#'   flattened filenames used by the companion download script.
#' @param cell_type One of the portal cell-type identifiers, such as `NK_cell`.
#' @param species `"Mouse"` or `"Human"`. Human symbols in both input modes
#'   are mapped to mouse genes using the portal's signature spreadsheets.
#'   This is a partial mapping; unmapped or ambiguous symbols are omitted.
#' @return An `irea_reference` object with source paths and checksums.
#' @export
PrepareIREAReference <- function(directory, cell_type, species = c("Mouse", "Human")) {
  species <- match.arg(species)
  if (!dir.exists(directory)) stop("Reference directory does not exist.", call. = FALSE)
  valid <- c("B_cell", "cDC1", "cDC2", "Langerhans", "Macrophage",
             "MigDC", "Monocyte", "Neutrophil", "NK_cell", "pDC",
             "T_cell_CD4", "T_cell_CD8", "T_cell_gd", "Treg")
  if (!cell_type %in% valid) stop("Unsupported IREA cell type: ", cell_type, call. = FALSE)
  find_file <- function(relative) {
    candidates <- file.path(directory, c(relative, gsub("/", "_", relative, fixed = TRUE)))
    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) stop("Missing reference file: ", relative, call. = FALSE)
    normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
  }
  object_path <- find_file(paste0("downloadableData/ligands-seurat-", cell_type, ".RDS"))
  suffix <- if (species == "Human") "_Human" else ""
  cytokine_path <- find_file(paste0("dataFiles/SuppTable3_Cytokine_Signatures", suffix, ".xlsx"))
  polarization_path <- find_file(paste0("dataFiles/SuppTable7_Polarization_Signatures", suffix, ".xlsx"))
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required to read portal signature tables.", call. = FALSE)
  x <- readRDS(object_path)
  if (!inherits(x, "Seurat")) stop("The reference RDS is not a Seurat object.", call. = FALSE)
  if (!all(c("sample", "polarization") %in% colnames(x@meta.data)))
    stop("Reference object lacks sample or polarization metadata.", call. = FALSE)
  cy <- as.data.frame(readxl::read_excel(cytokine_path, sheet = cell_type))
  polar_names <- c(B_cell="B cell", cDC1="cDC1", cDC2="cDC2",
    Langerhans="Langerhans", Macrophage="Macrophage", MigDC="MigDC",
    Monocyte="Monocyte", Neutrophil="Neutrophil", NK_cell="NK cell",
    pDC="pDC", T_cell_CD4="CD4+ T cell", T_cell_CD8="CD8+ T cell",
    T_cell_gd="", Treg="Treg")
  polar_sheet <- unname(polar_names[cell_type])
  # The portal workbook stores the gamma-delta sheet using an escaped Unicode
  # title in some downloads. Its fourth sheet is stable in both species files.
  if (cell_type == "T_cell_gd")
    polar_sheet <- readxl::excel_sheets(polarization_path)[[4]]
  if (is.na(polar_sheet) || !polar_sheet %in% readxl::excel_sheets(polarization_path))
    stop("No polarization signature sheet for ", cell_type, call. = FALSE)
  po <- as.data.frame(readxl::read_excel(polarization_path, sheet = polar_sheet))
  paths <- c(object = object_path, cytokine = cytokine_path, polarization = polarization_path)
  structure(list(object = x, cytokine = cy, polarization = po,
                 cell_type = cell_type, species = species, paths = paths,
                 checksum_md5 = unname(tools::md5sum(paths)),
                 provenance = "Immune Dictionary portal; mouse in-vivo lymph-node perturbation reference"),
            class = "irea_reference")
}

.irea_sample_names <- function(reference) {
  sort(setdiff(unique(as.character(reference$object@meta.data$sample)), "PBS"))
}

.irea_layer <- function(object, assay, layer) {
  available <- SeuratObject::Layers(object, assay=assay, search=NA)
  if (length(layer) != 1L || is.na(layer) || !layer %in% available)
    stop("Select one existing expression layer by its exact name; join split layers first when needed.",
         call.=FALSE)
  SeuratObject::LayerData(object, assay=assay, layer=layer)
}

.irea_input <- function(genes, matrix, object, group_by, case, control, assay, layer, contrast) {
  contrast_inputs <- NULL
  provided <- sum(!vapply(list(genes, matrix, object), is.null, logical(1)))
  if (provided != 1L) stop("Supply exactly one of genes, matrix, or object.", call. = FALSE)
  if (!is.null(genes)) {
    if (!is.character(genes)) stop("genes must be a character vector.", call. = FALSE)
    genes <- trimws(genes[!is.na(genes) & nzchar(trimws(genes))])
    duplicates <- sum(duplicated(genes))
    genes <- unique(genes)
    if (!length(genes)) stop("No usable input genes.", call. = FALSE)
    return(list(mode="gene_list", genes=genes, duplicate_count=duplicates))
  }
  if (!is.null(object)) {
    if (!inherits(object, "Seurat")) stop("object must be Seurat.", call. = FALSE)
    if (is.null(group_by) || is.null(case) || is.null(control) ||
        !group_by %in% colnames(object@meta.data))
      stop("A Seurat input needs group_by, case, and control.", call. = FALSE)
    if (length(case) != 1L || length(control) != 1L || anyNA(c(case, control)) ||
        identical(as.character(case), as.character(control)))
      stop("case and control must be distinct, nonmissing single group names.", call. = FALSE)
    dat <- .irea_layer(object, assay, layer)
    labels <- as.character(object@meta.data[colnames(dat), group_by])
    case_cells <- which(!is.na(labels) & labels == case)
    control_cells <- which(!is.na(labels) & labels == control)
    if (!length(case_cells) || !length(control_cells))
      stop("Both groups need cells in the selected assay/layer.", call. = FALSE)
    matrix <- Matrix::rowMeans(dat[, case_cells, drop=FALSE]) -
      Matrix::rowMeans(dat[, control_cells, drop=FALSE])
    names(matrix) <- rownames(dat)
  }
  if (is.character(matrix) && length(matrix) == 1L && file.exists(matrix)) {
    ext <- tolower(tools::file_ext(matrix))
    if (ext %in% c("xlsx","xls")) {
      if (!requireNamespace("readxl",quietly=TRUE))
        stop("Package 'readxl' is required for Excel matrix files.",call.=FALSE)
      matrix <- as.data.frame(readxl::read_excel(matrix))
    } else if (ext == "csv") matrix <- utils::read.csv(matrix,check.names=FALSE)
    else if (ext == "txt") matrix <- utils::read.delim(matrix,check.names=FALSE)
    else stop("Matrix file must be .txt, .csv, .xls, or .xlsx.",call.=FALSE)
  }
  if (is.data.frame(matrix)) {
    if (ncol(matrix) < 2L) stop("Matrix input needs gene and contrast columns.", call. = FALSE)
    if (is.null(contrast) && ncol(matrix) > 2L)
      stop("Select a contrast column from this matrix.",call.=FALSE)
    column <- if (is.null(contrast)) 2L else match(contrast,names(matrix))
    if (is.na(column) || column == 1L) stop("Unknown contrast column.",call.=FALSE)
    if (anyDuplicated(names(matrix)) || any(!nzchar(names(matrix))))
      stop("Matrix columns must have unique, nonempty names.",call.=FALSE)
    contrast <- names(matrix)[column]
    contrast_inputs <- lapply(matrix[-1], function(values) {
      stats::setNames(as.numeric(values), as.character(matrix[[1]]))
    })
    v <- as.numeric(matrix[[column]])
    names(v) <- as.character(matrix[[1]])
    matrix <- v
  }
  if (is.matrix(matrix)) {
    if (ncol(matrix) != 1L) stop("Select one contrast column from the matrix.", call. = FALSE)
    v <- as.numeric(matrix[,1]); names(v) <- rownames(matrix); matrix <- v
  }
  if (!is.numeric(matrix) || is.null(names(matrix)))
    stop("matrix must be a named numeric contrast vector or a two-column data frame.", call. = FALSE)
  keep <- !is.na(names(matrix)) & nzchar(names(matrix)) & is.finite(matrix)
  matrix <- matrix[keep]
  if (anyDuplicated(names(matrix))) stop("Matrix has duplicate gene names.", call. = FALSE)
  if (!length(matrix) || all(matrix == 0)) stop("No nonzero finite gene contrasts.", call. = FALSE)
  list(mode="projection", matrix=matrix, contrast=contrast, contrast_inputs=contrast_inputs)
}

.irea_wilcox <- function(a, b) {
  if (!length(a) || !length(b)) return(NA_real_)
  if (length(unique(c(a,b))) == 1L) return(1)
  # R 4.6's two-sided asymptotic path can use 1 - pnorm(z), losing
  # small upper-tail probabilities to cancellation. Request both tails
  # directly and retain the continuity/tie corrections of stats::wilcox.test.
  tails <- suppressWarnings(vapply(c("less", "greater"), function(alternative) {
    stats::wilcox.test(a, b, alternative=alternative, exact=FALSE,
                      correct=TRUE, digits.rank=Inf)$p.value
  }, numeric(1)))
  min(1, 2 * min(tails))
}

.irea_groups <- function(reference, cells, analysis) {
  meta <- reference$object@meta.data[cells, , drop=FALSE]
  sample <- as.character(meta$sample)
  group <- if (analysis == "cytokine_response") sample else as.character(meta$polarization)
  terms <- sort(setdiff(unique(group), c("PBS", "None", "", NA_character_)))
  baseline <- which(!is.na(sample) & sample == "PBS")
  if (!length(baseline)) stop("Reference has no required PBS baseline cells.", call. = FALSE)
  if (!length(terms)) stop("Reference has no target groups.", call. = FALSE)
  list(group=group, terms=terms, baseline=baseline)
}

.irea_score_table <- function(score, groups) {
  out <- do.call(rbind, lapply(groups$terms, function(term) {
    a <- score[which(!is.na(groups$group) & groups$group == term)]
    b <- score[groups$baseline]
    data.frame(term=term, effect=mean(a)-mean(b), p_value=.irea_wilcox(a,b),
               n_target=length(a), n_control=length(b))
  }))
  out$fdr <- stats::p.adjust(out$p_value, method="BH")
  out
}

.irea_radar_score <- function(effect, fdr, cutoff=0.05) {
  positive <- ifelse(is.finite(effect), pmax(effect, 0), 0)
  significant <- !is.na(fdr) & fdr < cutoff & positive > 0
  if (any(significant) && max(positive) > 0) positive/max(positive) else rep(0, length(effect))
}

.irea_gene_score <- function(reference, genes, analysis) {
  x <- reference$object
  d <- .irea_layer(x, "RNA", "data")
  matched <- intersect(genes, rownames(d))
  if (!length(matched)) stop("No input genes match the reference.", call. = FALSE)
  score <- Matrix::colSums(d[matched,,drop=FALSE])
  groups <- .irea_groups(reference, colnames(d), analysis)
  out <- .irea_score_table(score, groups)
  list(table=out, matched_genes=matched)
}

.irea_hypergeom <- function(reference, genes, analysis) {
  tab <- if (analysis == "cytokine_response") reference$cytokine else reference$polarization
  if (analysis == "cytokine_response") {
    if (!all(c("Cytokine_Str","Gene","FDR") %in% names(tab)))
      stop("Unrecognized cytokine signature columns.", call. = FALSE)
    term <- tab$Cytokine_Str; significant <- tab$FDR < 0.01
  } else {
    if (!all(c("Polarization","Gene","P_adj") %in% names(tab)))
      stop("Unrecognized polarization signature columns.", call. = FALSE)
    term <- tab$Polarization; significant <- tab$P_adj < 0.05
  }
  universe <- rownames(reference$object)
  query <- intersect(genes, universe)
  if (!length(query)) stop("No input genes match the signature universe.", call. = FALSE)
  sets <- split(tab$Gene[!is.na(significant) & significant], term[!is.na(significant) & significant])
  terms <- if (analysis == "cytokine_response") .irea_sample_names(reference) else sort(unique(term))
  out <- do.call(rbind, lapply(terms, function(t) {
    gs <- unique(intersect(sets[[t]], universe))
    overlap <- length(intersect(gs, query))
    p <- if (!length(gs)) 1 else stats::phyper(overlap-1, length(gs),
          length(universe)-length(gs), length(query), lower.tail=FALSE)
    data.frame(term=t, effect=overlap/length(query), p_value=p, overlap=overlap,
               set_size=length(gs), query_size=length(query), universe_size=length(universe))
  }))
  out$fdr <- stats::p.adjust(out$p_value, method="BH")
  list(table=out, matched_genes=query)
}

.irea_projection <- function(reference, contrast, analysis, gene_diff_cutoff) {
  x <- reference$object
  d <- .irea_layer(x, "RNA", "data")
  candidates <- intersect(names(contrast), rownames(d))
  matched <- candidates[Matrix::rowMeans(d[candidates,,drop=FALSE]) > gene_diff_cutoff]
  if (!length(matched)) stop("No genes survive matching and the difference cutoff.", call. = FALSE)
  v <- contrast[matched]
  if (sum(v^2) == 0) stop("Projection vector has zero magnitude.", call. = FALSE)
  sub <- d[matched,,drop=FALSE]
  norm <- sqrt(Matrix::colSums(sub^2))
  projection <- as.numeric(Matrix::crossprod(sub, v))/(norm * sqrt(sum(v^2)))
  projection[!is.finite(projection)] <- 0
  groups <- .irea_groups(reference, colnames(d), analysis)
  out <- .irea_score_table(projection, groups)
  list(table=out,matched_genes=matched)
}

.irea_map_human <- function(reference, input) {
  pairs <- unique(rbind(
    reference$cytokine[,c("Gene_Human","Gene")],
    reference$polarization[,c("Gene_Human","Gene")]))
  pairs <- pairs[!is.na(pairs$Gene_Human) & nzchar(pairs$Gene_Human) &
                 !is.na(pairs$Gene) & nzchar(pairs$Gene),,drop=FALSE]
  # A human symbol mapping to multiple mouse genes cannot be resolved by the
  # published signature tables alone; omit it rather than arbitrarily choose.
  counts <- table(pairs$Gene_Human)
  pairs <- pairs[counts[pairs$Gene_Human] == 1L,,drop=FALSE]
  lookup <- stats::setNames(pairs$Gene, pairs$Gene_Human)
  if (input$mode == "gene_list") {
    keep <- input$genes %in% names(lookup)
    mapped <- unique(unname(lookup[input$genes[keep]]))
    if (!length(mapped)) stop("No unambiguous human orthologues match the reference.",call.=FALSE)
    input$genes <- mapped
  } else {
    keep <- names(input$matrix) %in% names(lookup)
    mapped <- unname(lookup[names(input$matrix)[keep]])
    val <- input$matrix[keep]
    # Multiple human genes can map to one mouse symbol. Omit all such pairs.
    unique_mouse <- !duplicated(mapped) & !duplicated(mapped,fromLast=TRUE)
    val <- val[unique_mouse]; names(val) <- mapped[unique_mouse]
    if (!length(val)) stop("No unambiguous human orthologues match the reference.",call.=FALSE)
    input$matrix <- val
  }
  input$mapping <- pairs
  input
}

#' Analyse cytokine responses or immune-cell polarization
#'
#' Experimental reconstruction of IREA. Numerical equivalence to the portal
#' has not been established for all modes. Score and projection compare target
#' reference cells against PBS-treated reference cells, including polarization.
#' These P values describe reference-cell distributions, not biological-replicate
#' inference about the user's case and control samples.
#' By default FDR is adjusted jointly across reference terms and all contrast
#' columns in a supplied data frame or file, matching the adjustment family in
#' the portal's multiple-contrast matrix example. Matrix numerical concordance
#' has not been established. Reference RNA data and the
#' selected Seurat input layer must already be appropriately normalized.
#'
#' @param reference Output from [PrepareIREAReference()].
#' @param genes Character vector for gene-list analysis.
#' @param matrix A named gene contrast vector, two-column data frame, or
#'   one-column matrix, or path to a portal-format matrix file. Positive values
#'   mean higher in the case condition.
#' @param contrast Column name to analyse when the matrix has multiple contrasts.
#' @param object Optional Seurat object used to construct a mean-expression
#'   contrast from `case` and `control`.
#' @param group_by,case,control Metadata column and group names for `object`.
#' @param assay,layer Assay and exact expression layer name for `object`. Split
#'   layers must be joined first to compare cells across those layers.
#' @param analysis `"cytokine_response"` or `"cell_polarization"`.
#' @param method `"score"` or `"hypergeometric"` for gene lists;
#'   matrix input uses cosine projection.
#' @param gene_diff_cutoff Minimum mean reference-gene expression for projection.
#' @param fdr_scope `"all_contrasts"` (default) adjusts across all contrast
#'   columns supplied, then returns the selected `contrast`. `"selected_contrast"`
#'   adjusts only its reference terms. A vector or Seurat comparison has one
#'   contrast. This argument does not change gene-list calculations.
#' @return An `irea_result`, or a Seurat object with the result in
#'   `object@tools$IREA` when `object` is supplied.
#' @export
RunIREA <- function(reference, genes=NULL, matrix=NULL, object=NULL, contrast=NULL,
                    group_by=NULL, case=NULL, control=NULL, assay="RNA", layer="data",
                    analysis=c("cytokine_response","cell_polarization"),
                    method=c("score","hypergeometric"), gene_diff_cutoff=0.25,
                    fdr_scope=c("all_contrasts","selected_contrast")) {
  if (!inherits(reference,"irea_reference")) stop("reference must be an irea_reference.",call.=FALSE)
  analysis <- match.arg(analysis); method <- match.arg(method)
  fdr_scope <- match.arg(fdr_scope)
  if (!is.numeric(gene_diff_cutoff) || length(gene_diff_cutoff)!=1L ||
      !is.finite(gene_diff_cutoff) || gene_diff_cutoff < 0)
    stop("gene_diff_cutoff must be nonnegative.",call.=FALSE)
  input <- .irea_input(genes,matrix,object,group_by,case,control,assay,layer,contrast)
  original_genes <- if (input$mode == "gene_list") input$genes else names(input$matrix)
  if (reference$species == "Human") input <- .irea_map_human(reference,input)
  core <- if (input$mode == "gene_list") {
    if (method == "score") .irea_gene_score(reference,input$genes,analysis)
    else .irea_hypergeom(reference,input$genes,analysis)
  } else .irea_projection(reference,input$matrix,analysis,gene_diff_cutoff)
  adjusted_contrasts <- input$contrast
  fdr_family <- "reference_terms_within_selected_contrast"
  if (input$mode == "projection" && fdr_scope == "all_contrasts" &&
      length(input$contrast_inputs) > 1L) {
    all_cores <- lapply(names(input$contrast_inputs), function(name) {
      if (name == input$contrast) return(core)
      other <- .irea_input(NULL,input$contrast_inputs[[name]],NULL,NULL,NULL,NULL,
                           assay,layer,NULL)
      if (reference$species == "Human") other <- .irea_map_human(reference,other)
      .irea_projection(reference,other$matrix,analysis,gene_diff_cutoff)
    })
    sizes <- vapply(all_cores, function(x) nrow(x$table), integer(1))
    family <- rep(names(input$contrast_inputs), sizes)
    adjusted <- stats::p.adjust(unlist(lapply(all_cores,function(x) x$table$p_value)),"BH")
    core$table$fdr <- adjusted[family == input$contrast]
    adjusted_contrasts <- names(input$contrast_inputs)
    fdr_family <- "reference_terms_across_all_supplied_contrasts"
  }
  table <- core$table
  table$status <- ifelse(is.na(table$p_value),"not_computable",
                         ifelse(table$fdr < 0.05,"significant","not_significant"))
  if (analysis == "cell_polarization") {
    table$radar_score <- .irea_radar_score(table$effect, table$fdr)
  }
  result <- structure(list(table=table, matched_genes=core$matched_genes,
                           input_genes=original_genes,
                           gene_mapping=if (reference$species == "Human") input$mapping else NULL,
                           parameters=list(cell_type=reference$cell_type,species=reference$species,
                             analysis=analysis,mode=input$mode,method=if(input$mode=="gene_list") method else "projection",
                             gene_diff_cutoff=gene_diff_cutoff,
                             baseline=if(input$mode=="gene_list" && method=="hypergeometric")
                               "signature_universe" else "PBS",
                             fdr_family=fdr_family, fdr_scope=fdr_scope,
                             contrasts_adjusted=adjusted_contrasts),
                           reference=list(paths=reference$paths,checksum_md5=reference$checksum_md5,
                                          provenance=reference$provenance),
                           validation="method_reconstructed; portal numeric concordance not established"),
                      class="irea_result")
  if (!is.null(object)) {object@tools$IREA <- result; return(object)}
  result
}
