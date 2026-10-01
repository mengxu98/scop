irea_layer <- function(object, assay, layer) {
  available <- SeuratObject::Layers(object, assay = assay, search = NA)
  if (length(layer) != 1L || is.na(layer) || !layer %in% available) {
    stop("Select one existing expression layer by its exact name; join split layers first when needed.",
      call. = FALSE
    )
  }
  SeuratObject::LayerData(object, assay = assay, layer = layer)
}

irea_input <- function(object, contrast = NULL) {
  contrast_inputs <- NULL
  source_file <- NULL
  if (is.character(object)) {
    if (length(object) == 1L && !is.na(object) &&
      tolower(tools::file_ext(object)) %in% c("txt", "csv", "xls", "xlsx")) {
      if (!file.exists(object)) stop("Matrix file does not exist.", call. = FALSE)
      source_file <- normalizePath(object, winslash = "/", mustWork = TRUE)
      ext <- tolower(tools::file_ext(object))
      if (ext %in% c("xlsx", "xls")) {
        thisutils::check_r("readxl", install = FALSE, verbose = FALSE)
        object <- as.data.frame(readxl::read_excel(object))
      } else if (ext == "csv") {
        object <- utils::read.csv(object, check.names = FALSE)
      } else {
        object <- utils::read.delim(object, check.names = FALSE)
      }
    } else {
      genes <- trimws(object[!is.na(object) & nzchar(trimws(object))])
      duplicates <- sum(duplicated(genes))
      genes <- unique(genes)
      if (!length(genes)) stop("No usable input genes.", call. = FALSE)
      return(list(mode = "gene_list", genes = genes, duplicate_count = duplicates))
    }
  }
  if (is.data.frame(object)) {
    if (ncol(object) < 2L) stop("Matrix input needs gene and contrast columns.", call. = FALSE)
    if (is.null(contrast) && ncol(object) > 2L) {
      stop("Select a contrast column from this matrix.", call. = FALSE)
    }
    if (anyDuplicated(names(object)) || anyNA(names(object)) || any(!nzchar(names(object)))) {
      stop("Matrix columns must have unique, nonempty names.", call. = FALSE)
    }
    if (!all(vapply(object[-1], is.numeric, logical(1)))) {
      stop("Contrast columns must be numeric; factors and text are not accepted.", call. = FALSE)
    }
    if (!is.null(contrast) && (length(contrast) != 1L || is.na(contrast))) {
      stop("Select one contrast column by name.", call. = FALSE)
    }
    column <- if (is.null(contrast)) 2L else match(contrast, names(object))
    if (is.na(column) || column == 1L) stop("Unknown contrast column.", call. = FALSE)
    contrast <- names(object)[column]
    contrast_inputs <- lapply(object[-1], function(values) {
      stats::setNames(values, as.character(object[[1]]))
    })
    object <- contrast_inputs[[contrast]]
  }
  if (is.matrix(object)) {
    if (ncol(object) != 1L) stop("Select one contrast column from the matrix.", call. = FALSE)
    if (!is.numeric(object)) stop("Contrast matrix must be numeric.", call. = FALSE)
    if (is.null(contrast) && !is.null(colnames(object))) contrast <- colnames(object)[[1]]
    object <- stats::setNames(object[, 1], rownames(object))
  }
  if (!is.numeric(object) || is.null(names(object))) {
    stop("object must be genes, a named numeric contrast, a contrast table/file, or Seurat.", call. = FALSE)
  }
  keep <- !is.na(names(object)) & nzchar(names(object)) & is.finite(object)
  object <- object[keep]
  if (anyDuplicated(names(object))) stop("Matrix has duplicate gene names.", call. = FALSE)
  if (!length(object) || all(object == 0)) stop("No nonzero finite gene contrasts.", call. = FALSE)
  list(
    mode = "projection", matrix = object, contrast = contrast,
    contrast_inputs = contrast_inputs, source_file = source_file
  )
}

irea_wilcox <- function(a, b) {
  if (!length(a) || !length(b)) {
    return(NA_real_)
  }
  if (length(unique(c(a, b))) == 1L) {
    return(1)
  }
  tails <- suppressWarnings(vapply(c("less", "greater"), function(alternative) {
    stats::wilcox.test(a, b,
      alternative = alternative, exact = FALSE,
      correct = TRUE, digits.rank = Inf
    )$p.value
  }, numeric(1)))
  min(1, 2 * min(tails))
}

irea_groups <- function(reference, cells, analysis) {
  meta <- reference$object@meta.data[cells, , drop = FALSE]
  sample <- as.character(meta$sample)
  group <- if (analysis == "cytokine_response") sample else as.character(meta$polarization)
  terms <- sort(setdiff(unique(group), c("PBS", "None", "", NA_character_)))
  baseline <- which(!is.na(sample) & sample == "PBS")
  if (!length(baseline)) stop("Reference has no required PBS baseline cells.", call. = FALSE)
  if (!length(terms)) stop("Reference has no target groups.", call. = FALSE)
  list(group = group, terms = terms, baseline = baseline)
}

irea_score_table <- function(score, groups) {
  out <- do.call(rbind, lapply(groups$terms, function(term) {
    a <- score[which(!is.na(groups$group) & groups$group == term)]
    b <- score[groups$baseline]
    data.frame(
      term = term, effect = mean(a) - mean(b), p_value = irea_wilcox(a, b),
      n_target = length(a), n_control = length(b)
    )
  }))
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

irea_radar_score <- function(effect, fdr, cutoff = 0.05) {
  positive <- ifelse(is.finite(effect), pmax(effect, 0), 0)
  significant <- !is.na(fdr) & fdr < cutoff & positive > 0
  if (any(significant) && max(positive) > 0) positive / max(positive) else rep(0, length(effect))
}

irea_gene_score <- function(reference, genes, analysis) {
  x <- reference$object
  d <- irea_layer(x, "RNA", "data")
  matched <- intersect(genes, rownames(d))
  if (!length(matched)) stop("No input genes match the reference.", call. = FALSE)
  score <- Matrix::colSums(d[matched, , drop = FALSE])
  groups <- irea_groups(reference, colnames(d), analysis)
  out <- irea_score_table(score, groups)
  list(table = out, matched_genes = matched)
}

irea_hypergeom <- function(reference, genes, analysis) {
  tab <- if (analysis == "cytokine_response") reference$cytokine else reference$polarization
  if (analysis == "cytokine_response") {
    if (!all(c("Cytokine_Str", "Gene", "FDR") %in% names(tab))) {
      stop("Unrecognized cytokine signature columns.", call. = FALSE)
    }
    term <- tab$Cytokine_Str
    significant <- tab$FDR < 0.01
  } else {
    if (!all(c("Polarization", "Gene", "P_adj") %in% names(tab))) {
      stop("Unrecognized polarization signature columns.", call. = FALSE)
    }
    term <- tab$Polarization
    significant <- tab$P_adj < 0.05
  }
  universe <- rownames(reference$object)
  query <- intersect(genes, universe)
  if (!length(query)) stop("No input genes match the signature universe.", call. = FALSE)
  sets <- split(tab$Gene[!is.na(significant) & significant], term[!is.na(significant) & significant])
  terms <- if (analysis == "cytokine_response") sort(setdiff(unique(as.character(reference$object@meta.data$sample)), "PBS")) else sort(unique(term))
  out <- do.call(rbind, lapply(terms, function(t) {
    gs <- unique(intersect(sets[[t]], universe))
    overlap <- length(intersect(gs, query))
    p <- if (!length(gs)) {
      1
    } else {
      stats::phyper(overlap - 1, length(gs),
        length(universe) - length(gs), length(query),
        lower.tail = FALSE
      )
    }
    data.frame(
      term = t, effect = overlap / length(query), p_value = p, overlap = overlap,
      set_size = length(gs), query_size = length(query), universe_size = length(universe)
    )
  }))
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  list(table = out, matched_genes = query)
}

irea_projection <- function(reference, contrast, analysis, gene_diff_cutoff) {
  x <- reference$object
  d <- irea_layer(x, "RNA", "data")
  candidates <- intersect(names(contrast), rownames(d))
  matched <- candidates[Matrix::rowMeans(d[candidates, , drop = FALSE]) > gene_diff_cutoff]
  if (!length(matched)) stop("No genes survive matching and the difference cutoff.", call. = FALSE)
  v <- contrast[matched]
  if (sum(v^2) == 0) stop("Projection vector has zero magnitude.", call. = FALSE)
  sub <- d[matched, , drop = FALSE]
  norm <- sqrt(Matrix::colSums(sub^2))
  projection <- as.numeric(Matrix::crossprod(sub, v)) / (norm * sqrt(sum(v^2)))
  projection[!is.finite(projection)] <- 0
  groups <- irea_groups(reference, colnames(d), analysis)
  out <- irea_score_table(projection, groups)
  list(table = out, matched_genes = matched)
}

irea_map_human <- function(reference, input) {
  pairs <- unique(rbind(
    reference$cytokine[, c("Gene_Human", "Gene")],
    reference$polarization[, c("Gene_Human", "Gene")]
  ))
  pairs <- pairs[!is.na(pairs$Gene_Human) & nzchar(pairs$Gene_Human) &
    !is.na(pairs$Gene) & nzchar(pairs$Gene), , drop = FALSE]
  counts <- table(pairs$Gene_Human)
  pairs <- pairs[counts[pairs$Gene_Human] == 1L, , drop = FALSE]
  lookup <- stats::setNames(pairs$Gene, pairs$Gene_Human)
  if (input$mode == "gene_list") {
    keep <- input$genes %in% names(lookup)
    mapped <- unique(unname(lookup[input$genes[keep]]))
    if (!length(mapped)) stop("No unambiguous human orthologues match the reference.", call. = FALSE)
    input$genes <- mapped
  } else {
    keep <- names(input$matrix) %in% names(lookup)
    mapped <- unname(lookup[names(input$matrix)[keep]])
    val <- input$matrix[keep]
    unique_mouse <- !duplicated(mapped) & !duplicated(mapped, fromLast = TRUE)
    val <- val[unique_mouse]
    names(val) <- mapped[unique_mouse]
    if (!length(val)) stop("No unambiguous human orthologues match the reference.", call. = FALSE)
    input$matrix <- val
  }
  input$mapping <- pairs
  input
}

#' @title Analyse cytokine responses or immune-cell polarization
#'
#' @description Experimental reconstruction of published IREA score,
#' hypergeometric and cosine-projection analyses. Numerical equivalence to the
#' portal has not been established for all modes.
#' @details Score and projection compare target reference cells against
#' PBS-treated reference cells, including polarization. Their P values describe
#' reference-cell distributions, not biological-replicate inference about the
#' user's samples. Reference RNA data and the selected input layer must already
#' be normalized. Human symbols are mapped to mouse reference genes using the
#' signature tables; ambiguous and unmapped symbols are omitted.
#' @md
#' @param object Character gene vector, named numeric gene contrast vector,
#'   numeric one-column matrix (genes in rows), data frame with genes in its first
#'   column and numeric contrast columns, path to a .txt/.csv/.xls/.xlsx table,
#'   or a Seurat object. Positive contrasts mean higher expression in the case.
#'   Factors and text contrast columns are rejected. Missing/empty gene names and
#'   nonfinite contrasts are omitted; duplicated contrast gene names are rejected.
#' @param reference An `irea_reference` from [PrepareDB()] or [PrepareIREAReference()].
#' @param contrast Column name to analyse in a table with multiple contrasts.
#' @param group.by,case,control Metadata column and distinct group names for Seurat.
#' @param assay,layer Seurat assay and exact expression layer name. Split layers
#'   must be joined first to compare cells across those layers.
#' @param analysis `"cytokine_response"` or `"cell_polarization"`.
#' @param method `"score"` or `"hypergeometric"` for gene lists;
#'   numeric contrasts and Seurat inputs use cosine projection.
#' @param gene_diff_cutoff Nonnegative minimum mean reference-gene expression for
#'   projection. Genes must exceed this value to be retained.
#' @param fdr_scope `"all_contrasts"` adjusts BH jointly across reference terms
#'   and all supplied table columns before returning the selected contrast.
#'   `"selected_contrast"` adjusts its terms only. Vectors and Seurat inputs have
#'   one contrast. This option does not change gene-list calculations.
#' @param ... Arguments passed to the dispatched method.
#' @return An `irea_result` list, or the Seurat object with that result stored in
#'   `object@tools$IREA`. Other tools are preserved. The result contains:
#'   * `table`: one row per reference term, sorted by term. `effect` is the mean
#'     target-minus-PBS summed expression (score) or cosine projection difference,
#'     or the fraction of query genes in a significant signature (hypergeometric).
#'     `p_value` is a two-sided asymptotic Wilcoxon P value with continuity/tie
#'     correction, or an upper-tail hypergeometric P value; `fdr` is BH-adjusted.
#'     Empty target/control samples yield NA P/FDR. `status` uses FDR < 0.05.
#'     Score/projection include `n_target` and `n_control`; hypergeometric includes
#'     `overlap`, `set_size`, `query_size` and `universe_size`.
#'     Polarization adds `radar_score` between 0 and 1: zero if no positive effect has
#'     FDR < 0.05, otherwise each positive effect divided by the largest positive
#'     effect, with negative effects set to zero.
#'   * `matched_genes`, `input_genes`: gene identifiers in input order after
#'     matching/filtering and before mapping, respectively.
#'   * `gene_mapping`: human/mouse signature pairs for human input, otherwise NULL.
#'   * `parameters`: analysis choices, exact selected contrast and FDR family;
#'     Seurat results additionally record grouping, case/control and assay/layer.
#'   * `reference`: source paths, corresponding checksums and provenance.
#'   * `validation`: experimental numerical-validation status.
#' @inherit PrepareIREAReference references
#' @seealso [IREAPlot], [PrepareDB]
#' @export
#' @examples
#' \dontrun{
#' reference <- PrepareDB(
#'   db = "IREA_Macrophage", species = "Homo_sapiens"
#' )[["Homo_sapiens"]][["IREA_Macrophage"]]
#' data(panc8_sub)
#' counts <- SeuratObject::LayerData(panc8_sub, assay = "RNA", layer = "counts")
#' macrophage <- which(panc8_sub$celltype == "macrophage")
#' genes <- intersect(
#'   c("ISG15", "IFIT3", "BST2"),
#'   rownames(counts)[Matrix::rowSums(counts[, macrophage, drop = FALSE]) > 0]
#' )
#' result <- RunIREA(genes, reference = reference)
#' IREAPlot(result)
#' }
RunIREA <- function(object, ...) {
  UseMethod("RunIREA", object)
}

#' @rdname RunIREA
#' @export
RunIREA.default <- function(object, reference, contrast = NULL,
                            analysis = c("cytokine_response", "cell_polarization"),
                            method = c("score", "hypergeometric"), gene_diff_cutoff = 0.25,
                            fdr_scope = c("all_contrasts", "selected_contrast"), ...) {
  if (length(list(...))) stop("Unused RunIREA arguments.", call. = FALSE)
  if (!inherits(reference, "irea_reference")) stop("reference must be an irea_reference.", call. = FALSE)
  analysis <- match.arg(analysis)
  method <- match.arg(method)
  fdr_scope <- match.arg(fdr_scope)
  if (!is.numeric(gene_diff_cutoff) || length(gene_diff_cutoff) != 1L ||
    !is.finite(gene_diff_cutoff) || gene_diff_cutoff < 0) {
    stop("gene_diff_cutoff must be nonnegative.", call. = FALSE)
  }
  input <- irea_input(object, contrast)
  original_genes <- if (input$mode == "gene_list") input$genes else names(input$matrix)
  if (reference$species == "Human") input <- irea_map_human(reference, input)
  core <- if (input$mode == "gene_list") {
    if (method == "score") {
      irea_gene_score(reference, input$genes, analysis)
    } else {
      irea_hypergeom(reference, input$genes, analysis)
    }
  } else {
    irea_projection(reference, input$matrix, analysis, gene_diff_cutoff)
  }
  adjusted_contrasts <- input$contrast
  fdr_family <- "reference_terms_within_selected_contrast"
  if (input$mode == "projection" && fdr_scope == "all_contrasts" &&
    length(input$contrast_inputs) > 1L) {
    all_cores <- lapply(names(input$contrast_inputs), function(name) {
      if (name == input$contrast) {
        return(core)
      }
      other <- irea_input(input$contrast_inputs[[name]])
      if (reference$species == "Human") other <- irea_map_human(reference, other)
      irea_projection(reference, other$matrix, analysis, gene_diff_cutoff)
    })
    sizes <- vapply(all_cores, function(x) nrow(x$table), integer(1))
    family <- rep(names(input$contrast_inputs), sizes)
    adjusted <- stats::p.adjust(unlist(lapply(all_cores, function(x) x$table$p_value)), "BH")
    core$table$fdr <- adjusted[family == input$contrast]
    adjusted_contrasts <- names(input$contrast_inputs)
    fdr_family <- "reference_terms_across_all_supplied_contrasts"
  }
  table <- core$table
  table$status <- ifelse(is.na(table$p_value), "not_computable",
    ifelse(table$fdr < 0.05, "significant", "not_significant")
  )
  if (analysis == "cell_polarization") {
    table$radar_score <- irea_radar_score(table$effect, table$fdr)
  }
  result <- structure(
    list(
      table = table, matched_genes = core$matched_genes,
      input_genes = original_genes,
      gene_mapping = if (reference$species == "Human") input$mapping else NULL,
      parameters = list(
        cell_type = reference$cell_type, species = reference$species,
        analysis = analysis, mode = input$mode, method = if (input$mode == "gene_list") method else "projection",
        gene_diff_cutoff = gene_diff_cutoff,
        contrast = input$contrast, source_file = input$source_file,
        baseline = if (input$mode == "gene_list" && method == "hypergeometric") {
          "signature_universe"
        } else {
          "PBS"
        },
        fdr_family = fdr_family, fdr_scope = fdr_scope,
        contrasts_adjusted = adjusted_contrasts
      ),
      reference = list(
        paths = reference$paths, checksum_md5 = reference$checksum_md5,
        provenance = reference$provenance
      ),
      validation = "method_reconstructed; portal numeric concordance not established"
    ),
    class = "irea_result"
  )
  result
}

#' @rdname RunIREA
#' @export
RunIREA.Seurat <- function(object, reference, group.by, case, control,
                           assay = "RNA", layer = "data", ...) {
  if (length(group.by) != 1L || is.na(group.by) ||
    !group.by %in% colnames(object@meta.data)) {
    stop("A Seurat input needs an existing group.by metadata column.", call. = FALSE)
  }
  if (length(case) != 1L || length(control) != 1L || anyNA(c(case, control)) ||
    identical(as.character(case), as.character(control))) {
    stop("case and control must be distinct, nonmissing single group names.", call. = FALSE)
  }
  dat <- irea_layer(object, assay, layer)
  labels <- as.character(object@meta.data[colnames(dat), group.by])
  case_cells <- which(!is.na(labels) & labels == case)
  control_cells <- which(!is.na(labels) & labels == control)
  if (!length(case_cells) || !length(control_cells)) {
    stop("Both groups need cells in the selected assay/layer.", call. = FALSE)
  }
  values <- Matrix::rowMeans(dat[, case_cells, drop = FALSE]) -
    Matrix::rowMeans(dat[, control_cells, drop = FALSE])
  names(values) <- rownames(dat)
  result <- RunIREA.default(values, reference = reference, ...)
  result$parameters <- c(result$parameters, list(
    group.by = group.by,
    case = as.character(case), control = as.character(control), assay = assay, layer = layer
  ))
  object@tools$IREA <- result
  object
}
