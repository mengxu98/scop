#' @title Perform Gene Set Variation Analysis (GSVA)
#'
#' @md
#' @inheritParams RunEnrichment
#' @inheritParams CellDimPlot
#' @inheritParams standard_scop
#' @param group.by Name of metadata column to group cells by for averaging expression.
#' If provided, expression will be averaged within each group before GSVA analysis (cell-type level).
#' If `NULL`, GSVA is performed on each cell individually (single-cell level).
#' @param layer Data layer to use when `group.by = NULL`. Usually `"data"` for normalized or `"counts"` for count matrix.
#' Default is `"data"`.
#' @param assay_name Name of the assay to store GSVA scores when `group.by = NULL` and `new_assay = TRUE`.
#' Default is `"GSVA"`.
#' @param new_assay Whether to create a new assay for GSVA scores when `group.by = NULL`. Default is `TRUE`.
#' @param method The method to use for GSVA.
#' Options are `"gsva"`, `"ssgsea"`, `"zscore"`, or `"plage"`.
#' Default is `"gsva"`.
#' @param kcdf The kernel cumulative distribution function used for GSVA.
#' Options are `"Gaussian"` (for continuous data) or `"Poisson"` (for count data).
#' Default is `"Gaussian"`.
#' @param abs.ranking Logical indicating whether to use absolute ranking for GSVA.
#' Default is `FALSE`.
#' @param min.sz Minimum size of gene sets to be included in the analysis.
#' Default is `10`.
#' @param max.sz Maximum size of gene sets to be included in the analysis.
#' Default is `Inf`.
#' @param mx.diff Logical indicating whether to use the maximum difference method.
#' Default is `TRUE`.
#' @param tau Exponent for the GSVA method. Default is `1`.
#' @param ssgsea.norm Logical indicating whether to normalize SSGSEA scores.
#' Default is `TRUE`.
#'
#' @return
#' Returns the modified `Seurat` object. When `group.by` is provided, GSVA scores are stored in the `tools` slot.
#' When `group.by = NULL`, scores are stored in a new assay (if `new_assay = TRUE`) and in the `tools` slot.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#'
#' pancreas_sub <- RunGSVA(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   species = "Mus_musculus"
#' )
#' ht <- GSVAPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "heatmap",
#'   topTerm = 10,
#'   width = 1,
#'   height = 2
#' )
RunGSVA <- function(
  srt = NULL,
  assay = NULL,
  group.by = NULL,
  layer = "data",
  assay_name = "GSVA",
  new_assay = TRUE,
  db = "GO_BP",
  species = "Homo_sapiens",
  IDtype = "symbol",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  method = c("gsva", "ssgsea", "zscore", "plage"),
  kcdf = c("Gaussian", "Poisson"),
  abs.ranking = FALSE,
  min.sz = 10,
  max.sz = Inf,
  mx.diff = TRUE,
  tau = 1,
  ssgsea.norm = TRUE,
  verbose = TRUE
) {
  log_message("Start {.pkg GSVA} analysis", verbose = verbose)
  check_r("GSVA", verbose = FALSE)

  method <- match.arg(method)
  kcdf <- match.arg(kcdf)

  if (is.null(srt)) {
    log_message(
      "{.arg srt} must be provided for GSVA analysis",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  single_cell_mode <- is.null(group.by)
  if (single_cell_mode && isTRUE(new_assay) && is.null(assay_name)) {
    log_message(
      "{.arg assay_name} must be specified when {.arg group.by = NULL} and {.arg new_assay = TRUE}",
      message_type = "error"
    )
  }

  if (single_cell_mode) {
    log_message(
      "Single-cell GSVA mode: using expression matrix directly ...",
      verbose = verbose
    )
    expr <- GetAssayData5(srt, layer = layer, assay = assay)
    if (!inherits(expr, c("dgCMatrix", "matrix", "Matrix"))) {
      expr <- as_matrix(expr)
    }
    expr <- as_matrix(expr)
    expr <- expr[rowSums(expr) > 0, , drop = FALSE]
    if (nrow(expr) == 0 || ncol(expr) == 0) {
      log_message(
        "No expression values available for single-cell GSVA",
        message_type = "error"
      )
    }
    log_message(
      "Expression matrix: {.val {nrow(expr)}} genes x {.val {ncol(expr)}} cells",
      verbose = verbose
    )
  } else {
    if (length(group.by) != 1) {
      log_message(
        "{.arg group.by} must be a single metadata column",
        message_type = "error"
      )
    }
    if (!group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} {.val {group.by}} not found in meta.data",
        message_type = "error"
      )
    }
    log_message(
      "Averaging expression by {.val {group.by}} ...",
      verbose = verbose
    )
    features <- rownames(srt[[assay]])
    expr_avg_result <- Seurat::AggregateExpression(
      object = srt,
      features = features,
      assays = assay,
      group.by = group.by,
      verbose = FALSE
    )
    if (is.list(expr_avg_result) && length(expr_avg_result) > 0) {
      if (assay %in% names(expr_avg_result)) {
        expr <- expr_avg_result[[assay]]
      } else {
        expr <- expr_avg_result[[1]]
      }
    } else {
      expr <- expr_avg_result
    }
    if (is.null(dim(expr)) || length(dim(expr)) < 2) {
      log_message(
        "Failed to extract aggregated expression matrix. Result type: {.val {class(expr)}}, dimensions: {.val {dim(expr)}}",
        message_type = "error"
      )
    }
    if (inherits(expr, "dgCMatrix") || inherits(expr, "sparseMatrix")) {
      expr <- as.matrix(expr)
    }
    if (!is.matrix(expr)) {
      expr <- as.matrix(expr)
    }
    expr <- expr[rowSums(expr) > 0, , drop = FALSE]
    if (nrow(expr) == 0 || ncol(expr) == 0) {
      log_message(
        "No aggregated expression values available for {.val {group.by}}",
        message_type = "error"
      )
    }
    log_message(
      "Aggregated expression matrix: {.val {nrow(expr)}} genes x {.val {ncol(expr)}} groups",
      verbose = verbose
    )
  }

  if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species,
      db = db,
      db_update = db_update,
      db_version = db_version,
      db_IDtypes = IDtype,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    db <- "custom"
    db_list <- list()
    db_list[[species]] <- list()
    db_list[[species]][[db]] <- list()
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- "custom"
  }

  if (isTRUE(db_combine)) {
    log_message(
      "Create {.val Combined} database ...",
      verbose = verbose
    )
    TERM2GENE <- do.call(
      rbind,
      lapply(
        db_list[[species]],
        function(x) x[["TERM2GENE"]][, c("Term", IDtype)]
      )
    )
    TERM2NAME <- do.call(
      rbind,
      lapply(names(db_list[[species]]), function(x) {
        db_list[[species]][[x]][["TERM2NAME"]][["Name"]] <- paste0(
          db_list[[species]][[x]][["TERM2NAME"]][["Name"]],
          " [",
          x,
          "]"
        )
        db_list[[species]][[x]][["TERM2NAME"]][, c("Term", "Name")]
      })
    )
    version <- unlist(lapply(
      db_list[[species]],
      function(x) as.character(x[["version"]])
    ))
    version <- paste0(names(version), ":", version, collapse = ";")
    db <- "Combined"
    db_list[[species]][[db]] <- list()
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- unique(version)
  }

  gsva_results <- list()
  enrichment_results <- list()

  for (term in db) {
    log_message(
      "Processing database: {.val {term}} ...",
      verbose = verbose
    )

    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c(
      "Term",
      IDtype
    )]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]

    dup <- duplicated(TERM2GENE_tmp)
    na <- Matrix::rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[
      TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], ,
      drop = FALSE
    ]

    gene_sets <- split(
      TERM2GENE_tmp[[IDtype]],
      TERM2GENE_tmp[["Term"]]
    )
    gs_size <- lengths(gene_sets)
    min_size <- ifelse(term %in% unlimited_db, 1, max(minGSSize, min.sz))
    max_size <- ifelse(term %in% unlimited_db, Inf, min(maxGSSize, max.sz))
    gene_sets <- gene_sets[gs_size >= min_size & gs_size <= max_size]

    if (length(gene_sets) == 0) {
      log_message(
        "No gene sets found for {.val {term}} after filtering",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    genes_in_sets <- unique(unlist(gene_sets))
    overlap <- intersect(rownames(expr), genes_in_sets)
    log_message(
      "Initial overlap: {.val {length(overlap)}} genes out of {.val {nrow(expr)}} expression genes and {.val {length(genes_in_sets)}} genes in gene sets",
      verbose = verbose
    )

    if (length(overlap) == 0) {
      log_message(
        "No overlapping genes found for {.val {term}}",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    expr_filtered <- expr[overlap, , drop = FALSE]
    gene_sets_filtered <- lapply(gene_sets, function(gs) {
      intersect(gs, rownames(expr_filtered))
    })
    keep <- lengths(gene_sets_filtered) >= min_size
    gene_sets_filtered <- gene_sets_filtered[keep]

    if (length(gene_sets_filtered) == 0) {
      log_message(
        "No gene sets remaining after filtering for {.val {term}}",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    log_message(
      "Running GSVA for {.val {length(gene_sets_filtered)}} gene sets ...",
      verbose = verbose
    )

    if (method == "gsva") {
      param <- GSVA::gsvaParam(
        exprData = expr_filtered,
        geneSets = gene_sets_filtered,
        minSize = min_size,
        maxSize = max_size,
        kcdf = kcdf,
        tau = tau,
        maxDiff = mx.diff,
        absRanking = abs.ranking
      )
    } else if (method == "ssgsea") {
      param <- GSVA::ssgseaParam(
        exprData = expr_filtered,
        geneSets = gene_sets_filtered,
        minSize = min_size,
        maxSize = max_size,
        alpha = tau,
        normalize = ssgsea.norm,
        verbose = verbose
      )
    } else if (method == "zscore") {
      param <- GSVA::zscoreParam(
        exprData = expr_filtered,
        geneSets = gene_sets_filtered,
        minSize = min_size,
        maxSize = max_size
      )
    } else if (method == "plage") {
      param <- GSVA::plageParam(
        exprData = expr_filtered,
        geneSets = gene_sets_filtered,
        minSize = min_size,
        maxSize = max_size
      )
    }

    gsva_scores <- GSVA::gsva(
      param = param,
      verbose = verbose
    )

    if (inherits(gsva_scores, "SummarizedExperiment")) {
      gsva_scores <- SummarizedExperiment::assay(gsva_scores)
    }
    if (!is.matrix(gsva_scores)) {
      gsva_scores <- as.matrix(gsva_scores)
    }

    term_ids <- rownames(gsva_scores)
    term_names <- TERM2NAME_tmp[
      match(term_ids, TERM2NAME_tmp[["Term"]]),
      "Name"
    ]
    term_names[is.na(term_names)] <- term_ids[is.na(term_names)]
    term_names <- thisutils::capitalize(
      trimws(as.character(term_names)),
      force_tolower = TRUE
    )
    term_gene_ids <- vapply(
      term_ids,
      function(x) {
        genes <- TERM2GENE_tmp[TERM2GENE_tmp[["Term"]] %in% x, IDtype]
        genes <- unique(as.character(genes[!is.na(genes)]))
        paste0(genes, collapse = "/")
      },
      FUN.VALUE = character(1)
    )
    rownames(gsva_scores) <- term_names

    gsva_results[[term]] <- gsva_scores
    enrichment_results[[term]] <- do.call(rbind, lapply(colnames(gsva_scores), function(group) {
      data.frame(
        ID = term_ids,
        Description = term_names,
        geneID = term_gene_ids,
        Groups = group,
        Database = term,
        Version = as.character(db_list[[species]][[term]][["version"]]),
        GSVA_Score = as.numeric(gsva_scores[, group]),
        stringsAsFactors = FALSE
      )
    }))
  }

  if (length(gsva_results) == 0) {
    log_message(
      "No GSVA results generated",
      message_type = "error"
    )
  }

  gsva_scores_combined <- if (length(gsva_results) > 1) {
    do.call(rbind, gsva_results)
  } else {
    gsva_results[[1]]
  }
  enrichment <- do.call(rbind, enrichment_results)
  rownames(enrichment) <- NULL

  res <- list(
    enrichment = enrichment,
    results = gsva_results,
    scores = gsva_scores_combined,
    input = expr,
    group.by = group.by,
    method = method,
    kcdf = kcdf,
    db = unique(enrichment[["Database"]]),
    species = species
  )

  tool_name <- if (single_cell_mode) {
    paste("GSVA_cell", method, sep = "_")
  } else {
    paste("GSVA", group.by, method, sep = "_")
  }
  srt@tools[[tool_name]] <- res

  if (single_cell_mode) {
    scores_mat <- Matrix::t(as_matrix(gsva_scores_combined))
    scores_mat <- scores_mat[intersect(rownames(scores_mat), colnames(srt)), , drop = FALSE]
    scores_mat <- scores_mat[colnames(srt), , drop = FALSE]
    if (isTRUE(new_assay)) {
      srt[[assay_name]] <- Seurat::CreateAssayObject(
        counts = Matrix::t(as_matrix(scores_mat))
      )
      srt[[assay_name]] <- Seurat::AddMetaData(
        object = srt[[assay_name]],
        metadata = data.frame(
          termnames = rownames(gsva_scores_combined),
          row.names = rownames(gsva_scores_combined)
        )
      )
      log_message(
        "{.pkg GSVA} results stored in assay {.val {assay_name}} and tools slot {.val {tool_name}}",
        verbose = verbose
      )
    } else {
      srt <- Seurat::AddMetaData(object = srt, metadata = as.data.frame(scores_mat))
      log_message(
        "{.pkg GSVA} results stored in meta.data and tools slot {.val {tool_name}}",
        verbose = verbose
      )
    }
  } else {
    log_message(
      "{.pkg GSVA} results stored in {.code tools} slot: {.val {tool_name}}",
      verbose = verbose
    )
  }

  log_message(
    "{.pkg GSVA} analysis done",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}
