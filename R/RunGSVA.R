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
#' @param store_metadata Whether to also store single-cell GSVA scores in `meta.data`.
#' When `NULL`, custom `features` or `TERM2GENE` input is stored in `meta.data`
#' by default, while database-derived results stay assay-only when `new_assay = TRUE`.
#' @param method The method to use for GSVA.
#' Options are `"gsva"`, `"ssgsea"`, `"zscore"`, or `"plage"`.
#' Multiple methods can be supplied at once; in single-cell mode they will be
#' stored in method-suffixed assays such as `"GSVA_gsva"` and `"GSVA_ssgsea"`.
#' Default is `"gsva"`.
#' @param backend Scoring backend. `"cpp"` is the default and supports all
#' current `method` values. `"r"` uses the original [GSVA::gsva()]
#' implementation. `"cpp"` supports `method = "ssgsea"`,
#' `method = "zscore"`, `method = "plage"`, and `method = "gsva"` with
#' `kcdf = "Gaussian"`, `kcdf = "Poisson"`, or `kcdf = "none"`. Gaussian
#' GSVA uses the native C++ KDE/ranking kernel; the other GSVA kernels retain
#' the validated GSVA implementation. PLAGE scores are oriented to have non-negative
#' dot product with the gene set mean z-score so SVD signs are deterministic.
#' @param cpp_chunk_size Optional cell chunk size for C++ GSVA kernels. `NULL`
#' or `"auto"` automatically chunks large matrices to reduce peak dense
#' intermediate memory; positive values set the chunk size manually.
#' @param kcdf The kernel cumulative distribution function used for GSVA.
#' Options are `"Gaussian"` (for continuous data), `"Poisson"` (for count data),
#' or `"none"` (skip kernel estimation and use ranks directly).
#' When omitted, `backend = "cpp"` with `method = "gsva"` uses `"none"` for
#' faster single-cell scoring; explicit `"Gaussian"` or `"Poisson"` values are
#' still honored. Other backends and methods default to `"Gaussian"`.
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
#' When `group.by = NULL`, scores are stored in the `tools` slot, optionally in a new assay,
#' and optionally in `meta.data` for direct use with [FeatureDimPlot()] and [FeatureStatPlot()].
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
#'
#' features_all <- rownames(pancreas_sub)
#' pancreas_sub <- RunGSVA(
#'   pancreas_sub,
#'   features = list(
#'     A = features_all[1:20],
#'     B = features_all[21:40]
#'   ),
#'   method = c("gsva", "ssgsea")
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "GSVA_gsva_A",
#'   add_density = TRUE
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("GSVA_gsva_A", "GSVA_ssgsea_A"),
#'   group.by = "CellType",
#'   plot.by = "feature",
#'   plot_type = "violin",
#'   stack = TRUE,
#'   flip = TRUE
#' )
RunGSVA <- function(
  srt = NULL,
  assay = NULL,
  group.by = NULL,
  layer = "data",
  assay_name = "GSVA",
  new_assay = TRUE,
  store_metadata = NULL,
  db = "GO_BP",
  species = "Homo_sapiens",
  IDtype = "symbol",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  features = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  method = c("gsva", "ssgsea", "zscore", "plage"),
  backend = c("cpp", "r"),
  cpp_chunk_size = NULL,
  kcdf = c("Gaussian", "Poisson", "none"),
  abs.ranking = FALSE,
  min.sz = 10,
  max.sz = Inf,
  mx.diff = TRUE,
  tau = 1,
  ssgsea.norm = TRUE,
  verbose = TRUE,
  ...
) {
  log_message("Start {.pkg GSVA} analysis", verbose = verbose)

  kcdf_defaulted <- missing(kcdf)
  method <- unique(tolower(as.character(method)))
  backend <- match.arg(backend)
  valid_methods <- c("gsva", "ssgsea", "zscore", "plage")
  if (!all(method %in% valid_methods)) {
    log_message(
      "{.arg method} must be one of {.val {valid_methods}}",
      message_type = "error"
    )
  }
  single_cell_mode <- is.null(group.by)
  custom_input_supplied <- !is.null(features) || !is.null(TERM2GENE)
  if (is.null(store_metadata)) {
    store_metadata <- isTRUE(custom_input_supplied) || !isTRUE(new_assay)
  } else {
    if (!is.logical(store_metadata) || length(store_metadata) != 1L || is.na(store_metadata)) {
      log_message(
        "{.arg store_metadata} must be TRUE, FALSE, or NULL",
        message_type = "error"
      )
    }
  }
  if (length(method) > 1L) {
    for (method_i in method) {
      assay_name_i <- if (single_cell_mode) {
        paste(assay_name, method_i, sep = "_")
      } else {
        assay_name
      }
      srt <- RunGSVA(
        srt = srt,
        assay = assay,
        group.by = group.by,
        layer = layer,
        assay_name = assay_name_i,
        new_assay = new_assay,
        store_metadata = store_metadata,
        db = db,
        species = species,
        IDtype = IDtype,
        db_update = db_update,
        db_version = db_version,
        db_combine = db_combine,
        convert_species = convert_species,
        Ensembl_version = Ensembl_version,
        mirror = mirror,
        features = features,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        unlimited_db = unlimited_db,
        method = method_i,
        backend = backend,
        cpp_chunk_size = cpp_chunk_size,
        kcdf = if (isTRUE(kcdf_defaulted) && identical(backend, "cpp") && identical(method_i, "gsva")) "none" else kcdf,
        abs.ranking = abs.ranking,
        min.sz = min.sz,
        max.sz = max.sz,
        mx.diff = mx.diff,
        tau = tau,
        ssgsea.norm = ssgsea.norm,
        verbose = verbose,
        ...
      )
    }
    if (single_cell_mode && isTRUE(new_assay) && !assay_name %in% SeuratObject::Assays(srt)) {
      first_assay_name <- paste(assay_name, method[[1]], sep = "_")
      if (first_assay_name %in% SeuratObject::Assays(srt)) {
        srt[[assay_name]] <- srt[[first_assay_name]]
      }
    }
    return(srt)
  }
  method <- match.arg(method)
  kcdf <- match.arg(kcdf)
  if (isTRUE(kcdf_defaulted) && identical(backend, "cpp") && identical(method, "gsva")) {
    kcdf <- "none"
  }
  if (identical(backend, "cpp")) {
    if (!method %in% c("gsva", "ssgsea", "zscore", "plage")) {
      log_message(
        "{.arg backend = 'cpp'} currently supports {.arg method = 'gsva'}, {.arg method = 'ssgsea'}, {.arg method = 'zscore'}, and {.arg method = 'plage'} only",
        message_type = "error"
      )
    }
  } else {
    gene_set_scoring_require_namespace("GSVA")
  }

  if (is.null(srt)) {
    log_message(
      "{.arg srt} must be provided for GSVA analysis",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
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
    expr_row_sums <- if (inherits(expr, "Matrix")) Matrix::rowSums(expr) else rowSums(expr)
    expr <- expr[expr_row_sums > 0, , drop = FALSE]
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
    assay_features <- rownames(srt[[assay]])
    expr <- NULL
  }

  if (!is.null(features)) {
    db <- "custom"
    custom_db <- create_custom_db_list_from_features(
      species = species,
      db = db,
      features = features,
      IDtype = IDtype,
      version = "custom"
    )
    db_list <- custom_db[["db_list"]]
    TERM2GENE <- custom_db[["TERM2GENE"]]
    TERM2NAME <- custom_db[["TERM2NAME"]]
  } else if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species,
      db = db,
      db_update = db_update,
      db_version = db_version,
      db_IDtypes = IDtype,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      ...
    )
  } else {
    db <- "custom"
    custom_db <- create_custom_db_list(
      species = species,
      db = db,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      IDtype = IDtype,
      version = "custom"
    )
    db_list <- custom_db[["db_list"]]
    TERM2GENE <- custom_db[["TERM2GENE"]]
    TERM2NAME <- custom_db[["TERM2NAME"]]
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

  if (!single_cell_mode) {
    log_message(
      "Averaging expression by {.val {group.by}} ...",
      verbose = verbose
    )
    expr_avg_result <- Seurat::AggregateExpression(
      object = srt,
      features = assay_features,
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
    expr_row_sums <- if (inherits(expr, "Matrix")) Matrix::rowSums(expr) else rowSums(expr)
    expr <- expr[expr_row_sums > 0, , drop = FALSE]
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

    keep_term_gene <- !duplicated(TERM2GENE_tmp) &
      stats::complete.cases(TERM2GENE_tmp)
    TERM2GENE_tmp <- TERM2GENE_tmp[keep_term_gene, , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[
      TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], ,
      drop = FALSE
    ]

    min_size <- ifelse(term %in% unlimited_db, 1, max(minGSSize, min.sz))
    max_size <- ifelse(term %in% unlimited_db, Inf, min(maxGSSize, max.sz))

    term_ids_vector <- as.character(TERM2GENE_tmp[["Term"]])
    gene_ids_vector <- as.character(TERM2GENE_tmp[[IDtype]])
    term_levels <- unique(term_ids_vector)
    term_index <- match(term_ids_vector, term_levels)
    term_sizes <- tabulate(term_index, nbins = length(term_levels))
    kept_terms <- term_levels[term_sizes >= min_size & term_sizes <= max_size]
    keep_initial <- term_ids_vector %in% kept_terms
    gene_sets <- split(
      gene_ids_vector[keep_initial],
      term_ids_vector[keep_initial]
    )

    if (length(gene_sets) == 0) {
      log_message(
        "No gene sets found for {.val {term}} after filtering",
        message_type = "warning",
        verbose = verbose
      )
      next
    }

    feature_names <- rownames(expr)
    feature_idx <- match(gene_ids_vector, feature_names)
    keep_overlap <- keep_initial & !is.na(feature_idx)
    overlap <- feature_names[sort(unique(feature_idx[keep_overlap]))]
    n_genes_in_sets <- length(unique(gene_ids_vector[keep_initial]))
    log_message(
      "Initial overlap: {.val {length(overlap)}} genes out of {.val {nrow(expr)}} expression genes and {.val {n_genes_in_sets}} genes in gene sets",
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
    gene_sets_filtered <- split(
      gene_ids_vector[keep_overlap],
      term_ids_vector[keep_overlap]
    )
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

    if (identical(backend, "cpp")) {
      if (identical(method, "gsva")) {
        gsva_scores <- t(run_gsva_scores(
          expr_counts = expr_filtered,
          gene_sets = gene_sets_filtered,
          kcdf = kcdf,
          min_gs_size = min_size,
          max_gs_size = max_size,
          max_diff = mx.diff,
          abs_ranking = abs.ranking,
          tau = tau,
          chunk_size = cpp_chunk_size
        ))
      } else if (identical(method, "ssgsea")) {
        gsva_scores <- t(run_ssgsea_scores(
          expr_counts = expr_filtered,
          gene_sets = gene_sets_filtered,
          min_gs_size = min_size,
          max_gs_size = max_size,
          alpha = tau,
          normalize = ssgsea.norm
        ))
      } else if (identical(method, "zscore")) {
        sparse_expr <- inherits(expr_filtered, "dgCMatrix")
        sparse_standardize_full <- sparse_expr &&
          gene_set_scoring_zscore_sparse_standardize_full()
        gsva_scores <- t(run_zscore_scores(
          expr_counts = expr_filtered,
          gene_sets = gene_sets_filtered,
          min_gs_size = min_size,
          max_gs_size = max_size,
          sparse_standardize = sparse_expr && !sparse_standardize_full,
          sparse_standardize_full = sparse_standardize_full
        ))
      } else {
        gsva_scores <- t(run_plage_scores(
          expr_counts = expr_filtered,
          gene_sets = gene_sets_filtered,
          min_gs_size = min_size,
          max_gs_size = max_size,
          dense_standardize = gene_set_scoring_plage_dense_standardize()
        ))
      }
      if (!is.matrix(gsva_scores)) {
        gsva_scores <- as_matrix(gsva_scores)
      }
      if (identical(method, "plage")) {
        gsva_scores <- orient_plage_scores(
          scores = gsva_scores,
          expr = expr_filtered,
          gene_sets = gene_sets_filtered
        )
      }
    } else {
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
      if (identical(method, "plage")) {
        gsva_scores <- orient_plage_scores(
          scores = gsva_scores,
          expr = expr_filtered,
          gene_sets = gene_sets_filtered
        )
      }
    }

    term_ids <- rownames(gsva_scores)
    term_name_lookup <- stats::setNames(
      as.character(TERM2NAME_tmp[["Name"]]),
      TERM2NAME_tmp[["Term"]]
    )
    term_names <- unname(term_name_lookup[term_ids])
    term_names[is.na(term_names)] <- term_ids[is.na(term_names)]
    term_names <- capitalize(
      trimws(as.character(term_names)),
      force_tolower = TRUE
    )
    term_gene_lookup <- gene_sets
    term_gene_ids <- vapply(
      term_ids,
      function(x) {
        genes <- term_gene_lookup[[x]]
        if (is.null(genes)) {
          genes <- gene_sets_filtered[[x]]
        }
        genes <- unique(as.character(genes[!is.na(genes)]))
        paste0(genes, collapse = "/")
      },
      FUN.VALUE = character(1)
    )
    rownames(gsva_scores) <- term_names

    gsva_results[[term]] <- gsva_scores
    n_terms <- length(term_ids)
    n_groups <- ncol(gsva_scores)
    enrichment_results[[term]] <- data.frame(
      ID = rep(term_ids, times = n_groups),
      Description = rep(term_names, times = n_groups),
      geneID = rep(term_gene_ids, times = n_groups),
      Groups = rep(colnames(gsva_scores), each = n_terms),
      Database = term,
      Version = as.character(db_list[[species]][[term]][["version"]]),
      GSVA_Score = as.numeric(gsva_scores),
      stringsAsFactors = FALSE
    )
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
    backend = backend,
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
    colnames(scores_mat) <- gene_set_scoring_make_score_colnames(
      feature_names = rownames(gsva_scores_combined),
      prefix = assay_name,
      fallback = "GSVA"
    )
    if (isTRUE(new_assay)) {
      srt[[assay_name]] <- Seurat::CreateAssayObject(
        data = Matrix::t(as_matrix(scores_mat))
      )
      srt[[assay_name]] <- Seurat::AddMetaData(
        object = srt[[assay_name]],
        metadata = data.frame(
          termnames = rownames(gsva_scores_combined),
          row.names = rownames(gsva_scores_combined)
        )
      )
    }
    if (isTRUE(store_metadata)) {
      srt <- Seurat::AddMetaData(object = srt, metadata = as.data.frame(scores_mat))
    }
    if (isTRUE(new_assay) && isTRUE(store_metadata)) {
      log_message(
        "{.pkg GSVA} results stored in assay {.val {assay_name}}, meta.data, and tools slot {.val {tool_name}}",
        verbose = verbose
      )
    } else if (isTRUE(new_assay)) {
      log_message(
        "{.pkg GSVA} results stored in assay {.val {assay_name}} and tools slot {.val {tool_name}}",
        verbose = verbose
      )
    } else if (isTRUE(store_metadata)) {
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
