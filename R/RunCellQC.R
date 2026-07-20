#' @title Run doublet-calling for single cell RNA-seq data.
#'
#' @md
#' @inheritParams standard_scop
#' @param assay The name of the assay to be used for doublet-calling.
#' Default is `"RNA"`.
#' @param db_method Method used for doublet-calling.
#' Can be one of `"scDblFinder"`, `"Scrublet"`, `"DoubletDetection"`,
#' `"scds_cxds"`, `"scds_bcds"`, `"scds_hybrid"`.
#' @param db_rate The expected doublet rate.
#' Default is calculated as `ncol(srt) / 1000 * 0.01`.
#' @param data_type Optional precomputed result from [CheckDataType] for the
#' input assay. Primarily used internally to avoid repeated scans of the same
#' count matrix across nested QC calls.
#' @param ... Additional arguments to be passed to the corresponding doublet-calling method.
#'
#' @return Returns a Seurat object with the doublet prediction results and prediction scores stored in the meta.data.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDoubletCalling(
#'   pancreas_sub,
#'   db_method = "scDblFinder"
#' )
#' table(pancreas_sub$db.scDblFinder_class)
#' head(pancreas_sub@meta.data[, c(
#'   "db.scDblFinder_class", "db.scDblFinder_score"
#' )])
RunDoubletCalling <- function(
  srt,
  assay = "RNA",
  db_rate = ncol(srt) / 1000 * 0.01,
  db_method = "scDblFinder",
  data_type = NULL,
  ...,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- data_type %||% CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  db_methods <- c(
    "scDblFinder",
    "Scrublet",
    "DoubletDetection",
    "scds_cxds",
    "scds_bcds",
    "scds_hybrid"
  )
  if (db_method %in% db_methods) {
    methods <- unlist(strsplit(db_method, "_"))
    method1 <- methods[1]
    method2 <- methods[2]
    if (is.na(method2)) {
      args1 <- mget(names(formals()), sys.frame(sys.nframe()))
      args2 <- as.list(match.call())
    } else {
      args1 <- c(
        mget(names(formals()), sys.frame(sys.nframe())),
        method = method2
      )
      args2 <- c(as.list(match.call()), method = method2)
    }
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("db_method", "...")]
    tryCatch(
      expr = {
        srt <- do.call(
          what = paste0("db_", method1),
          args = args1
        )
      },
      error = function(e) {
        err_msg <- conditionMessage(e)
        err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
        err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
        log_message(err_msg, message_type = "error")
      }
    )
    return(srt)
  } else {
    log_message(
      "{.arg db_method} must be one of {.val {db_methods}}",
      message_type = "error"
    )
  }
}

#' @title Run doublet-calling with scDblFinder
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param ... Additional arguments to be passed to [scDblFinder::scDblFinder()].
#'
#' @export
db_scDblFinder <- function(
  srt,
  assay = "RNA",
  db_rate = ncol(srt) / 1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE
) {
  log_message("Running {.pkg scDblFinder}", verbose = verbose)
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- data_type %||% CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  check_r("scDblFinder", verbose = FALSE)
  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  sce <- scDblFinder::scDblFinder(sce, dbr = db_rate, verbose = FALSE, ...)
  srt[["db.scDblFinder_score"]] <- sce[["scDblFinder.score"]]
  srt[["db.scDblFinder_class"]] <- sce[["scDblFinder.class"]]
  return(srt)
}

#' @title Run doublet-calling with scds
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param method The method to be used for doublet-calling.
#' Options are `"hybrid"`, `"cxds"`, or `"bcds"`.
#' @param ... Additional arguments passed to the selected `scds` method.
#'
#' @export
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.scds_hybrid_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.scds_hybrid_score"
#' )
db_scds <- function(
  srt,
  assay = "RNA",
  db_rate = ncol(srt) / 1000 * 0.01,
  method = c("hybrid", "cxds", "bcds"),
  data_type = NULL,
  ...,
  verbose = TRUE
) {
  log_message("Running {.pkg scds} with method {.val {method}}", verbose = verbose)
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- data_type %||% CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  check_r("scds", verbose = FALSE)
  method <- match.arg(method)
  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  sce <- switch(
    method,
    cxds = scds::cxds(sce, ...),
    bcds = scds::bcds(sce, ...),
    hybrid = scds::cxds_bcds_hybrid(sce, ...)
  )
  score_col <- paste0(method, "_score")
  srt[[paste0("db.scds_", method, "_score")]] <- sce[[score_col]]
  ntop <- ceiling(db_rate * ncol(sce))
  db_qc <- names(sort(
    srt[[paste0("db.scds_", method, "_score"), drop = TRUE]],
    decreasing = TRUE
  )[1:ntop])
  srt[[paste0("db.scds_", method, "_class")]] <- "singlet"
  srt[[paste0("db.scds_", method, "_class")]][db_qc, ] <- "doublet"
  return(srt)
}

#' @title Run doublet-calling with Scrublet
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param ... Additional arguments to be passed to [scrublet.Scrublet](https://github.com/swolock/scrublet).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- db_Scrublet(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.Scrublet_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.Scrublet_score"
#' )
#' }
db_Scrublet <- function(
  srt,
  assay = "RNA",
  db_rate = ncol(srt) / 1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE
) {
  user_args <- list(...)
  log_message("Running {.pkg Scrublet}", verbose = verbose)
  PrepareEnv(modules = "scrublet")
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- data_type %||% CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  counts <- GetAssayData5(
    object = srt,
    assay = assay,
    layer = "counts"
  )
  raw_counts <- Matrix::t(as_matrix(counts))
  check_python("scrublet", verbose = FALSE)
  scr <- reticulate::import("scrublet")
  scrub_doublets_args <- c(
    "synthetic_doublet_umi_subsampling",
    "use_approx_neighbors",
    "distance_metric",
    "get_doublet_neighbor_parents",
    "min_counts",
    "min_cells",
    "min_gene_variability_pctl",
    "log_transform",
    "mean_center",
    "normalize_variance",
    "n_prin_comps",
    "svd_solver"
  )
  score_args <- user_args[intersect(names(user_args), scrub_doublets_args)]
  integer_score_args <- c(
    "min_counts",
    "min_cells",
    "min_gene_variability_pctl",
    "n_prin_comps"
  )
  for (arg in intersect(names(score_args), integer_score_args)) {
    score_args[[arg]] <- as.integer(score_args[[arg]])
  }
  constructor_args <- user_args[setdiff(names(user_args), names(score_args))]
  scrub <- do.call(
    scr$Scrublet,
    c(list(raw_counts, expected_doublet_rate = db_rate), constructor_args)
  )
  res <- do.call(scrub$scrub_doublets, score_args)
  doublet_scores <- res[[1]]
  predicted_doublets <- res[[2]]
  scrublet_threshold <- NA_real_
  if (reticulate::py_has_attr(scrub, "threshold_")) {
    scrublet_threshold <- suppressWarnings(
      as.numeric(reticulate::py_to_r(scrub$threshold_))
    )[1]
  }
  parameters <- list()
  detected_doublet_rate <- mean(predicted_doublets)

  srt[["db.Scrublet_score"]] <- doublet_scores
  srt[["db.Scrublet_threshold"]] <- rep(scrublet_threshold, ncol(srt))
  srt[["db.Scrublet_detected_rate"]] <- rep(detected_doublet_rate, ncol(srt))
  srt[["db.Scrublet_class"]] <- sapply(
    predicted_doublets,
    function(i) {
      switch(as.character(i),
        "FALSE" = "singlet",
        "TRUE" = "doublet"
      )
    }
  )
  srt@tools[["Scrublet"]] <- list(
    expected_doublet_rate = db_rate,
    threshold = scrublet_threshold,
    detected_doublet_rate = detected_doublet_rate,
    parameters = parameters
  )
  log_message(
    "{.pkg Scrublet} doublet calling completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt)
}

#' @title Run doublet-calling with DoubletDetection
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param cores The number of CPU cores to use for `doubletdetection`.
#' Default is `1`.
#' @param ... Additional arguments to be passed to [doubletdetection.BoostClassifier](https://github.com/JonathanShor/DoubletDetection).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- db_DoubletDetection(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.DoubletDetection_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.DoubletDetection_score"
#' )
#' }
db_DoubletDetection <- function(
  srt,
  assay = "RNA",
  db_rate = ncol(srt) / 1000 * 0.01,
  cores = 1,
  data_type = NULL,
  ...,
  verbose = TRUE
) {
  log_message("Running {.pkg DoubletDetection}", verbose = verbose)
  # Sys.setenv(NUMBA_NUM_THREADS = "1")
  # Sys.setenv(NUMBA_DISABLE_JIT = "0")
  PrepareEnv(modules = "doubletdetection")

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- data_type %||% CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  user_args <- list(...)

  if (!"n_jobs" %in% names(user_args)) {
    user_args[["n_jobs"]] <- as.integer(cores)
  }
  if (!"n_iters" %in% names(user_args)) {
    user_args[["n_iters"]] <- as.integer(5)
  }
  if (!"standard_scaling" %in% names(user_args)) {
    user_args[["standard_scaling"]] <- TRUE
  }

  check_python("doubletdetection")
  check_python("louvain")

  doubletdetection <- reticulate::import("doubletdetection")
  counts <- GetAssayData5(
    object = srt,
    assay = assay,
    layer = "counts"
  )
  clf <- do.call(doubletdetection$BoostClassifier, user_args)

  clf_fit <- clf$fit(Matrix::t(counts))

  labels <- clf_fit$predict()
  scores <- clf_fit$doublet_score()

  labels <- as.integer(reticulate::py_to_r(labels))
  scores <- as.numeric(reticulate::py_to_r(scores))

  n_cells <- ncol(srt)
  if (length(labels) != n_cells) {
    log_message(
      "Length of labels ({.val {length(labels)}}) does not match number of cells",
      message_type = "error"
    )
  }
  if (length(scores) != n_cells) {
    log_message(
      "Length of scores ({.val {length(scores)}}) does not match number of cells",
      message_type = "error"
    )
  }

  cell_names <- colnames(srt)

  if (!is.null(names(scores)) && !identical(names(scores), cell_names)) {
    scores <- scores[cell_names]
  }
  if (!is.null(names(labels)) && !identical(names(labels), cell_names)) {
    labels <- labels[cell_names]
  }

  scores <- unname(scores)
  labels <- unname(labels)

  class_map <- c(
    "0" = "singlet",
    "1" = "doublet"
  )
  class_labels <- unname(class_map[as.character(labels)])
  idx_unknown <- is.na(class_labels)
  if (any(idx_unknown)) {
    log_message(
      "Found {.val {sum(idx_unknown)}} cells with unexpected DoubletDetection labels; setting to {.val doublet}",
      message_type = "warning",
      verbose = verbose
    )
    class_labels[idx_unknown] <- "doublet"
  }
  class_labels <- factor(class_labels, levels = c("singlet", "doublet"))

  srt@meta.data$db.DoubletDetection_score <- scores
  srt@meta.data$db.DoubletDetection_class <- class_labels

  log_message(
    "{.pkg DoubletDetection} doublet calling completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt)
}

#' @title Run ambient RNA decontamination with decontX
#'
#' @md
#' @inheritParams standard_scop
#' @param assay The name of the assay to be used for decontamination.
#' Default is `"RNA"`.
#' @param group.by Cell cluster labels passed to [decontX::decontX()].
#' Can be `NULL`, a meta.data column name, or a vector aligned to cells.
#' Default is `NULL`.
#' @param batch Batch labels passed to [decontX::decontX()].
#' Can be `NULL`, a meta.data column name, or a vector aligned to cells.
#' Default is `NULL`.
#' @param background Optional background / empty-droplet input passed to [decontX::decontX()].
#' Can be a `Seurat` object, `SingleCellExperiment`, or count matrix.
#' Default is `NULL`.
#' @param background_assay Assay name used when `background` is a `Seurat` object
#' or `SingleCellExperiment`.
#' Default is `NULL`, which falls back to `assay` for `Seurat` background
#' and `"counts"` for `SingleCellExperiment` background.
#' @param bg_batch Batch labels for `background` passed to [decontX::decontX()].
#' Can be `NULL`, a metadata column name, or a vector aligned to the background droplets.
#' Default is `NULL`.
#' @param assay_name Name of the assay used to store decontaminated counts.
#' Default is `"decontXcounts"`.
#' @param store_assay Whether to store decontaminated counts as a new assay.
#' Default is `TRUE`.
#' @param round_counts Whether to round decontaminated counts before creating the assay.
#' Default is `FALSE`.
#' @param data_type Optional precomputed result from [CheckDataType()] for the
#' input assay. Primarily used internally to avoid repeated scans of the same
#' count matrix across nested QC calls.
#' @param ... Additional arguments passed to [decontX::decontX()].
#'
#' @return Returns a Seurat object with decontX contamination estimates stored in the meta.data,
#' and optional decontaminated counts stored in a new assay.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDecontX(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#'
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = "decontX_contamination"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "decontX_contamination"
#' )
RunDecontX <- function(
  srt,
  assay = "RNA",
  group.by = NULL,
  batch = NULL,
  background = NULL,
  background_assay = NULL,
  bg_batch = NULL,
  assay_name = "decontXcounts",
  store_assay = TRUE,
  round_counts = FALSE,
  data_type = NULL,
  seed = 11,
  ...,
  verbose = TRUE
) {
  log_message("Running {.pkg decontX}", verbose = verbose)

  .decontx_get_meta <- function(object) {
    if (inherits(object, "Seurat")) {
      return(object@meta.data)
    }
    if (
      inherits(object, "SingleCellExperiment") ||
        inherits(object, "SummarizedExperiment")
    ) {
      return(as.data.frame(SummarizedExperiment::colData(object)))
    }
    NULL
  }

  .decontx_resolve_vector <- function(object, value, arg_name) {
    if (is.null(value)) {
      return(NULL)
    }
    meta <- .decontx_get_meta(object)
    cells <- colnames(object)

    if (
      !is.null(meta) &&
        is.character(value) &&
        length(value) == 1 &&
        value %in% colnames(meta)
    ) {
      return(meta[[value]])
    }

    if (!is.null(names(value))) {
      idx <- match(cells, names(value))
      if (all(!is.na(idx))) {
        return(unname(value[idx]))
      }
    }

    if (length(value) == length(cells)) {
      return(unname(value))
    }

    log_message(
      "{.arg {arg_name}} must be NULL, a metadata column name, or a vector aligned to the cells",
      message_type = "error"
    )
  }

  .decontx_as_input <- function(
    object,
    assay_name = NULL,
    arg_name = "background"
  ) {
    if (is.null(object)) {
      return(NULL)
    }
    if (inherits(object, "Seurat")) {
      assay_use <- assay_name %||% assay
      if (isFALSE(assay_use %in% SeuratObject::Assays(object))) {
        log_message(
          "{.arg {arg_name}} does not contain assay {.val {assay_use}}",
          message_type = "error"
        )
      }
      return(Seurat::as.SingleCellExperiment(object, assay = assay_use))
    }
    if (inherits(object, "SingleCellExperiment")) {
      return(object)
    }
    if (inherits(object, "Matrix") || is.matrix(object)) {
      return(object)
    }
    log_message(
      "{.arg {arg_name}} must be a {.cls Seurat}, {.cls SingleCellExperiment}, or count matrix",
      message_type = "error"
    )
  }

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (isFALSE(assay %in% SeuratObject::Assays(srt))) {
    log_message(
      "{.arg srt} does not contain {.arg {assay}} assay",
      message_type = "error"
    )
  }

  status <- data_type %||% CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "warning",
      verbose = verbose
    )
  }

  check_r(c("decontX", "SingleCellExperiment"), verbose = FALSE)

  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  group_use <- .decontx_resolve_vector(
    srt,
    group.by,
    "group.by"
  )
  if (!is.null(group_use) && !is.factor(group_use)) {
    group_use <- factor(group_use)
  }
  batch_use <- .decontx_resolve_vector(srt, batch, "batch")

  background_is_seurat <- inherits(background, "Seurat")
  background_is_sce <- inherits(background, "SingleCellExperiment")
  background_input <- .decontx_as_input(
    object = background,
    assay_name = background_assay,
    arg_name = "background"
  )
  bg_assay_name_use <- NULL
  if (background_is_seurat) {
    bg_assay_name_use <- "counts"
  } else if (background_is_sce) {
    bg_assay_name_use <- background_assay %||% "counts"
  }
  if (is.null(background_input) && !is.null(bg_batch)) {
    log_message(
      "{.arg bg_batch} requires {.arg background} to be provided.",
      message_type = "error"
    )
  }
  bg_batch_use <- .decontx_resolve_vector(
    background_input,
    bg_batch,
    "bg_batch"
  )

  user_args <- list(...)
  decontx_args <- c(
    list(
      x = sce,
      assayName = "counts",
      z = group_use,
      batch = batch_use,
      background = background_input,
      bgBatch = bg_batch_use,
      seed = seed
    ),
    user_args
  )
  if (is.null(decontx_args[["verbose"]])) {
    decontx_args[["verbose"]] <- FALSE
  }
  if (inherits(background_input, "SingleCellExperiment")) {
    decontx_args[["bgAssayName"]] <- bg_assay_name_use %||% "counts"
  }

  sce <- do.call(decontX::decontX, args = decontx_args)

  metadata_to_add <- data.frame(
    decontX_contamination = sce[["decontX_contamination"]],
    decontX_clusters = sce[["decontX_clusters"]],
    row.names = colnames(srt)
  )
  srt <- Seurat::AddMetaData(
    object = srt,
    metadata = metadata_to_add
  )
  srt@tools[["decontX"]] <- S4Vectors::metadata(sce)[["decontX"]]

  decontX_contamination_values <- sce[["decontX_contamination"]]
  if (
    length(decontX_contamination_values) > 0 &&
      any(!is.na(decontX_contamination_values))
  ) {
    decontX_contamination_summary <- sprintf(
      "%.4f / %.4f / %.4f",
      stats::median(decontX_contamination_values, na.rm = TRUE),
      mean(decontX_contamination_values, na.rm = TRUE),
      max(decontX_contamination_values, na.rm = TRUE)
    )
    log_message(
      "decontX contamination (median/mean/max): {decontX_contamination_summary}",
      verbose = verbose
    )
  }

  if (isTRUE(store_assay)) {
    decontx_counts <- decontX::decontXcounts(sce)
    if (isTRUE(round_counts)) {
      decontx_counts <- round(decontx_counts)
    }
    if (!inherits(decontx_counts, "Matrix")) {
      decontx_counts <- Matrix::Matrix(decontx_counts, sparse = TRUE)
    }
    srt[[assay_name]] <- Seurat::CreateAssayObject(counts = decontx_counts)
    log_message("decontX assay stored as {assay_name}", verbose = verbose)
  }

  log_message(
    "{.pkg decontX} decontamination completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt)
}

#' @title Run cell-level quality control
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams RunDoubletCalling
#' @inheritParams standard_scop
#' @param split.by Name of a meta.data column used to split the object before QC.
#' Default is `NULL`. When specified, QC and doublet-calling are performed
#' separately within each split object and merged back afterward.
#' @param return_filtered Logical indicating whether to return a cell-filtered Seurat object.
#' Default is `FALSE`.
#' @param qc_metrics A character vector specifying the quality control metrics to be applied.
#' Available metrics are `"doublets"`, `"decontX"`, `"atac"`, `"outlier"`, `"umi"`, `"gene"`,
#' `"mito"`, `"ribo"`, `"hb"`, `"ribo_mito_ratio"`, `"species"`, and any
#' rule name supplied in `qc_features`.
#' Default is `c("doublets", "decontX", "outlier", "umi", "gene", "mito", "ribo", "ribo_mito_ratio", "species")`.
#' For `ChromatinAssay`, if `.arg qc_metrics` is not supplied, the default is `"atac"`.
#' @param db_method Method used for doublet-calling.
#' Can be one of `"scDblFinder"`, `"Scrublet"`, `"DoubletDetection"`,
#' `"scds_cxds"`, `"scds_bcds"`, `"scds_hybrid"`.
#' The resulting doublet labels are aggregated afterward into `db_qc` and do not
#' affect the thresholds used by the other QC metrics.
#' @param outlier_threshold A character vector specifying the outlier threshold.
#' Default is `c("log10_nCount:lower:2.5", "log10_nCount:higher:5", "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5")`.
#' @param db_coefficient The coefficient used to calculate the doublet rate.
#' Default is `0.01`. Doublet rate is calculated as `ncol(srt) / 1000 * db_coefficient`.
#' @param decontX_threshold Optional contamination threshold used to filter cells
#' after running [RunDecontX()]. Cells with `decontX_contamination` greater than
#' this value are marked as failed in `decontX_qc`. Default is `NULL`, which
#' computes decontX results without filtering cells by contamination.
#' @param group.by Group labels passed to [RunDecontX()] when `"decontX"`
#' is included in `qc_metrics`. Can be `NULL`, a meta.data column name, or a vector
#' aligned to cells. Default is `NULL`.
#' @param decontX_batch Batch labels passed to [RunDecontX()] when `"decontX"`
#' is included in `qc_metrics`. Default is `NULL`.
#' @param decontX_background Optional background / empty-droplet input passed to
#' [RunDecontX()] when `"decontX"` is included in `qc_metrics`.
#' Default is `NULL`.
#' @param decontX_background_assay Assay name used when `decontX_background` is a
#' `Seurat` object or `SingleCellExperiment`. Default is `NULL`.
#' @param decontX_bg_batch Batch labels for `decontX_background` passed to
#' [RunDecontX()]. Default is `NULL`.
#' @param decontX_assay_name Name of the assay used to store decontaminated counts
#' from [RunDecontX()]. Default is `"decontXcounts"`.
#' @param decontX_store_assay Whether to store decontaminated counts as a new assay
#' when running [RunDecontX()]. Default is `FALSE`.
#' @param decontX_round_counts Whether to round decontaminated counts before creating
#' the assay in [RunDecontX()]. Default is `TRUE`.
#' @param decontX_args A named list of additional advanced arguments passed to
#' [RunDecontX()] when `"decontX"` is included in `qc_metrics`.
#' Explicit `decontX_*` parameters are preferred for common options and take
#' precedence when both are supplied.
#' Default is `list()`.
#' @param atac_args A named list of additional arguments passed to [RunATACQC()]
#' when `"atac"` is included in `qc_metrics`. Threshold arguments from
#' [RunATACQC()] are used to label failed cells in `atac_qc`, but filtering is
#' deferred to [RunCellQC()]. Default is `list()`.
#' @param outlier_n Minimum number of outlier metrics that meet the conditions for determining outlier cells.
#' Default is `1`.
#' @param UMI_threshold UMI number threshold.
#' Cells that exceed this threshold will be considered as kept.
#' Default is `3000`.
#' @param gene_threshold Gene number threshold.
#' Cells that exceed this threshold will be considered as kept.
#' Default is `1000`.
#' @param mito_threshold Percentage of UMI counts of mitochondrial genes.
#' Cells that exceed this threshold will be considered as discarded.
#' Default is `20`.
#' @param mito_pattern Regex patterns to match the mitochondrial genes.
#' Default is `c("MT-", "Mt-", "mt-")`.
#' @param mito_gene A defined mitochondrial genes.
#' If features provided, will ignore the `mito_pattern` matching.
#' Default is `NULL`.
#' @param ribo_threshold Percentage of UMI counts of ribosomal genes.
#' Cells that exceed this threshold will be considered as discarded.
#' Default is `50`.
#' @param ribo_pattern Regex patterns to match the ribosomal genes.
#' Default is `c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$", "rp[sl]\\d+\\w{0,1}\\d*$")`.
#' @param ribo_gene A defined ribosomal genes. If features provided, will ignore the `ribo_pattern` matching.
#' Default is `NULL`.
#' @param ribo_mito_ratio_range A numeric vector specifying the range of ribosomal/mitochondrial gene expression ratios for ribo_mito_ratio outlier cells.
#' Default is `c(1, Inf)`.
#' @param species Species used as the suffix of the QC metrics.
#' The first is the species of interest.
#' Default is `NULL`.
#' @param species_gene_prefix Species gene prefix used to calculate QC metrics for each species.
#' Default is `NULL`.
#' @param species_percent Percentage of UMI counts of the first species.
#' Cells that exceed this threshold will be considered as kept.
#' Default is `95`.
#' @param hb_range A numeric vector specifying the accepted percentage range for
#' hemoglobin features. Values outside the range fail `hb_qc` when `"hb"` is
#' included in `qc_metrics`. Default is `c(0, 5)`.
#' @param hb_pattern Regex patterns used to match hemoglobin features. The
#' defaults match common human, mouse, and lower-case symbols while excluding
#' `HBP`/`hbp` genes. Default is `c("HB[^P]", "Hb[^p]", "hb[^p]")`.
#' @param hb_gene A defined set of hemoglobin features. When supplied, these
#' features are used instead of `hb_pattern`. Default is `NULL`.
#' @param qc_features A named list defining additional feature-percentage QC
#' rules. Each rule must contain exactly one of `features` (a character vector)
#' or `pattern` (one or more regex patterns), plus a numeric `range` of length
#' two. Providing a rule computes `percent.<name>`; the rule filters cells and
#' creates `<name>_qc` only when `<name>` is included in `qc_metrics`.
#' Patterns are matched directly against assay feature names and do not receive
#' species prefixes. Rule names cannot collide with built-in QC columns.
#' Default is `list()`.
#'
#' @return Returns Seurat object with the QC results stored in the meta.data layer.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCellQC(
#'   pancreas_sub,
#'   qc_metrics = c("umi", "gene", "example_set"),
#'   qc_features = list(
#'     example_set = list(
#'       features = head(rownames(pancreas_sub), 5),
#'       range = c(0, 0)
#'     )
#'   )
#' )
#' head(pancreas_sub@meta.data[, c("percent.hb", "percent.example_set")])
#' # Hemoglobin expression can be biological signal in erythroid samples;
#' # include "hb" in qc_metrics only when HB-based filtering is appropriate.
#'
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "umi_qc", "gene_qc", "example_set_qc"
#'   ),
#'   plot_type = "upset",
#'   stat_level = "Fail"
#' )
RunCellQC <- function(
  srt,
  assay = "RNA",
  split.by = NULL,
  group.by = NULL,
  return_filtered = FALSE,
  qc_metrics = c(
    "doublets",
    "decontX",
    "atac",
    "outlier",
    "umi",
    "gene",
    "mito",
    "ribo",
    "ribo_mito_ratio",
    "species"
  ),
  db_method = "scDblFinder",
  db_rate = NULL,
  db_coefficient = 0.01,
  decontX_threshold = NULL,
  decontX_batch = NULL,
  decontX_background = NULL,
  decontX_background_assay = NULL,
  decontX_bg_batch = NULL,
  decontX_assay_name = "decontXcounts",
  decontX_store_assay = FALSE,
  decontX_round_counts = TRUE,
  decontX_args = list(),
  atac_args = list(),
  outlier_threshold = c(
    "log10_nCount:lower:2.5",
    "log10_nCount:higher:5",
    "log10_nFeature:lower:2.5",
    "log10_nFeature:higher:5",
    "featurecount_dist:lower:2.5"
  ),
  outlier_n = 1,
  UMI_threshold = 3000,
  gene_threshold = 1000,
  mito_threshold = 20,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  ribo_threshold = 50,
  ribo_pattern = c(
    "RP[SL]\\d+\\w{0,1}\\d*$",
    "Rp[sl]\\d+\\w{0,1}\\d*$",
    "rp[sl]\\d+\\w{0,1}\\d*$"
  ),
  ribo_gene = NULL,
  ribo_mito_ratio_range = c(1, Inf),
  species = NULL,
  species_gene_prefix = NULL,
  species_percent = 95,
  seed = 11,
  verbose = TRUE,
  hb_range = c(0, 5),
  hb_pattern = c("HB[^P]", "Hb[^p]", "hb[^p]"),
  hb_gene = NULL,
  qc_features = list()
) {
  log_message(
    "Running cell-level quality control",
    message_type = "running",
    verbose = verbose
  )
  set.seed(seed)

  group.by_missing <- missing(group.by)
  decontX_batch_missing <- missing(decontX_batch)
  decontX_background_missing <- missing(decontX_background)
  decontX_background_assay_missing <- missing(decontX_background_assay)
  decontX_bg_batch_missing <- missing(decontX_bg_batch)
  decontX_assay_name_missing <- missing(decontX_assay_name)
  decontX_store_assay_missing <- missing(decontX_store_assay)
  decontX_round_counts_missing <- missing(decontX_round_counts)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (isFALSE(assay %in% SeuratObject::Assays(srt))) {
    log_message(
      "{.arg srt} does not contain {.arg {assay}} assay",
      message_type = "error"
    )
  }
  if (length(species) != length(species_gene_prefix)) {
    log_message(
      "{.arg species_gene_prefix} must be the same length as {.arg species}.",
      message_type = "error"
    )
  }
  if (length(species) == 0) {
    species <- species_gene_prefix <- NULL
  }
  if (
    !is.null(decontX_threshold) &&
      (!is.numeric(decontX_threshold) || length(decontX_threshold) != 1)
  ) {
    log_message(
      "{.arg decontX_threshold} must be NULL or a single numeric value.",
      message_type = "error"
    )
  }
  if (!is.list(decontX_args)) {
    log_message(
      "{.arg decontX_args} must be a named list.",
      message_type = "error"
    )
  }
  if (!is.list(atac_args)) {
    log_message(
      "{.arg atac_args} must be a named list.",
      message_type = "error"
    )
  }
  is_atac_assay <- inherits(srt[[assay]], "ChromatinAssay")
  if (isTRUE(is_atac_assay) && missing(qc_metrics)) {
    qc_metrics <- "atac"
  }
  if (!is.character(qc_metrics) || anyNA(qc_metrics)) {
    log_message(
      "{.arg qc_metrics} must be a character vector without missing values",
      message_type = "error"
    )
  }
  if (
    isTRUE(is_atac_assay) &&
      ("hb" %in% qc_metrics || !is.null(hb_gene) || length(qc_features) > 0)
  ) {
    log_message(
      "{.arg hb} and {.arg qc_features} require an RNA-like assay, not a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  assay_features <- rownames(Seurat::GetAssay(srt, assay = assay))
  if (!is.list(qc_features)) {
    stop("qc_features must be a named list", call. = FALSE)
  }
  if (length(qc_features) > 0) {
    rule_names <- names(qc_features)
    if (
      is.null(rule_names) ||
        anyNA(rule_names) ||
        any(!nzchar(rule_names)) ||
        anyDuplicated(rule_names)
    ) {
      stop("qc_features must have unique, non-empty names", call. = FALSE)
    }
    if (any(make.names(rule_names) != rule_names)) {
      stop("qc_features names must be syntactically valid R names", call. = FALSE)
    }
    species_suffix <- c("", paste0(".", species))
    built_in_columns <- c(
      unlist(lapply(species_suffix, function(suffix) {
        paste0(
          c("percent.mito", "percent.ribo", "percent.hb", "percent.genome"),
          suffix
        )
      }), use.names = FALSE),
      "db_qc", "atac_qc", "decontX_qc", "outlier_qc", "umi_qc",
      "gene_qc", "mito_qc", "ribo_qc", "hb_qc", "ribo_mito_ratio_qc",
      "species_qc", "CellQC"
    )
    generated_columns <- c(
      paste0("percent.", rule_names),
      paste0(rule_names, "_qc")
    )
    conflicting_columns <- intersect(generated_columns, built_in_columns)
    if (length(conflicting_columns) > 0) {
      stop(
        sprintf(
          "qc_features names generate columns that conflict with built-in QC columns: %s",
          paste(conflicting_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    qc_features_validated <- vector("list", length(qc_features))
    names(qc_features_validated) <- rule_names
    for (rule_name in rule_names) {
      rule <- qc_features[[rule_name]]
      if (!is.list(rule)) {
        stop(
          sprintf("qc_features rule '%s' must be a list", rule_name),
          call. = FALSE
        )
      }
      unknown_fields <- setdiff(names(rule), c("features", "pattern", "range"))
      if (length(unknown_fields) > 0) {
        stop(
          sprintf(
            "qc_features rule '%s' contains unsupported fields: %s",
            rule_name,
            paste(unknown_fields, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      has_features <- !is.null(rule[["features"]])
      has_pattern <- !is.null(rule[["pattern"]])
      if (has_features == has_pattern) {
        stop(
          sprintf(
            "qc_features rule '%s' must provide exactly one of features or pattern",
            rule_name
          ),
          call. = FALSE
        )
      }
      rule_range <- cellqc_validate_range(
        rule[["range"]],
        sprintf("qc_features[['%s']]$range", rule_name)
      )
      if (has_features) {
        requested <- rule[["features"]]
        if (
          !is.character(requested) ||
            length(requested) == 0 ||
            anyNA(requested) ||
            any(!nzchar(requested))
        ) {
          stop(
            sprintf(
              "qc_features rule '%s' features must be a non-empty character vector",
              rule_name
            ),
            call. = FALSE
          )
        }
        requested <- unique(requested)
        resolved <- intersect(requested, assay_features)
        missing_features <- setdiff(requested, resolved)
        if (length(missing_features) > 0 && length(resolved) > 0) {
          warning(
            sprintf(
              "qc_features rule '%s' ignored missing features: %s",
              rule_name,
              paste(missing_features, collapse = ", ")
            ),
            call. = FALSE
          )
        }
      } else {
        patterns <- rule[["pattern"]]
        if (
          !is.character(patterns) ||
            length(patterns) == 0 ||
            anyNA(patterns) ||
            any(!nzchar(patterns))
        ) {
          stop(
            sprintf(
              "qc_features rule '%s' pattern must be a non-empty character vector",
              rule_name
            ),
            call. = FALSE
          )
        }
        resolved <- tryCatch(
          unique(unlist(lapply(
            patterns,
            function(pattern) grep(pattern, assay_features, value = TRUE)
          ))),
          error = function(e) {
            stop(
              sprintf(
                "qc_features rule '%s' has an invalid pattern: %s",
                rule_name,
                conditionMessage(e)
              ),
              call. = FALSE
            )
          }
        )
      }
      if (length(resolved) == 0) {
        stop(
          sprintf(
            "qc_features rule '%s' did not match any assay features",
            rule_name
          ),
          call. = FALSE
        )
      }
      qc_features_validated[[rule_name]] <- list(
        features = resolved,
        range = rule_range
      )
    }
    qc_features <- qc_features_validated
  }
  hb_range <- cellqc_validate_range(hb_range, "hb_range")
  if (
    !is.null(hb_gene) &&
      (
        !is.character(hb_gene) ||
          length(hb_gene) == 0 ||
          anyNA(hb_gene) ||
          any(!nzchar(hb_gene))
      )
  ) {
    log_message(
      "{.arg hb_gene} must be NULL or a non-empty character vector",
      message_type = "error"
    )
  }
  if (is.null(hb_gene)) {
    if (
      !is.character(hb_pattern) ||
        length(hb_pattern) == 0 ||
        anyNA(hb_pattern) ||
        any(!nzchar(hb_pattern))
    ) {
      log_message(
        "{.arg hb_pattern} must be a non-empty character vector",
        message_type = "error"
      )
    }
    tryCatch(
      lapply(hb_pattern, grep, x = assay_features),
      error = function(e) {
        stop(
          sprintf("hb_pattern contains an invalid regex: %s", conditionMessage(e)),
          call. = FALSE
        )
      }
    )
  }
  if (!is.null(hb_gene)) {
    hb_gene <- unique(hb_gene)
    hb_gene_found <- intersect(hb_gene, assay_features)
    hb_gene_missing <- setdiff(hb_gene, hb_gene_found)
    if (length(hb_gene_found) == 0) {
      log_message(
        "{.arg hb_gene} did not match any assay features",
        message_type = "error"
      )
    }
    if (length(hb_gene_missing) > 0) {
      warning(
        sprintf(
          "hb_gene ignored missing features: %s",
          paste(hb_gene_missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    hb_gene <- hb_gene_found
  }
  qc_metrics_available <- c(
    "doublets", "decontX", "atac", "outlier", "umi", "gene",
    "mito", "ribo", "hb", "ribo_mito_ratio", "species",
    names(qc_features)
  )
  if (any(!qc_metrics %in% qc_metrics_available)) {
    log_message(
      "{.arg qc_metrics} contains unsupported values: {.val {setdiff(qc_metrics, qc_metrics_available)}}",
      message_type = "error"
    )
  }
  active_qc_features <- intersect(names(qc_features), qc_metrics)
  decontX_assay_name_use <- decontX_args[["assay_name"]] %||% "decontXcounts"
  if (!decontX_assay_name_missing) {
    decontX_assay_name_use <- decontX_assay_name
  }
  decontX_store_assay_use <- decontX_args[["store_assay"]] %||% TRUE
  if (!decontX_store_assay_missing) {
    decontX_store_assay_use <- decontX_store_assay
  }
  decontX_round_counts_use <- decontX_args[["round_counts"]] %||% FALSE
  if (!decontX_round_counts_missing) {
    decontX_round_counts_use <- decontX_round_counts
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "warning",
      verbose = verbose
    )
  }
  srt <- cellqc_initialize_count_metrics(srt, assay = assay)
  srt_raw <- srt
  if (!is.null(split.by)) {
    split_cells <- split(colnames(srt), srt@meta.data[[split.by]])
  } else {
    split_cells <- list(colnames(srt))
  }
  srt_list <- vector("list", length(split_cells))
  names(srt_list) <- names(split_cells)

  for (i in seq_along(split_cells)) {
    if (!is.null(split.by)) {
      srt <- srt_raw[, split_cells[[i]]]
      log_message("Running QC for {.val {srt@meta.data[[split.by]][1]}}", verbose = verbose)
    } else {
      srt <- srt_raw
    }
    ntotal <- ncol(srt)

    db_qc <- c()
    if ("doublets" %in% qc_metrics) {
      if (!is.null(db_method)) {
        db_rate_use <- db_rate
        if (is.null(db_rate_use)) {
          db_rate_use <- ncol(srt) / 1000 * db_coefficient
        }
        if (db_rate_use >= 1) {
          log_message(
            "The db_rate is equal to or greater than 1",
            message_type = "error"
          )
        }
        for (dbm in db_method) {
          srt <- do.call(
            what = RunDoubletCalling,
            args = list(
              srt = srt,
              db_method = dbm,
              db_rate = db_rate_use,
              data_type = status
            )
          )
          db_qc <- unique(c(
            db_qc,
            colnames(srt)[
              srt[[paste0("db.", dbm, "_class"), drop = TRUE]] == "doublet"
            ]
          ))
        }
      }
    }

    atac_qc <- c()
    if ("atac" %in% qc_metrics) {
      if (!inherits(srt[[assay]], "ChromatinAssay")) {
        log_message(
          "Skip {.val atac} QC because {.arg assay = '{assay}'} is not a {.cls ChromatinAssay}",
          message_type = "warning",
          verbose = verbose
        )
      } else {
        atac_args_use <- atac_args
        min_pct_reads_in_peaks <- atac_args_use[["min_pct_reads_in_peaks"]] %||% NULL
        min_TSS_enrichment <- atac_args_use[["min_TSS_enrichment"]] %||% NULL
        max_nucleosome_signal <- atac_args_use[["max_nucleosome_signal"]] %||% NULL
        max_blacklist_ratio <- atac_args_use[["max_blacklist_ratio"]] %||% NULL
        atac_args_use[["min_pct_reads_in_peaks"]] <- NULL
        atac_args_use[["min_TSS_enrichment"]] <- NULL
        atac_args_use[["max_nucleosome_signal"]] <- NULL
        atac_args_use[["max_blacklist_ratio"]] <- NULL
        atac_args_use[["srt"]] <- srt
        atac_args_use[["assay"]] <- assay
        if (is.null(atac_args_use[["verbose"]])) {
          atac_args_use[["verbose"]] <- TRUE
        }
        srt <- do.call(
          what = RunATACQC,
          args = atac_args_use
        )
        atac_fail <- rep(FALSE, ncol(srt))
        names(atac_fail) <- colnames(srt)
        if (
          !is.null(min_pct_reads_in_peaks) &&
            "pct_reads_in_peaks" %in% colnames(srt@meta.data)
        ) {
          atac_fail <- atac_fail | srt$pct_reads_in_peaks < min_pct_reads_in_peaks
        }
        if (
          !is.null(min_TSS_enrichment) &&
            "TSS.enrichment" %in% colnames(srt@meta.data)
        ) {
          atac_fail <- atac_fail | srt$TSS.enrichment < min_TSS_enrichment
        }
        if (
          !is.null(max_nucleosome_signal) &&
            "nucleosome_signal" %in% colnames(srt@meta.data)
        ) {
          atac_fail <- atac_fail | srt$nucleosome_signal > max_nucleosome_signal
        }
        if (
          !is.null(max_blacklist_ratio) &&
            "blacklist_ratio" %in% colnames(srt@meta.data)
        ) {
          atac_fail <- atac_fail | (
            !is.na(srt$blacklist_ratio) &
              srt$blacklist_ratio > max_blacklist_ratio
          )
        }
        atac_qc <- names(atac_fail)[atac_fail]
      }
    }

    decontX_qc <- c()
    if ("decontX" %in% qc_metrics) {
      decontX_args_use <- decontX_args
      if (!group.by_missing) {
        decontX_args_use[["group.by"]] <- group.by
      }
      if (!decontX_batch_missing) {
        decontX_args_use[["batch"]] <- decontX_batch
      }
      if (!decontX_background_missing) {
        decontX_args_use[["background"]] <- decontX_background
      }
      if (!decontX_background_assay_missing) {
        decontX_args_use[["background_assay"]] <- decontX_background_assay
      }
      if (!decontX_bg_batch_missing) {
        decontX_args_use[["bg_batch"]] <- decontX_bg_batch
      }
      if (
        !decontX_assay_name_missing ||
          is.null(decontX_args_use[["assay_name"]])
      ) {
        decontX_args_use[["assay_name"]] <- decontX_assay_name_use
      }
      if (
        !decontX_store_assay_missing ||
          is.null(decontX_args_use[["store_assay"]])
      ) {
        decontX_args_use[["store_assay"]] <- decontX_store_assay_use
      }
      if (
        !decontX_round_counts_missing ||
          is.null(decontX_args_use[["round_counts"]])
      ) {
        decontX_args_use[["round_counts"]] <- decontX_round_counts_use
      }
      if (is.null(decontX_args_use[["seed"]])) {
        decontX_args_use[["seed"]] <- seed
      }
      decontX_args_use[["srt"]] <- srt
      decontX_args_use[["assay"]] <- assay
      decontX_args_use[["data_type"]] <- status
      srt <- do.call(
        what = RunDecontX,
        args = decontX_args_use
      )

      if (!is.null(decontX_threshold)) {
        decontX_qc <- colnames(srt)[which(
          srt[["decontX_contamination", drop = TRUE]] > decontX_threshold
        )]
      }
    }

    assay_genes <- rownames(Seurat::GetAssay(srt, assay = assay))
    counts_mat <- GetAssayData5(srt, assay = assay, layer = "counts")
    for (rule_name in names(qc_features)) {
      srt[[paste0("percent.", rule_name)]] <- cellqc_feature_percentage(
        counts = counts_mat,
        features = qc_features[[rule_name]][["features"]]
      )
    }
    outlier_qc <- c()
    for (n in 1:length(species)) {
      if (n == 0) {
        break
      }
      sp <- species[n]
      prefix <- species_gene_prefix[n]
      sp_genes <- assay_genes[grep(pattern = paste0("^", prefix), x = assay_genes)]
      nCount <- srt[[paste0(
        c(paste0("nCount_", assay), sp),
        collapse = "."
      )]] <- Matrix::colSums(
        counts_mat[sp_genes, , drop = FALSE]
      )
      nFeature <- srt[[paste0(
        c(paste0("nFeature_", assay), sp),
        collapse = "."
      )]] <- Matrix::colSums(
        counts_mat[sp_genes, , drop = FALSE] >
          0
      )
      percent.mito <- srt[[paste0(
        c("percent.mito", sp),
        collapse = "."
      )]] <- Seurat::PercentageFeatureSet(
        object = srt,
        assay = assay,
        pattern = paste0(
          "(",
          paste0("^", prefix, "-*", mito_pattern),
          ")",
          collapse = "|"
        ),
        features = mito_gene
      )
      percent.ribo <- srt[[paste0(
        c("percent.ribo", sp),
        collapse = "."
      )]] <- Seurat::PercentageFeatureSet(
        object = srt,
        assay = assay,
        pattern = paste0(
          "(",
          paste0("^", prefix, "-*", ribo_pattern),
          ")",
          collapse = "|"
        ),
        features = ribo_gene
      )
      if (!isTRUE(is_atac_assay)) {
        hb_features_use <- if (!is.null(hb_gene)) {
          intersect(hb_gene, sp_genes)
        } else {
          hb_pattern_use <- paste0(
            "(",
            paste0("^", prefix, "-*", hb_pattern),
            ")",
            collapse = "|"
          )
          assay_genes[grep(pattern = hb_pattern_use, x = assay_genes)]
        }
        if (n == 1 && "hb" %in% qc_metrics && length(hb_features_use) == 0) {
          log_message(
            "The hemoglobin feature definition did not match any assay features",
            message_type = "error"
          )
        }
        srt[[paste0(c("percent.hb", sp), collapse = ".")]] <-
          cellqc_feature_percentage(
            counts = counts_mat,
            features = hb_features_use
          )
      }
      percent.genome <- srt[[paste0(
        c("percent.genome", sp),
        collapse = "."
      )]] <- Seurat::PercentageFeatureSet(
        object = srt,
        assay = assay,
        pattern = paste0("^", prefix)
      )

      ribo.mito.ratio <- srt[[
        paste0(c("percent.ribo", sp), collapse = "."),
        drop = TRUE
      ]] /
        srt[[paste0(c("percent.mito", sp), collapse = "."), drop = TRUE]]
      ribo.mito.ratio[is.na(ribo.mito.ratio)] <- 1
      srt[[paste0(c("ribo.mito.ratio", sp), collapse = ".")]] <- ribo.mito.ratio

      if (n == 1) {
        if ("outlier" %in% qc_metrics) {
          log10_nFeature <- srt[[paste0(
            c(paste0("log10_nFeature_", assay), sp),
            collapse = "."
          )]] <- log10(nFeature)
          log10_nCount <- srt[[paste0(
            c(paste0("log10_nCount_", assay), sp),
            collapse = "."
          )]] <- log10(nCount)
          log10_nCount[is.infinite(log10_nCount)] <- NA
          log10_nFeature[is.infinite(log10_nFeature)] <- NA
          mod <- stats::loess(log10_nFeature ~ log10_nCount)
          pred <- stats::predict(
            mod,
            newdata = data.frame(log10_nCount = log10_nCount)
          )
          featurecount_dist <- srt[[paste0(
            c("featurecount_dist", sp),
            collapse = "."
          )]] <- log10_nFeature - pred

          var <- sapply(strsplit(outlier_threshold, ":"), function(x) x[[1]])
          var_valid <- var %in%
            colnames(srt@meta.data) |
            sapply(var, FUN = function(x) exists(x, where = environment()))
          if (any(!var_valid)) {
            log_message(
              "Variable {.val {names(var_valid)[!var_valid]}} is not found in the srt object",
              message_type = "error"
            )
          }
          outlier <- lapply(
            strsplit(outlier_threshold, ":"),
            function(m) {
              colnames(srt)[is_outlier(
                get(m[1]),
                nmads = as.numeric(m[3]),
                type = m[2]
              )]
            }
          )
          names(outlier) <- outlier_threshold
          outlier_tb <- table(unlist(outlier))
          outlier_qc <- c(
            outlier_qc,
            names(outlier_tb)[outlier_tb >= outlier_n]
          )
          for (nm in names(outlier)) {
            srt[[make.names(nm)]] <- colnames(srt) %in% outlier[[nm]]
          }
        }
      }
    }

    umi_qc <- gene_qc <- mito_qc <- ribo_qc <- hb_qc <-
      ribo_mito_ratio_qc <- species_qc <- c()
    if ("umi" %in% qc_metrics) {
      umi_qc <- colnames(srt)[which(
        srt[[
          paste0(c(paste0("nCount_", assay), species[1]), collapse = "."),
          drop = TRUE
        ]] <
          UMI_threshold
      )]
    }
    if ("gene" %in% qc_metrics) {
      gene_qc <- colnames(srt)[which(
        srt[[
          paste0(c(paste0("nFeature_", assay), species[1]), collapse = "."),
          drop = TRUE
        ]] <
          gene_threshold
      )]
    }
    if ("mito" %in% qc_metrics) {
      mito_qc <- colnames(srt)[which(
        srt[[
          paste0(c("percent.mito", species[1]), collapse = "."),
          drop = TRUE
        ]] >
          mito_threshold
      )]
    }
    if ("ribo" %in% qc_metrics) {
      ribo_qc <- colnames(srt)[which(
        srt[[
          paste0(c("percent.ribo", species[1]), collapse = "."),
          drop = TRUE
        ]] >
          ribo_threshold
      )]
    }
    if ("hb" %in% qc_metrics) {
      hb_percent <- srt[[
        paste0(c("percent.hb", species[1]), collapse = "."),
        drop = TRUE
      ]]
      hb_qc <- colnames(srt)[which(
        hb_percent < hb_range[1] |
          hb_percent > hb_range[2]
      )]
    }
    if ("ribo_mito_ratio" %in% qc_metrics) {
      ribo_mito_ratio_qc <- colnames(srt)[which(
        srt[[
          paste0(c("ribo.mito.ratio", species[1]), collapse = "."),
          drop = TRUE
        ]] <
          ribo_mito_ratio_range[1] |
          srt[[
            paste0(c("ribo.mito.ratio", species[1]), collapse = "."),
            drop = TRUE
          ]] >
            ribo_mito_ratio_range[2]
      )]
    }
    if ("species" %in% qc_metrics) {
      species_qc <- colnames(srt)[which(
        srt[[
          paste0(c("percent.genome", species[1]), collapse = "."),
          drop = TRUE
        ]] <
          species_percent
      )]
    }
    feature_qc_cells <- setNames(
      vector("list", length(active_qc_features)),
      active_qc_features
    )
    for (rule_name in active_qc_features) {
      feature_percent <- srt[[paste0("percent.", rule_name), drop = TRUE]]
      feature_range <- qc_features[[rule_name]][["range"]]
      feature_qc_cells[[rule_name]] <- colnames(srt)[which(
        feature_percent < feature_range[1] |
          feature_percent > feature_range[2]
      )]
    }

    CellQC <- unique(
      c(
        db_qc,
        atac_qc,
        decontX_qc,
        outlier_qc,
        umi_qc,
        gene_qc,
        mito_qc,
        ribo_qc,
        hb_qc,
        ribo_mito_ratio_qc,
        species_qc,
        unlist(feature_qc_cells, use.names = FALSE)
      )
    )
    qc_summary <- c(
      "{cli::symbol$record} Total cells: {.pkg {ntotal}}\n",
      "{cli::symbol$circle_filled} {.pkg {ntotal - length(CellQC)}} cells remained\n",
      "{cli::symbol$circle} {.pkg {length(CellQC)}} cells filtered out:\n",
      if ("doublets" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(db_qc)}} potential doublets\n"
      },
      if ("atac" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(atac_qc)}} ATAC QC failed cells\n"
      },
      if ("decontX" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(decontX_qc)}} high-contamination cells\n"
      },
      if ("outlier" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(outlier_qc)}} outlier cells\n"
      },
      if ("umi" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(umi_qc)}} low-UMI cells\n"
      },
      if ("gene" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(gene_qc)}} low-gene cells\n"
      },
      if ("mito" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(mito_qc)}} high-mito cells\n"
      },
      if ("ribo" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(ribo_qc)}} high-ribo cells\n"
      },
      if ("hb" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(hb_qc)}} HB feature-QC outlier cells\n"
      },
      if ("ribo_mito_ratio" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(ribo_mito_ratio_qc)}} ribo_mito_ratio outlier cells\n"
      },
      if ("species" %in% qc_metrics) {
        "{cli::symbol$circle}   {.pkg {length(species_qc)}} species-contaminated cells"
      },
      unname(vapply(active_qc_features, function(rule_name) {
        paste0(
          "{cli::symbol$circle}   ",
          length(feature_qc_cells[[rule_name]]),
          " ",
          rule_name,
          " feature-QC outlier cells\n"
        )
      }, character(1)))
    )
    qc_summary <- Filter(Negate(is.null), qc_summary)
    do.call(
      what = log_message,
      args = c(
        as.list(qc_summary),
        message_type = "success"
      )
    )

    qc_nm <- c(
      if ("doublets" %in% qc_metrics) {
        "db_qc"
      },
      if ("atac" %in% qc_metrics) {
        "atac_qc"
      },
      if ("decontX" %in% qc_metrics) {
        "decontX_qc"
      },
      if ("outlier" %in% qc_metrics) {
        "outlier_qc"
      },
      if ("umi" %in% qc_metrics) {
        "umi_qc"
      },
      if ("gene" %in% qc_metrics) {
        "gene_qc"
      },
      if ("mito" %in% qc_metrics) {
        "mito_qc"
      },
      if ("ribo" %in% qc_metrics) {
        "ribo_qc"
      },
      if ("hb" %in% qc_metrics) {
        "hb_qc"
      },
      if ("ribo_mito_ratio" %in% qc_metrics) {
        "ribo_mito_ratio_qc"
      },
      if ("species" %in% qc_metrics) {
        "species_qc"
      },
      "CellQC"
    )
    for (qc in qc_nm) {
      srt[[qc]] <- ifelse(
        colnames(srt) %in% get(qc),
        "Fail",
        "Pass"
      )
      srt[[qc]] <- factor(
        srt[[qc, drop = TRUE]],
        levels = c("Pass", "Fail")
      )
    }
    feature_qc_nm <- paste0(active_qc_features, "_qc")
    for (i_rule in seq_along(active_qc_features)) {
      rule_name <- active_qc_features[[i_rule]]
      qc_col <- feature_qc_nm[[i_rule]]
      srt[[qc_col]] <- factor(
        ifelse(
          colnames(srt) %in% feature_qc_cells[[rule_name]],
          "Fail",
          "Pass"
        ),
        levels = c("Pass", "Fail")
      )
    }
    qc_nm <- c(qc_nm, feature_qc_nm)

    if (return_filtered) {
      srt <- srt[, srt$CellQC == "Pass"]
      srt@meta.data[, intersect(qc_nm, colnames(srt@meta.data))] <- NULL
    }
    srt_list[[i]] <- srt
  }
  cells <- unlist(lapply(srt_list, colnames))
  srt_raw <- srt_raw[, cells]
  meta.data <- do.call(
    rbind.data.frame,
    unname(lapply(srt_list, function(x) x@meta.data))
  )
  srt_raw <- Seurat::AddMetaData(
    srt_raw,
    metadata = meta.data
  )

  if ("decontX" %in% qc_metrics && isTRUE(decontX_store_assay_use)) {
    decontX_assay_name <- decontX_assay_name_use
    decontX_mats <- unname(lapply(srt_list, function(x) {
      if (decontX_assay_name %in% SeuratObject::Assays(x)) {
        GetAssayData5(
          object = x,
          assay = decontX_assay_name,
          layer = "counts"
        )
      } else {
        NULL
      }
    }))
    decontX_mats <- Filter(Negate(is.null), decontX_mats)
    if (length(decontX_mats) > 0) {
      features_all <- unique(unlist(lapply(decontX_mats, rownames)))
      decontX_mats <- lapply(decontX_mats, function(mat) {
        if (identical(rownames(mat), features_all)) {
          return(mat)
        }
        mat_full <- Matrix::Matrix(
          0,
          nrow = length(features_all),
          ncol = ncol(mat),
          sparse = TRUE
        )
        rownames(mat_full) <- features_all
        colnames(mat_full) <- colnames(mat)
        mat_full[rownames(mat), colnames(mat)] <- mat
        mat_full
      })
      decontX_mat <- do.call(cbind, decontX_mats)
      decontX_mat <- decontX_mat[, colnames(srt_raw), drop = FALSE]
      srt_raw[[decontX_assay_name]] <- Seurat::CreateAssayObject(
        counts = decontX_mat
      )
    }
  }

  if ("decontX" %in% qc_metrics) {
    decontX_tools <- unname(lapply(srt_list, function(x) x@tools[["decontX"]]))
    if (length(decontX_tools) == 1) {
      srt_raw@tools[["decontX"]] <- decontX_tools[[1]]
    } else if (length(decontX_tools) > 1) {
      names(decontX_tools) <- names(srt_list)
      srt_raw@tools[["decontX"]] <- decontX_tools
    }
  }
  if ("doublets" %in% qc_metrics && "Scrublet" %in% db_method) {
    scrublet_tools <- unname(lapply(srt_list, function(x) {
      x@tools[["Scrublet"]]
    }))
    scrublet_tools <- Filter(Negate(is.null), scrublet_tools)
    if (length(scrublet_tools) == 1) {
      srt_raw@tools[["Scrublet"]] <- scrublet_tools[[1]]
    } else if (length(scrublet_tools) > 1) {
      names(scrublet_tools) <- names(srt_list)[seq_along(scrublet_tools)]
      srt_raw@tools[["Scrublet"]] <- scrublet_tools
    }
  }

  return(srt_raw)
}

cellqc_initialize_count_metrics <- function(srt, assay) {
  ncount_col <- paste0("nCount_", assay)
  nfeature_col <- paste0("nFeature_", assay)
  need_ncount <- !ncount_col %in% colnames(srt@meta.data)
  need_nfeature <- !nfeature_col %in% colnames(srt@meta.data)
  if (need_ncount || need_nfeature) {
    counts <- GetAssayData5(srt, assay = assay, layer = "counts")
    if (need_ncount) {
      srt@meta.data[[ncount_col]] <- Matrix::colSums(counts)
    }
    if (need_nfeature) {
      srt@meta.data[[nfeature_col]] <- Matrix::colSums(counts > 0)
    }
  }
  srt
}

cellqc_validate_range <- function(x, arg) {
  if (
    !is.numeric(x) ||
      length(x) != 2 ||
      anyNA(x) ||
      x[[1]] > x[[2]]
  ) {
    stop(
      sprintf("%s must be a numeric vector of length two with lower <= upper", arg),
      call. = FALSE
    )
  }
  as.numeric(x)
}

cellqc_feature_percentage <- function(counts, features) {
  totals <- Matrix::colSums(counts)
  feature_totals <- Matrix::colSums(counts[features, , drop = FALSE])
  percent <- numeric(length(totals))
  nonzero <- totals > 0
  percent[nonzero] <- feature_totals[nonzero] / totals[nonzero] * 100
  names(percent) <- colnames(counts)
  percent
}
