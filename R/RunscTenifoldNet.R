#' @title Run scTenifoldNet network comparison
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param object A `Seurat` object or a raw count matrix with genes in rows and
#' cells in columns.
#' @param y A second raw count matrix. Required when `object` is a matrix and
#' ignored when `object` is a `Seurat` object.
#' @param group.by Metadata column used to split a `Seurat` object into the two
#' conditions being compared.
#' @param condition1,condition2 Condition labels from `group.by`. If omitted,
#' the first two group levels are used.
#' @param assay,layer Assay and layer used as the count matrix when `object` is
#' a `Seurat` object.
#' @param features Optional genes to retain before running the comparison.
#' @param qc Whether to apply scTenifoldNet-style quality control.
#' @param qc_min_library_size,qc_remove_outlier_cells,qc_min_pct,qc_max_mt_ratio
#' Quality-control parameters forwarded to `scTenifoldNet::scQC()`.
#' @param nc_nNet,nc_nCells,nc_nComp,nc_symmetric,nc_scaleScores,nc_q Network
#' construction parameters forwarded to `scTenifoldNet::makeNetworks()`.
#' @param td_K,td_nDecimal,td_maxIter,td_maxError Tensor decomposition
#' parameters forwarded to `scTenifoldNet::tensorDecomposition()`.
#' @param ma_nDim Manifold-alignment dimension forwarded to
#' `scTenifoldNet::manifoldAlignment()`.
#' @param cores Number of cores forwarded to `scTenifoldNet::scTenifoldNet()`.
#' @param store_networks Whether to keep tensor networks in the stored result
#' when `object` is a `Seurat` object.
#' @param store_manifold Whether to keep manifold-alignment coordinates in the
#' stored result when `object` is a `Seurat` object.
#' @param tool_name Name of the `object@tools` entry when `object` is a `Seurat`
#' object.
#'
#' @return A scTenifoldNet result list for matrix input, or a `Seurat` object
#' with results stored in `object@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' counts <- GetAssayData5(pancreas_sub, assay = "RNA", layer = "counts")
#' detected <- names(sort(Matrix::rowSums(counts > 0), decreasing = TRUE))
#' features_use <- head(detected, 300)
#'
#' pancreas_sub <- RunscTenifoldNet(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   condition1 = "Ductal",
#'   condition2 = "Endocrine",
#'   features = features_use,
#'   qc = FALSE,
#'   nc_nNet = 3,
#'   nc_nCells = 200,
#'   td_maxIter = 200,
#'   ma_nDim = 2,
#'   store_networks = FALSE,
#'   store_manifold = TRUE
#' )
#'
#' dr <- pancreas_sub@tools$scTenifoldNet$diffRegulation
#' head(dr)
#'
#' scTenifoldNetPlot(pancreas_sub, plot_type = "effect")
RunscTenifoldNet <- function(
  object,
  y = NULL,
  group.by = NULL,
  condition1 = NULL,
  condition2 = NULL,
  assay = NULL,
  layer = "counts",
  features = NULL,
  qc = TRUE,
  qc_min_library_size = 1000,
  qc_remove_outlier_cells = TRUE,
  qc_min_pct = 0.05,
  qc_max_mt_ratio = 0.1,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_symmetric = FALSE,
  nc_scaleScores = TRUE,
  nc_q = 0.05,
  td_K = 3,
  td_nDecimal = 1,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  ma_nDim = 30,
  cores = 1,
  store_networks = TRUE,
  store_manifold = TRUE,
  tool_name = "scTenifoldNet",
  verbose = TRUE
) {
  numeric_params <- list(
    nc_q = c(value = nc_q, lower = 0, upper = 1),
    td_maxError = c(value = td_maxError, lower = 0, upper = 1),
    qc_min_pct = c(value = qc_min_pct, lower = 0, upper = 1),
    qc_max_mt_ratio = c(value = qc_max_mt_ratio, lower = 0, upper = 1)
  )
  for (param_name in names(numeric_params)) {
    param <- numeric_params[[param_name]]
    if (
      length(param[["value"]]) != 1L ||
        !is.numeric(param[["value"]]) ||
        is.na(param[["value"]]) ||
        param[["value"]] < param[["lower"]] ||
        param[["value"]] > param[["upper"]]
    ) {
      log_message(
        "{.arg {param_name}} must be a single number between {.val {param[['lower']]}} and {.val {param[['upper']]}}",
        message_type = "error"
      )
    }
  }
  integer_params <- list(
    nc_nNet = c(value = nc_nNet, lower = 1),
    nc_nCells = c(value = nc_nCells, lower = 1),
    nc_nComp = c(value = nc_nComp, lower = 2),
    td_K = c(value = td_K, lower = 1),
    td_nDecimal = c(value = td_nDecimal, lower = 0),
    td_maxIter = c(value = td_maxIter, lower = 1),
    ma_nDim = c(value = ma_nDim, lower = 1),
    cores = c(value = cores, lower = 1)
  )
  for (param_name in names(integer_params)) {
    param <- integer_params[[param_name]]
    if (
      length(param[["value"]]) != 1L ||
        !is.numeric(param[["value"]]) ||
        is.na(param[["value"]]) ||
        param[["value"]] < param[["lower"]] ||
        param[["value"]] != as.integer(param[["value"]])
    ) {
      log_message(
        "{.arg {param_name}} must be a single integer greater than or equal to {.val {param[['lower']]}}",
        message_type = "error"
      )
    }
  }

  is_seurat <- inherits(object, "Seurat")
  assay <- assay %||% if (is_seurat) SeuratObject::DefaultAssay(object) else NULL
  input_summary <- list(type = if (is_seurat) "Seurat" else class(object)[[1]])

  if (is_seurat) {
    if (is.null(group.by) || !group.by %in% colnames(object[[]])) {
      log_message(
        "{.arg group.by} must identify a metadata column when {.arg object} is a {.cls Seurat} object",
        message_type = "error"
      )
    }
    groups <- object[[group.by]][, 1]
    condition_levels <- levels(groups) %||% unique(as.character(groups))
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
    if (
      is.null(condition1) ||
        is.null(condition2) ||
        identical(condition1, condition2)
    ) {
      log_message(
        "{.arg condition1} and {.arg condition2} must identify two different groups",
        message_type = "error"
      )
    }
    cells_x <- rownames(object[[]])[as.character(groups) == condition1]
    cells_y <- rownames(object[[]])[as.character(groups) == condition2]
    if (length(cells_x) == 0L || length(cells_y) == 0L) {
      log_message(
        "Both selected conditions must contain cells",
        message_type = "error"
      )
    }
    count_matrix <- GetAssayData5(object, assay = assay, layer = layer)
    x <- count_matrix[, cells_x, drop = FALSE]
    y <- count_matrix[, cells_y, drop = FALSE]
    input_summary$group.by <- group.by
    input_summary$condition1 <- condition1
    input_summary$condition2 <- condition2
    input_summary$cells <- stats::setNames(
      c(length(cells_x), length(cells_y)),
      c(condition1, condition2)
    )
  } else {
    x <- object
    if (is.null(y)) {
      log_message(
        "{.arg y} is required when {.arg object} is a matrix",
        message_type = "error"
      )
    }
  }

  if (is.null(rownames(x)) || is.null(rownames(y))) {
    log_message(
      "Both input matrices must contain gene names as row names",
      message_type = "error"
    )
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X_cell_", seq_len(ncol(x)))
  }
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("Y_cell_", seq_len(ncol(y)))
  }

  if (!is.null(features)) {
    features <- unique(as.character(features))
    shared_features <- Reduce(intersect, list(features, rownames(x), rownames(y)))
    if (length(shared_features) == 0L) {
      log_message(
        "No requested {.arg features} are present in both inputs",
        message_type = "error"
      )
    }
    missing_features <- setdiff(features, shared_features)
    if (length(missing_features) > 0L) {
      log_message(
        "Ignoring {.val {length(missing_features)}} requested features absent from at least one input",
        message_type = "warning",
        verbose = verbose
      )
    }
    x <- x[shared_features, , drop = FALSE]
    y <- y[shared_features, , drop = FALSE]
  }

  if (!inherits(x, "Matrix")) {
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  }
  if (!inherits(y, "Matrix")) {
    y <- Matrix::Matrix(as.matrix(y), sparse = TRUE)
  }
  x <- methods::as(x, "dgCMatrix")
  y <- methods::as(y, "dgCMatrix")

  if (nrow(x) <= nc_nComp || nrow(y) <= nc_nComp) {
    log_message(
      "{.arg nc_nComp} must be lower than the number of genes in both inputs",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg scTenifoldNet} comparison using upstream implementation",
    verbose = verbose
  )

  check_r("cailab-tamu/scTenifoldNet", verbose = FALSE)
  result <- get_namespace_fun("scTenifoldNet", "scTenifoldNet")(
    X = x,
    Y = y,
    qc = qc,
    qc_minLibSize = as.integer(qc_min_library_size),
    qc_removeOutlierCells = isTRUE(qc_remove_outlier_cells),
    qc_minPCT = qc_min_pct,
    qc_maxMTratio = qc_max_mt_ratio,
    nc_nNet = as.integer(nc_nNet),
    nc_nCells = as.integer(nc_nCells),
    nc_nComp = as.integer(nc_nComp),
    nc_symmetric = isTRUE(nc_symmetric),
    nc_scaleScores = isTRUE(nc_scaleScores),
    nc_q = nc_q,
    td_K = as.integer(td_K),
    td_nDecimal = as.integer(td_nDecimal),
    td_maxIter = as.integer(td_maxIter),
    td_maxError = td_maxError,
    ma_nDim = as.integer(ma_nDim),
    nCores = as.integer(cores)
  )

  if (!is.list(result) || is.null(result$diffRegulation)) {
    log_message(
      "{.pkg scTenifoldNet} did not return a valid {.field diffRegulation} table",
      message_type = "error"
    )
  }
  result$qc_summary <- result$qc_summary %||% list(applied = qc)

  parameters <- list(
    assay = assay,
    layer = layer,
    features = features,
    qc = qc,
    qc_min_library_size = qc_min_library_size,
    qc_remove_outlier_cells = qc_remove_outlier_cells,
    qc_min_pct = qc_min_pct,
    qc_max_mt_ratio = qc_max_mt_ratio,
    nc_nNet = nc_nNet,
    nc_nCells = nc_nCells,
    nc_nComp = nc_nComp,
    nc_symmetric = nc_symmetric,
    nc_scaleScores = nc_scaleScores,
    nc_q = nc_q,
    td_K = td_K,
    td_nDecimal = td_nDecimal,
    td_maxIter = td_maxIter,
    td_maxError = td_maxError,
    ma_nDim = ma_nDim,
    cores = cores
  )

  if (!is_seurat) {
    result$parameters <- parameters
    result$input_summary <- input_summary
    return(result)
  }

  stored_result <- result
  if (!isTRUE(store_networks)) {
    stored_result$tensorNetworks <- NULL
  }
  if (!isTRUE(store_manifold)) {
    stored_result$manifoldAlignment <- NULL
  }
  object@tools[[tool_name]] <- list(
    diffRegulation = result$diffRegulation,
    result = stored_result,
    qc_summary = result$qc_summary %||% list(applied = qc),
    input_summary = input_summary,
    parameters = c(
      parameters,
      list(
        store_networks = store_networks,
        store_manifold = store_manifold
      )
    )
  )

  log_message(
    "{.pkg scTenifoldNet} results stored in {.code object@tools[[{tool_name}]]}",
    message_type = "success",
    verbose = verbose
  )
  object
}
