#' @title Run SpaNorm spatial normalization
#'
#' @description
#' Normalize spatial transcriptomics counts with the optional Bioconductor
#' `SpaNorm` backend and store the normalized expression in a new Seurat assay.
#'
#' @md
#' @inheritParams RunSpatialVariableFeatures
#' @inheritParams thisutils::log_message
#' @param new_assay Name of the assay used to store SpaNorm-normalized data.
#' @param tool_name Name used to store detailed SpaNorm results in `srt@tools`.
#' @param store_spe Whether to store the backend `SpatialExperiment` returned
#' by `SpaNorm`.
#' @param ... Additional arguments passed to `SpaNorm::SpaNorm()`, such as
#' `sample.p`.
#'
#' @return A `Seurat` object with SpaNorm-normalized expression stored in
#' `new_assay`. When `store_results = TRUE`, parameters, coordinates, features,
#' cells, and optional backend output are stored in `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- subset(
#'   visium_human_pancreas_sub,
#'   cells = colnames(visium_human_pancreas_sub)[1:80],
#'   features = rownames(visium_human_pancreas_sub)[1:300]
#' )
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#'
#' SpatialSpotPlot(
#'   spatial,
#'   features = rownames(spatial)[1:2],
#'   assay = "Spatial",
#'   layer = "data",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' if (
#'   requireNamespace("SpaNorm", quietly = TRUE) &&
#'     requireNamespace("SpatialExperiment", quietly = TRUE) &&
#'     identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
#' ) {
#'   spatial <- RunSpaNorm(
#'     spatial,
#'     assay = "Spatial",
#'     layer = "counts",
#'     coord.cols = c("x", "y"),
#'     new_assay = "SpaNorm",
#'     store_spe = FALSE,
#'     sample.p = 0.25,
#'     verbose = FALSE
#'   )
#'
#'   SpatialSpotPlot(
#'     spatial,
#'     features = rownames(spatial[["SpaNorm"]])[1:2],
#'     assay = "SpaNorm",
#'     layer = "data",
#'     overlay_image = FALSE,
#'     coord.cols = c("x", "y")
#'   )
#' }
RunSpaNorm <- function(
  srt,
  assay = NULL,
  layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  new_assay = "SpaNorm",
  tool_name = "SpaNorm",
  store_results = TRUE,
  store_spe = FALSE,
  verbose = TRUE,
  ...
) {
  log_message(
    "Running SpaNorm spatial normalization",
    message_type = "running",
    verbose = verbose
  )
  spanorm_validate_srt(srt)
  spanorm_assert_string(new_assay, "new_assay")
  spanorm_assert_string(tool_name, "tool_name")
  spanorm_assert_flag(store_results, "store_results")
  spanorm_assert_flag(store_spe, "store_spe")

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }

  input <- spanorm_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    image = image,
    coord.cols = coord.cols
  )
  backend_args <- list(...)
  spanorm_validate_named_args(backend_args)

  spe <- spanorm_make_spe(
    counts = input$counts,
    coords = input$coords
  )
  log_message(
    "Run {.pkg SpaNorm} with {.val {nrow(input$counts)}} features and {.val {ncol(input$counts)}} spots",
    verbose = verbose
  )
  result <- spanorm_run_backend(spe, backend_args)
  normalized <- spanorm_extract_logcounts(
    result = result,
    features = rownames(input$counts),
    cells = colnames(input$counts)
  )
  suppressWarnings(
    srt[[new_assay]] <- SeuratObject::CreateAssayObject(data = normalized)
  )

  if (isTRUE(store_results)) {
    tool <- list(
      parameters = list(
        assay = assay,
        layer = layer,
        image = image,
        coord.cols = coord.cols,
        new_assay = new_assay,
        tool_name = tool_name,
        store_results = store_results,
        store_spe = store_spe,
        backend_args = backend_args
      ),
      coords = input$coords,
      features = rownames(input$counts),
      cells = colnames(input$counts)
    )
    if (isTRUE(store_spe)) {
      tool$spe <- result
    }
    srt@tools[[tool_name]] <- tool
  }

  log_message(
    "{.pkg SpaNorm} normalized expression stored in assay {.val {new_assay}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

spanorm_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spanorm_prepare_input <- function(
  srt,
  assay,
  layer,
  image = NULL,
  coord.cols = c("col", "row")
) {
  coords <- spatial_dim_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    overlay_image = FALSE
  )$data
  cells <- colnames(srt)[colnames(srt) %in% rownames(coords)]
  if (length(cells) == 0L) {
    log_message(
      "No spatial coordinates match spots in {.arg srt}",
      message_type = "error"
    )
  }
  coords <- coords[cells, , drop = FALSE]
  keep <- is.finite(coords$x) & is.finite(coords$y)
  if (!all(keep)) {
    cells <- cells[keep]
    coords <- coords[cells, , drop = FALSE]
  }
  if (length(cells) < 3L) {
    log_message(
      "At least three spots with finite coordinates are required for {.fn RunSpaNorm}",
      message_type = "error"
    )
  }

  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  counts <- counts[, cells, drop = FALSE]
  if (!inherits(counts, "Matrix")) {
    counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "CsparseMatrix")
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }
  counts@x[!is.finite(counts@x) | counts@x < 0] <- 0
  counts <- Matrix::drop0(counts)
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    log_message(
      "Input expression matrix is empty after coordinate matching",
      message_type = "error"
    )
  }
  list(counts = counts, coords = coords)
}

spanorm_make_spe <- function(counts, coords) {
  check_r(c("SpatialExperiment", "SummarizedExperiment", "S4Vectors"), verbose = FALSE)
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = as.matrix(coords[, c("x", "y"), drop = FALSE])
  )
}

spanorm_run_backend <- function(spe, backend_args = list()) {
  check_r("SpaNorm", verbose = FALSE)
  spanorm_fun <- get_namespace_fun("SpaNorm", "SpaNorm")
  do.call(spanorm_fun, c(list(spe), backend_args))
}

spanorm_extract_logcounts <- function(result, features, cells) {
  if (!methods::is(result, "SummarizedExperiment")) {
    log_message(
      "{.pkg SpaNorm} must return a {.cls SpatialExperiment} or {.cls SummarizedExperiment} object",
      message_type = "error"
    )
  }
  assay_names <- SummarizedExperiment::assayNames(result)
  if (!"logcounts" %in% assay_names) {
    log_message(
      "{.pkg SpaNorm} result does not contain a {.val logcounts} assay",
      message_type = "error"
    )
  }
  mat <- SummarizedExperiment::assay(result, "logcounts")
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "CsparseMatrix")
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- features
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- cells
  }
  if (!all(features %in% rownames(mat)) || !all(cells %in% colnames(mat))) {
    log_message(
      "{.pkg SpaNorm} logcounts could not be matched to input features and spots",
      message_type = "error"
    )
  }
  mat[features, cells, drop = FALSE]
}

spanorm_assert_string <- function(x, arg) {
  if (
    is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
  ) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spanorm_assert_flag <- function(x, arg) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    log_message(
      "{.arg {arg}} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spanorm_validate_named_args <- function(x) {
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message(
      "{.arg ...} must contain named arguments only",
      message_type = "error"
    )
  }
  invisible(TRUE)
}
