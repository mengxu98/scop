#' @title Run CytoTRACE 2
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams thisutils::parallelize_fun
#' @inheritParams standard_scop
#' @inheritParams GroupHeatmap
#' @param object An object.
#' This can be a Seurat object or a matrix-like object (genes as rows, cells as columns).
#' @param species The species of the input data.
#' Currently supported values are `"human"` and `"mouse"`.
#' Default is `"human"`.
#' @param batch_size The number of cells to process at once,
#' including subsampling for KNN smoothing.
#' No subsampling if `NULL`.
#' Default is `10000` (recommended for input data size > 10K cells).
#' @param smooth_batch_size The number of cells to subsample further within the batch_size for the smoothing by diffusion step of the pipeline.
#' No subsampling if `NULL`.
#' Default is `1000` (recommended for input data size > 1K cells).
#' @param parallelize_models Whether to run the prediction function on models in parallel on multiple threads.
#' Default is `TRUE`.
#' @param parallelize_smoothing Whether to run the smoothing function on subsamples in parallel on multiple threads.
#' Default is `TRUE`.
#' @param ... Additional arguments to be passed to `CytoTRACE2::cytotrace2`.
#'
#' @rdname RunCytoTRACE
#' @export
#'
#' @return
#' When the input is a Seurat object,
#' the function returns a Seurat object with the following metadata columns added:
#' \itemize{
#'   \item \code{CytoTRACE2_Score}: The final predicted cellular potency score (0-1)
#'   \item \code{CytoTRACE2_Potency}: The final predicted cellular potency category
#'     (Differentiated, Unipotent, Oligopotent, Multipotent, Pluripotent, Totipotent)
#'   \item \code{CytoTRACE2_Relative}: The predicted relative order (normalized to 0-1)
#'   \item \code{preKNN_CytoTRACE2_Score}: The potency score before KNN smoothing
#'   \item \code{preKNN_CytoTRACE2_Potency}: The potency category before KNN smoothing
#' }
#'
#' When the input is a matrix or data.frame,
#' the function returns a data.frame with the same columns as above,
#' with cell IDs as row names.
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCytoTRACE(pancreas_sub, species = "mouse")
#' CytoTRACEPlot(pancreas_sub, group.by = "CellType")
RunCytoTRACE <- function(object, ...) {
  UseMethod(generic = "RunCytoTRACE", object = object)
}

#' @rdname RunCytoTRACE
#' @method RunCytoTRACE Seurat
#' @export
RunCytoTRACE.Seurat <- function(
    object,
    assay = NULL,
    layer = c("counts", "data"),
    species = c("human", "mouse"),
    batch_size = 10000,
    smooth_batch_size = 1000,
    parallelize_models = TRUE,
    parallelize_smoothing = TRUE,
    cores = 1,
    seed = 11,
    verbose = TRUE,
    ...) {
  log_message(
    "Running {.pkg CytoTRACE2}",
    message_type = "running",
    verbose = verbose
  )

  if (!check_pkg_status("CytoTRACE2")) {
    log_message(
      "Package {.pkg CytoTRACE2} is not installed. Installing from GitHub...",
      message_type = "info",
      verbose = verbose
    )
    pak::pak("digitalcytometry/cytotrace2/cytotrace2_r")
  }

  layer <- match.arg(layer)
  species <- match.arg(species)

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)

  cytotrace2_args <- list(
    input = object,
    species = species,
    is_seurat = TRUE,
    slot_type = layer,
    batch_size = batch_size,
    smooth_batch_size = smooth_batch_size,
    parallelize_models = parallelize_models,
    parallelize_smoothing = parallelize_smoothing,
    ncores = cores,
    seed = seed
  )

  cytotrace2_args <- c(cytotrace2_args, list(...))

  result <- do.call(
    get_namespace_fun("CytoTRACE2", "cytotrace2"),
    cytotrace2_args
  )

  result <- Seurat::LogSeuratCommand(object = result)

  log_message(
    "{.pkg CytoTRACE2} computed successfully",
    message_type = "success",
    verbose = verbose
  )

  return(result)
}

#' @rdname RunCytoTRACE
#' @method RunCytoTRACE default
#' @export
RunCytoTRACE.default <- function(
    object,
    species = c("human", "mouse"),
    batch_size = 10000,
    smooth_batch_size = 1000,
    parallelize_models = TRUE,
    parallelize_smoothing = TRUE,
    cores = 1,
    seed = 11,
    verbose = TRUE,
    ...) {
  log_message(
    "Running {.pkg CytoTRACE2}",
    message_type = "running",
    verbose = verbose
  )

  if (!check_pkg_status("CytoTRACE2")) {
    log_message(
      "Package {.pkg CytoTRACE2} is not installed. Installing from GitHub...",
      message_type = "info",
      verbose = verbose
    )
    pak::pak("digitalcytometry/cytotrace2/cytotrace2_r")
  }

  species <- match.arg(species)

  if (inherits(object, c("matrix", "Matrix"))) {
    object <- as.data.frame(object)
  }

  if (!is.data.frame(object)) {
    log_message(
      "{.arg object} must be a {.cls Seurat}, matrix, or data.frame",
      message_type = "error"
    )
  }

  if (is.null(rownames(object))) {
    log_message(
      "{.arg object} must have row names (gene names)",
      message_type = "error"
    )
  }
  if (is.null(colnames(object))) {
    log_message(
      "{.arg object} must have column names (cell IDs)",
      message_type = "error"
    )
  }

  cytotrace2_args <- list(
    input = object,
    species = species,
    is_seurat = FALSE,
    batch_size = batch_size,
    smooth_batch_size = smooth_batch_size,
    parallelize_models = parallelize_models,
    parallelize_smoothing = parallelize_smoothing,
    ncores = cores,
    seed = seed
  )

  cytotrace2_args <- c(cytotrace2_args, list(...))

  result <- do.call(
    get_namespace_fun("CytoTRACE2", "cytotrace2"),
    cytotrace2_args
  )

  log_message(
    "{.pkg CytoTRACE2} computed successfully",
    message_type = "success",
    verbose = verbose
  )

  return(result)
}
