#' @title Get expression data from Assay5 or Seurat object
#'
#' @description
#' A re-implementation of the [SeuratObject::GetAssayData] function to compatible with Assay5 objects.
#'
#' @md
#' @inheritParams SeuratObject::GetAssayData
#' @param verbose Whether to print messages. Default is `TRUE`.
#' @param ... Additional arguments passed to [SeuratObject::GetAssayData].
#'
#' @return A matrix or data frame containing the assay data.
#' @export
#' @seealso [SeuratObject::GetAssayData]
#'
#' @examples
#' data(pancreas_sub)
#' GetAssayData5(
#'   pancreas_sub,
#'   layer = "counts",
#'   assay = "RNA"
#' )[1:5, 1:5]
#'
#' data(panc8_sub)
#' GetAssayData5(
#'   panc8_sub,
#'   layer = "counts",
#'   assay = "RNA"
#' )[1:5, 1:5]
GetAssayData5 <- function(object, ...) {
  UseMethod(generic = "GetAssayData5", object = object)
}

#' @rdname GetAssayData5
#' @method GetAssayData5 Seurat
#' @export
GetAssayData5.Seurat <- function(
    object,
    layer = "counts",
    assay = NULL,
    verbose = TRUE,
    ...) {
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  assay_obj <- Seurat::GetAssay(
    object = object,
    assay = assay,
    ...
  )
  data <- GetAssayData5(
    object = assay_obj,
    layer = layer,
    verbose = verbose,
    ...
  )
  return(data)
}

#' @param join_layers Whether to join layers if the object is an Assay5 object. Default is `TRUE`.
#' @rdname GetAssayData5
#' @method GetAssayData5 Assay5
#' @export
GetAssayData5.Assay5 <- function(
    object,
    layer = "counts",
    join_layers = TRUE,
    verbose = TRUE,
    ...) {
  if (join_layers) {
    if (verbose) {
      warning_key <- "assay5_join_layers_warning"
      last_warning_time <- if (exists(warning_key, envir = .scop_env)) {
        get(warning_key, envir = .scop_env)
      } else {
        as.numeric(Sys.time()) - 30000
      }

      current_time <- as.numeric(Sys.time())
      warning_interval <- getOption("scop.warning_interval", 28800)
      if (current_time - last_warning_time > warning_interval) {
        log_message(
          "The input data is a {.cls Assay5} object\n",
          "{.fn SeuratObject::JoinLayers} will be used to combine the layers\n",
          cli::col_grey(
            "This warning will be shown only once every {warning_interval / 3600} hours"
          ), "\n",
          cli::col_grey(
            "To change this interval, set the scop.warning_interval option"
          ),
          message_type = "warning"
        )
        assign(warning_key, current_time, envir = .scop_env)
      }
    }

    object <- SeuratObject::JoinLayers(object)
    data <- SeuratObject::GetAssayData(
      object,
      layer = layer,
      ...
    )
  } else {
    data <- SeuratObject::LayerData(
      object,
      layer = layer,
      ...
    )
  }

  return(data)
}

#' @rdname GetAssayData5
#' @method GetAssayData5 Assay
#' @export
GetAssayData5.Assay <- function(
    object,
    layer = "counts",
    verbose = TRUE,
    ...) {
  log_message(
    "Get expression data from {.cls Assay} object",
    verbose = verbose
  )
  SeuratObject::GetAssayData(
    object,
    layer = layer,
    ...
  )
}
