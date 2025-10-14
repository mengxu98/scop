#' @title Get expression data from `Assay5` or Seurat object
#'
#' @description
#' A re-implementation of the [SeuratObject::GetAssayData] function to compatible with Assay5 objects.
#'
#' @md
#' @inheritParams SeuratObject::GetAssayData
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
    ...
  )
  return(data)
}

#' @rdname GetAssayData5
#' @method GetAssayData5 Assay5
#' @export
GetAssayData5.Assay5 <- function(
    object,
    layer = "counts",
    ...) {
  object <- SeuratObject::JoinLayers(object)
  data <- SeuratObject::GetAssayData(
    object,
    layer = layer,
    ...
  )

  return(data)
}

#' @rdname GetAssayData5
#' @method GetAssayData5 Assay
#' @export
GetAssayData5.Assay <- function(
    object,
    layer = "counts",
    ...) {
  SeuratObject::GetAssayData(
    object,
    layer = layer,
    ...
  )
}
