#' @title GetAssayData5
#'
#' @description A re-implementation of the GetAssayData function to compatible with Assay5 objects.
#'
#' @md
#' @inheritParams SeuratObject::GetAssayData
#' @param verbose Logical, whether to print messages.
#' @param ... Additional arguments passed to [SeuratObject::GetAssayData].
#'
#' @return A matrix or data frame containing the assay data.
#' @export
#' @seealso [SeuratObject::GetAssayData]
#'
#' @examples
#' data("pancreas_sub")
#' GetAssayData5(
#'   pancreas_sub,
#'   layer = "counts",
#'   assay = "RNA"
#' )[1:5, 1:5]
#'
#' data("panc8_sub")
#' GetAssayData5(
#'   panc8_sub,
#'   layer = "counts",
#'   assay = "RNA"
#' )[1:5, 1:5]
GetAssayData5 <- function(
    object,
    layer = "counts",
    assay = NULL,
    verbose = TRUE,
    ...) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  assay_obj <- Seurat::GetAssay(object, assay = assay)

  if (inherits(assay_obj, "Assay5")) {
    if (verbose) {
      warning(
        "The input data is a 'Assay5' object. The 'SeuratObject::JoinLayers' function will be used to combine the layers.",
        immediate. = TRUE
      )
      object <- SeuratObject::JoinLayers(object)
    }
  }

  data <- SeuratObject::GetAssayData(
    object,
    layer = layer,
    assay = assay,
    ...
  )

  return(data)
}
