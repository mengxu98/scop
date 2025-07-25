#' @title GetFeaturesData
#' @description Get the data from the \pkg{Seurat} object
#'
#' @param object A object
#' @param ... Additional arguments passed to the method
#' @return data
#' @export
GetFeaturesData <- function(object, ...) {
  UseMethod(generic = "GetFeaturesData", object = object)
}

#' @param assay Assay to get data from
#' @export
#' @rdname GetFeaturesData
#' @method GetFeaturesData Seurat
#' @examples
#' data(pancreas_sub)
#' features <- GetFeaturesData(pancreas_sub)
#' head(features)
GetFeaturesData.Seurat <- function(
    object,
    assay = NULL,
    ...) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  assay_obj <- Seurat::GetAssay(
    object,
    assay = assay
  )
  return(GetFeaturesData(assay_obj))
}

#' @rdname GetFeaturesData
#' @method GetFeaturesData Assay
#' @export
GetFeaturesData.Assay <- function(
    object,
    ...) {
  return(object@meta.features)
}

#' @rdname GetFeaturesData
#' @method GetFeaturesData Assay5
#' @export
GetFeaturesData.Assay5 <- function(
    object,
    ...) {
  misc <- object@misc
  meta.features <- data.frame(
    row.names = rownames(object)
  ) |> as.data.frame()

  if (is.null(misc)) {
    return(meta.features)
  } else {
    if ("meta.features" %in% names(misc)) {
      return(as.data.frame(misc[["meta.features"]]))
    } else {
      return(meta.features)
    }
  }
}

#' @title AddFeaturesData
#' @description Add the data to the Seurat object
#' @param object A object
#' @param assay Assay to add data to
#' @param features Features to add data to
#' @param ... Additional arguments passed to the method
#' @return data
#' @export
AddFeaturesData <- function(object, ...) {
  UseMethod(
    generic = "AddFeaturesData",
    object = object
  )
}

#' @rdname AddFeaturesData
#' @method AddFeaturesData Seurat
#' @export
#' @examples
#' data(pancreas_sub)
#' features <- GetFeaturesData(pancreas_sub)
#' pancreas_sub <- AddFeaturesData(pancreas_sub, features)
AddFeaturesData.Seurat <- function(
    object,
    features,
    assay = NULL,
    ...) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  assay_obj <- Seurat::GetAssay(
    object,
    assay = assay
  )
  assay_obj <- AddFeaturesData(assay_obj, features, ...)
  object[[assay]] <- assay_obj

  return(object)
}

#' @rdname AddFeaturesData
#' @method AddFeaturesData Assay
#' @export
AddFeaturesData.Assay <- function(
    object,
    features,
    ...) {
  object@meta.features <- features
  return(object)
}

#' @rdname AddFeaturesData
#' @method AddFeaturesData Assay5
#' @export
AddFeaturesData.Assay5 <- function(
    object,
    features,
    ...) {
  object@misc$meta.features <- features
  return(object)
}
