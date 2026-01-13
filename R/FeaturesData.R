#' @title Get features data
#'
#' @description Get the data from the `Assay`, `Assay5` or `Seurat` object.
#'
#' @md
#' @inheritParams standard_scop
#' @param object A `Assay`, `Assay5` or `Seurat` object.
#' @param ... Additional arguments passed to the method.
#'
#' @return A data frame containing the features data.
#' @export
GetFeaturesData <- function(object, ...) {
  UseMethod(generic = "GetFeaturesData", object = object)
}

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
  return(object[[]])
}

#' @title Add features data
#'
#' @description
#' Add features data to the `Assay`, `Assay5` or `Seurat` object.
#'
#' @md
#' @inheritParams GetFeaturesData
#' @param features Features data to add.
#' @param ... Additional arguments passed to the method.
#'
#' @return A `Assay`, `Assay5` or `Seurat` object.
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
  object@meta.features <- check_features_data(object, features)
  return(object)
}

#' @rdname AddFeaturesData
#' @method AddFeaturesData Assay5
#' @export
AddFeaturesData.Assay5 <- function(
    object,
    features,
    ...) {
  object[[]] <- check_features_data(object, features)
  return(object)
}

check_features_data <- function(object, features) {
  features_rownames <- rownames(object)
  features_rownames_add <- rownames(features)
  if (is.null(features_rownames_add)) {
    cli::cli_abort(
      "{.arg features} must have rownames"
    )
  }
  if (length(features_rownames) != length(features_rownames_add)) {
    cli::cli_abort(
      "{.arg features} must have the same number of rownames as the object"
    )
  }
  if (any(!features_rownames_add %in% features_rownames)) {
    cli::cli_abort(
      "{.arg features} must have the same rownames as the object"
    )
  }
  features <- features[features_rownames, , drop = FALSE]
  return(features)
}
