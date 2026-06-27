#' @export
ScaleData.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  split.by = NULL,
  model.use = "linear",
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  if (
    !is.null(vars.to.regress) ||
      !is.null(split.by) ||
      !identical(model.use, "linear") ||
      !isFALSE(use.umi) ||
      !isTRUE(do.scale) ||
      !isTRUE(do.center)
  ) {
    stop("ScaleData.Seurat supports linear scaling without regression, split.by, or use.umi.", call. = FALSE)
  }
  assay_name <- if (is.null(assay)) {
    SeuratObject::DefaultAssay(object)
  } else {
    assay[1L]
  }
  assay_obj <- methods::slot(object, "assays")[[assay_name]]
  if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "StdAssay")) {
    data_mat <- methods::slot(assay_obj, "data")
    if (!inherits(data_mat, "dgCMatrix")) {
      data_mat <- tryCatch(
        methods::as(data_mat, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(data_mat)) {
        stop("ScaleData.Seurat requires data convertible to dgCMatrix.", call. = FALSE)
      }
    }
    features <- if (is.null(features)) {
      SeuratObject::VariableFeatures(object)
    } else {
      features
    }
    if (length(features) == 0) {
      features <- rownames(assay_obj)
    }

    all_genes <- rownames(assay_obj)
    features <- intersect(features, all_genes)
    features <- features[order(match(features, all_genes))]
    idx <- match(features, all_genes) - 1L

    result <- scale_sparse_full(data_mat, idx, scale.max)
    dimnames(result) <- list(features, colnames(object))
    methods::slot(assay_obj, "scale.data") <- result
    methods::slot(object, "assays")[[assay_name]] <- assay_obj
    return(object)
  }
  if (
    !inherits(assay_obj, "Assay5") ||
      length(SeuratObject::Layers(assay_obj, search = "data")) == 0L
  ) {
    stop("ScaleData.Seurat requires an Assay5 object with a data layer.", call. = FALSE)
  }
  data_layers <- SeuratObject::Layers(assay_obj, search = "data")
  data_layers <- data_layers[grepl("^data(\\.|$)", data_layers)]
  if (length(data_layers) == 0L) {
    stop("ScaleData.Seurat requires an Assay5 object with a data layer.", call. = FALSE)
  }
  layers <- methods::slot(assay_obj, "layers")
  data_mat <- if (length(data_layers) == 1L) {
    layers[[data_layers]]
  } else {
    cells_map <- methods::slot(assay_obj, "cells")
    data_mats <- lapply(data_layers, function(data_layer) {
      mat <- layers[[data_layer]]
      if (is.null(colnames(mat))) {
        colnames(mat) <- rownames(cells_map)[cells_map[, data_layer, drop = TRUE]]
      }
      mat
    })
    do.call(cbind, unname(data_mats))
  }
  if (!inherits(data_mat, "dgCMatrix")) {
    data_mat <- tryCatch(
      methods::as(data_mat, "dgCMatrix"),
      error = function(e) NULL
    )
    if (is.null(data_mat)) {
      stop("ScaleData.Seurat requires data convertible to dgCMatrix.", call. = FALSE)
    }
  }
  if (length(data_layers) > 1L) {
    data_mat <- data_mat[, colnames(object), drop = FALSE]
  }
  features <- if (is.null(features)) {
    SeuratObject::VariableFeatures(object)
  } else {
    features
  }
  if (length(features) == 0) {
    features <- rownames(assay_obj)
  }

  all_genes <- rownames(assay_obj)
  features <- intersect(features, all_genes)
  features <- features[order(match(features, all_genes))]
  idx <- match(features, all_genes) - 1L

  result <- scale_sparse_full(data_mat, idx, scale.max)
  dimnames(result) <- list(features, colnames(object))

  assay_obj@layers[["scale.data"]] <- result
  if (!"scale.data" %in% colnames(assay_obj@cells)) {
    cm <- assay_obj@cells
    methods::slot(assay_obj@cells, ".Data") <- cbind(
      cm,
      matrix(
        TRUE,
        nrow = nrow(cm),
        ncol = 1,
        dimnames = list(rownames(cm), "scale.data")
      )
    )
    fm <- assay_obj@features
    methods::slot(assay_obj@features, ".Data") <- cbind(
      fm,
      matrix(
        rownames(fm) %in% features,
        nrow = nrow(fm),
        ncol = 1,
        dimnames = list(rownames(fm), "scale.data")
      )
    )
  }
  methods::slot(object, "assays")[[assay_name]] <- assay_obj
  object
}

#' Scale expression data
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return The input object with scaled data.
#' @export
ScaleData <- function(object, ...) {
  UseMethod("ScaleData")
}

#' @export
ScaleData.default <- function(object, ...) {
  stop("ScaleData supports Seurat objects.", call. = FALSE)
}
