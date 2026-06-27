#' @export
NormalizeData.Seurat <- function(
  object,
  assay = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  block.size = NULL,
  verbose = TRUE,
  ...
) {
  dots_call <- match.call(expand.dots = FALSE)$...
  if (
    !identical(normalization.method, "LogNormalize") ||
      !is.numeric(scale.factor) ||
      length(scale.factor) != 1L ||
      !is.finite(scale.factor) ||
      !is.numeric(margin) ||
      length(margin) != 1L ||
      !is.finite(margin) ||
      margin != 1 ||
      !is.null(block.size) ||
      length(dots_call) != 0L
  ) {
    stop("NormalizeData.Seurat supports LogNormalize, margin = 1, and no extra arguments.", call. = FALSE)
  }

  assay <- if (is.null(assay)) {
    SeuratObject::DefaultAssay(object)
  } else {
    if (!is.character(assay) || length(assay) != 1L || is.na(assay)) {
      stop("NormalizeData.Seurat requires a single assay name.", call. = FALSE)
    }
    assay
  }
  assays <- methods::slot(object, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' is not present in object.", assay), call. = FALSE)
  }
  assay_obj <- assays[[assay]]
  if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "StdAssay")) {
    counts <- methods::slot(assay_obj, "counts")
    if (!inherits(counts, "dgCMatrix")) {
      counts <- tryCatch(
        methods::as(counts, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(counts)) {
        stop("NormalizeData.Seurat requires counts convertible to dgCMatrix.", call. = FALSE)
      }
    }
    data_mat <- counts
    data_mat@x <- counts@x + 0
    log_normalize_dgc(data_mat, scale.factor, 100L)
    methods::slot(assay_obj, "data") <- data_mat
    assays[[assay]] <- assay_obj
    methods::slot(object, "assays") <- assays
    return(SeuratObject::LogSeuratCommand(object))
  }
  assay_slots <- methods::slotNames(assay_obj)
  if (
    !inherits(assay_obj, "StdAssay") ||
      !all(c("layers", "cells", "features") %in% assay_slots)
  ) {
    stop("NormalizeData.Seurat requires an assay with layers, cells, and features slots.", call. = FALSE)
  }
  counts_layers <- tryCatch(
    SeuratObject::Layers(assay_obj, search = "counts"),
    error = function(e) NULL
  )
  if (
    length(counts_layers) == 0L ||
      anyNA(counts_layers) ||
      !all(grepl("^counts(\\.|$)", counts_layers))
  ) {
    stop("NormalizeData.Seurat requires counts layers named 'counts' or 'counts.*'.", call. = FALSE)
  }
  layers <- methods::slot(assay_obj, "layers")
  data_layers <- sub("^counts", "data", counts_layers)
  for (i in seq_along(counts_layers)) {
    counts_layer <- counts_layers[[i]]
    data_layer <- data_layers[[i]]
    counts <- layers[[counts_layer]]
    if (is.null(counts)) {
      stop("NormalizeData.Seurat requires a counts layer.", call. = FALSE)
    }
    if (!inherits(counts, "dgCMatrix")) {
      counts <- tryCatch(
        methods::as(counts, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(counts)) {
        stop("NormalizeData.Seurat requires counts convertible to dgCMatrix.", call. = FALSE)
      }
    }

    data_mat <- counts
    data_mat@x <- counts@x + 0
    log_normalize_dgc(data_mat, scale.factor, 100L)
    layers[[data_layer]] <- data_mat
  }
  methods::slot(assay_obj, "layers") <- layers

  cm <- methods::slot(assay_obj, "cells")
  missing_data_layers <- setdiff(data_layers, colnames(cm))
  if (length(missing_data_layers) > 0L) {
    cm_data <- methods::slot(cm, ".Data")
    for (i in seq_along(counts_layers)) {
      counts_layer <- counts_layers[[i]]
      data_layer <- data_layers[[i]]
      if (!data_layer %in% missing_data_layers) {
        next
      }
      cm_data <- cbind(cm_data, cm_data[, counts_layer, drop = FALSE])
      colnames(cm_data)[ncol(cm_data)] <- data_layer
    }
    methods::slot(cm, ".Data") <- cm_data
    methods::slot(assay_obj, "cells") <- cm
  }
  fm <- methods::slot(assay_obj, "features")
  missing_data_layers <- setdiff(data_layers, colnames(fm))
  if (length(missing_data_layers) > 0L) {
    fm_data <- methods::slot(fm, ".Data")
    for (i in seq_along(counts_layers)) {
      counts_layer <- counts_layers[[i]]
      data_layer <- data_layers[[i]]
      if (!data_layer %in% missing_data_layers) {
        next
      }
      fm_data <- cbind(fm_data, fm_data[, counts_layer, drop = FALSE])
      colnames(fm_data)[ncol(fm_data)] <- data_layer
    }
    methods::slot(fm, ".Data") <- fm_data
    methods::slot(assay_obj, "features") <- fm
  }
  assays[[assay]] <- assay_obj
  methods::slot(object, "assays") <- assays
  SeuratObject::LogSeuratCommand(object)
}

#' Normalize a single-cell object
#'
#' @param object Object to normalize.
#' @param ... Passed to methods.
#'
#' @return A normalized object of the same class.
#' @export
NormalizeData <- function(object, ...) {
  UseMethod("NormalizeData")
}

#' @export
NormalizeData.default <- function(object, ...) {
  stop("NormalizeData supports Seurat objects.", call. = FALSE)
}
