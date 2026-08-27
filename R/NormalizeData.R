sct_norm_sparse <- function(counts, normalization.method, scale.factor, margin) {
  if (identical(normalization.method, "RC")) {
    counts@x <- counts@x + 0
    sizes <- Matrix::colSums(counts)
    sizes[sizes == 0] <- 1
    counts@x <- counts@x * rep.int(
      scale.factor / as.numeric(sizes),
      diff(counts@p)
    )
    return(counts)
  }
  if (identical(normalization.method, "CLR")) {
    if (as.integer(margin) == 2L) {
      denom <- exp(Matrix::colSums(log1p(counts)) / nrow(counts))
      counts@x <- log1p(counts@x / rep.int(denom, diff(counts@p)))
    } else {
      denom <- exp(Matrix::rowSums(log1p(counts)) / ncol(counts))
      counts@x <- log1p(counts@x / denom[counts@i + 1L])
    }
    return(counts)
  }
  norm <- counts
  norm@x <- counts@x + 0
  log_normalize_dgc(norm, scale.factor, 100L)
  norm
}

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
    !normalization.method %in% c("LogNormalize", "CLR", "RC") ||
      !is.numeric(scale.factor) ||
      length(scale.factor) != 1L ||
      !is.finite(scale.factor) ||
      scale.factor <= 0 ||
      !is.numeric(margin) ||
      length(margin) != 1L ||
      !is.finite(margin) ||
      !(margin %in% c(1, 2)) ||
      (identical(normalization.method, "LogNormalize") &&
        as.integer(margin) != 1L)
  ) {
    log_message(
      "NormalizeData.Seurat supports LogNormalize (margin = 1), CLR, and RC with margin 1 or 2.",
      message_type = "error"
    )
  }
  if (!is.null(block.size)) {
    message("NormalizeData: block.size is ignored; normalization runs in a single pass.")
  }

  assay <- if (is.null(assay)) {
    SeuratObject::DefaultAssay(object)
  } else {
    if (!is.character(assay) || length(assay) != 1L || is.na(assay)) {
      log_message("NormalizeData.Seurat requires a single assay name.", message_type = "error")
    }
    assay
  }
  assays <- methods::slot(object, "assays")
  if (!assay %in% names(assays)) {
    log_message(sprintf("Assay '%s' is not present in object.", assay), message_type = "error")
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
        log_message("NormalizeData.Seurat requires counts convertible to dgCMatrix.", message_type = "error")
      }
    }
    data_mat <- sct_norm_sparse(counts, normalization.method, scale.factor, margin)
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
    log_message("NormalizeData.Seurat requires an assay with layers, cells, and features slots.", message_type = "error")
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
    log_message("NormalizeData.Seurat requires counts layers named 'counts' or 'counts.*'.", message_type = "error")
  }
  layers <- methods::slot(assay_obj, "layers")
  data_layers <- sub("^counts", "data", counts_layers)
  for (i in seq_along(counts_layers)) {
    counts_layer <- counts_layers[[i]]
    data_layer <- data_layers[[i]]
    counts <- layers[[counts_layer]]
    if (is.null(counts)) {
      log_message("NormalizeData.Seurat requires a counts layer.", message_type = "error")
    }
    if (!inherits(counts, "dgCMatrix")) {
      counts <- tryCatch(
        methods::as(counts, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(counts)) {
        log_message("NormalizeData.Seurat requires counts convertible to dgCMatrix.", message_type = "error")
      }
    }

    data_mat <- sct_norm_sparse(counts, normalization.method, scale.factor, margin)
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
  log_message("NormalizeData supports Seurat objects.", message_type = "error")
}
