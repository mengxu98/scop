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
  dots <- list(...)
  delegate <- function() {
    do.call(
      utils::getFromNamespace("NormalizeData.Seurat", "Seurat"),
      c(
        list(
          object = object,
          assay = assay,
          normalization.method = normalization.method,
          scale.factor = scale.factor,
          margin = margin,
          block.size = block.size,
          verbose = verbose
        ),
        dots
      )
    )
  }
  if (length(dots) > 0L && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    return(delegate())
  }
  if (
    !is.character(normalization.method) ||
      length(normalization.method) != 1L || is.na(normalization.method) ||
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
    return(delegate())
  }
  if (!is.null(block.size)) {
    return(delegate())
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
  seurat_fallback <- delegate
  layer_arg <- dots[["layer"]] %||% "counts"
  save_arg <- dots[["save"]] %||% "data"
  remaining_dots <- dots[setdiff(names(dots), c("layer", "save"))]
  if (length(remaining_dots) > 0L || inherits(assay_obj, "SCTAssay")) {
    return(seurat_fallback())
  }
  if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "StdAssay")) {
    if (!identical(layer_arg, "counts") || !identical(save_arg, "data")) {
      return(seurat_fallback())
    }
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
    object[[assay]] <- assay_obj
    return(SeuratObject::LogSeuratCommand(object))
  }
  assay_slots <- methods::slotNames(assay_obj)
  if (
    !inherits(assay_obj, "StdAssay") ||
      !all(c("layers", "cells", "features") %in% assay_slots)
  ) {
    log_message("NormalizeData.Seurat requires an assay with layers, cells, and features slots.", message_type = "error")
  }
  source_patterns <- unique(layer_arg)
  source_layers <- tryCatch(
    SeuratObject::Layers(assay_obj, search = source_patterns),
    error = function(e) NULL
  )
  if (length(source_layers) == 0L || anyNA(source_layers)) {
    return(seurat_fallback())
  }
  save_layers <- save_arg
  if (length(save_layers) != length(source_layers)) {
    save_layers <- make.unique(gsub(
      pattern = source_patterns,
      replacement = save_arg,
      x = source_layers
    ))
  }
  layers <- methods::slot(assay_obj, "layers")
  for (i in seq_along(source_layers)) {
    source_layer <- source_layers[[i]]
    save_layer <- save_layers[[i]]
    counts <- layers[[source_layer]]
    if (is.null(counts)) {
      return(seurat_fallback())
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
    layers[[save_layer]] <- data_mat
  }
  methods::slot(assay_obj, "layers") <- layers

  cm <- methods::slot(assay_obj, "cells")
  missing_data_layers <- setdiff(save_layers, colnames(cm))
  if (length(missing_data_layers) > 0L) {
    cm_data <- methods::slot(cm, ".Data")
    for (i in seq_along(source_layers)) {
      source_layer <- source_layers[[i]]
      save_layer <- save_layers[[i]]
      if (!save_layer %in% missing_data_layers) {
        next
      }
      cm_data <- cbind(cm_data, cm_data[, source_layer, drop = FALSE])
      colnames(cm_data)[ncol(cm_data)] <- save_layer
    }
    methods::slot(cm, ".Data") <- cm_data
    methods::slot(assay_obj, "cells") <- cm
  }
  fm <- methods::slot(assay_obj, "features")
  missing_data_layers <- setdiff(save_layers, colnames(fm))
  if (length(missing_data_layers) > 0L) {
    fm_data <- methods::slot(fm, ".Data")
    for (i in seq_along(source_layers)) {
      source_layer <- source_layers[[i]]
      save_layer <- save_layers[[i]]
      if (!save_layer %in% missing_data_layers) {
        next
      }
      fm_data <- cbind(fm_data, fm_data[, source_layer, drop = FALSE])
      colnames(fm_data)[ncol(fm_data)] <- save_layer
    }
    methods::slot(fm, ".Data") <- fm_data
    methods::slot(assay_obj, "features") <- fm
  }
  object[[assay]] <- assay_obj
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
NormalizeData.default <- function(
  object,
  normalization.method = c("LogNormalize", "CLR", "RC"),
  scale.factor = 10000,
  cmargin = 2L,
  margin = 1L,
  verbose = TRUE,
  ...
) {
  method <- normalization.method[[1L]]
  native <- inherits(object, "dgCMatrix") &&
    method %in% c("LogNormalize", "CLR", "RC") &&
    is.numeric(scale.factor) && length(scale.factor) == 1L &&
    is.finite(scale.factor) && scale.factor > 0 &&
    is.numeric(margin) && length(margin) == 1L && margin %in% c(1, 2) &&
    (!identical(method, "LogNormalize") || identical(as.integer(cmargin), 2L))
  if (native) {
    return(sct_norm_sparse(object, method, scale.factor, margin))
  }
  do.call(
    utils::getFromNamespace("NormalizeData.default", "Seurat"),
    c(
      list(
        object = object,
        normalization.method = normalization.method,
        scale.factor = scale.factor,
        cmargin = cmargin,
        margin = margin,
        verbose = verbose
      ),
      list(...)
    )
  )
}
