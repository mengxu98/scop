variable_features_vst_from_stats <- function(
  mu,
  variance,
  nnz_per_row,
  n_cells,
  var_std,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  ...
) {
  nfeatures <- length(mu)
  var_expected <- numeric(nfeatures)
  not.const <- variance > 0
  log_mean_nc <- log10(mu[not.const])
  log_var_nc <- log10(variance[not.const])
  simple_loess <- utils::getFromNamespace("simpleLoess", "stats")
  if (length(log_mean_nc) > 0L) {
    fit_result <- simple_loess(
      y = log_var_nc,
      x = matrix(log_mean_nc, ncol = 1),
      weights = rep.int(1, length(log_mean_nc)),
      span = span,
      degree = 2L,
      parametric = FALSE,
      drop.square = FALSE,
      normalize = FALSE,
      statistics = "approximate",
      surface = "interpolate",
      cell = 0.2,
      iterations = 1L,
      iterTrace = FALSE,
      trace.hat = "approximate"
    )
    var_expected[not.const] <- 10^fit_result$fitted
  }

  sd_vec <- sqrt(var_expected)
  vmax <- if (is.null(clip)) sqrt(n_cells) else clip
  if (is.null(var_std)) {
    log_message("variable_features_vst_from_stats requires standardized variance.", message_type = "error")
  }

  hvf.info <- SeuratObject::EmptyDF(n = nfeatures)
  hvf.info$mean <- mu
  hvf.info$variance <- variance
  hvf.info$variance.expected <- var_expected
  hvf.info$variance.standardized <- var_std
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA
  vf <- head(order(var_std, decreasing = TRUE), n = nselect)
  hvf.info$variable[vf] <- TRUE
  hvf.info$rank[vf] <- seq_along(vf)
  hvf.info
}

variable_features_vst_sparse <- function(
  data,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  ...
) {
  nfeatures <- nrow(data)
  n_cells <- ncol(data)
  mv <- sparse_row_mean_var(
    p = data@p,
    i = data@i,
    x = data@x,
    nrow = nfeatures,
    ncol = n_cells
  )
  mu <- mv$mean
  variance <- mv$variance
  nnz_per_row <- mv$nnz

  var_expected <- numeric(nfeatures)
  not.const <- variance > 0
  log_mean_nc <- log10(mu[not.const])
  log_var_nc <- log10(variance[not.const])
  if (length(log_mean_nc) > 0L) {
    simple_loess <- utils::getFromNamespace("simpleLoess", "stats")
    fit_result <- simple_loess(
      y = log_var_nc,
      x = matrix(log_mean_nc, ncol = 1),
      weights = rep.int(1, length(log_mean_nc)),
      span = span,
      degree = 2L,
      parametric = FALSE,
      drop.square = FALSE,
      normalize = FALSE,
      statistics = "approximate",
      surface = "interpolate",
      cell = 0.2,
      iterations = 1L,
      iterTrace = FALSE,
      trace.hat = "approximate"
    )
    var_expected[not.const] <- 10^fit_result$fitted
  }
  sd_vec <- sqrt(var_expected)
  vmax <- if (is.null(clip)) sqrt(n_cells) else clip
  var_std <- sparse_row_var_std(
    p = data@p,
    i = data@i,
    x = data@x,
    nrow = nfeatures,
    ncol = n_cells,
    mu = mu,
    sd = sd_vec,
    vmax = vmax,
    nnzPerRow = nnz_per_row
  )

  variable_features_vst_from_stats(
    mu = mu,
    variance = variance,
    nnz_per_row = nnz_per_row,
    n_cells = n_cells,
    var_std = var_std,
    nselect = nselect,
    span = span,
    clip = clip,
    verbose = verbose,
    ...
  )
}

variable_features_vst_sparse_layers <- function(
  layers,
  nselect = 2000L,
  span = 0.3,
  clip = NULL,
  verbose = TRUE,
  ...
) {
  if (!length(layers)) {
    log_message("No sparse layers were supplied.", message_type = "error")
  }
  nfeatures <- nrow(layers[[1L]])
  feature_names <- rownames(layers[[1L]])
  for (layer in layers) {
    if (!inherits(layer, "dgCMatrix")) {
      log_message("variable_features_vst_sparse_layers requires dgCMatrix layers.", message_type = "error")
    }
    if (nrow(layer) != nfeatures || !identical(rownames(layer), feature_names)) {
      log_message("All sparse layers must have identical feature rows.", message_type = "error")
    }
  }
  mv <- sparse_row_mean_var_dgc_list(layers, nfeatures)
  mu <- mv$mean
  variance <- mv$variance
  nnz_per_row <- mv$nnz
  n_cells <- mv$ncol

  var_expected <- numeric(nfeatures)
  not.const <- variance > 0
  log_mean_nc <- log10(mu[not.const])
  log_var_nc <- log10(variance[not.const])
  if (length(log_mean_nc) > 0L) {
    simple_loess <- utils::getFromNamespace("simpleLoess", "stats")
    fit_result <- simple_loess(
      y = log_var_nc,
      x = matrix(log_mean_nc, ncol = 1),
      weights = rep.int(1, length(log_mean_nc)),
      span = span,
      degree = 2L,
      parametric = FALSE,
      drop.square = FALSE,
      normalize = FALSE,
      statistics = "approximate",
      surface = "interpolate",
      cell = 0.2,
      iterations = 1L,
      iterTrace = FALSE,
      trace.hat = "approximate"
    )
    var_expected[not.const] <- 10^fit_result$fitted
  }
  sd_vec <- sqrt(var_expected)
  vmax <- if (is.null(clip)) sqrt(n_cells) else clip
  var_std <- sparse_row_var_std_dgc_list(
    mats = layers,
    nrow = nfeatures,
    mu = mu,
    sd = sd_vec,
    vmax = vmax
  )

  variable_features_vst_from_stats(
    mu = mu,
    variance = variance,
    nnz_per_row = nnz_per_row,
    n_cells = n_cells,
    var_std = var_std,
    nselect = nselect,
    span = span,
    clip = clip,
    verbose = verbose,
    ...
  )
}

#' @export
FindVariableFeatures.StdAssay <- function(
  object,
  method = NULL,
  nfeatures = 2000L,
  layer = NULL,
  span = 0.3,
  clip = NULL,
  key = NULL,
  verbose = TRUE,
  selection.method = "vst",
  ...
) {
  if (!identical(selection.method, "vst")) {
    log_message("FindVariableFeatures.StdAssay supports selection.method = 'vst'.", message_type = "error")
  }
  if (inherits(object, "Assay") && !inherits(object, "StdAssay")) {
    data <- methods::slot(object, "counts")
    if (!inherits(data, "dgCMatrix")) {
      data <- tryCatch(
        methods::as(data, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(data)) {
        log_message("FindVariableFeatures.StdAssay requires counts convertible to dgCMatrix.", message_type = "error")
      }
    }
    hvf.info <- variable_features_vst_sparse(
      data,
      nselect = nfeatures,
      span = span,
      clip = clip,
      verbose = verbose
    )
    colnames(hvf.info) <- paste("vf_vst_counts", colnames(hvf.info), sep = "_")
    rownames(hvf.info) <- rownames(object)
    meta_features <- methods::slot(object, "meta.features")
    meta_features <- meta_features[
      ,
      setdiff(colnames(meta_features), c("var.features", "var.features.rank", names(hvf.info))),
      drop = FALSE
    ]
    methods::slot(object, "meta.features") <- cbind(meta_features, hvf.info)
    SeuratObject::VariableFeatures(object) <- rownames(hvf.info)[
      order(hvf.info$vf_vst_counts_rank, na.last = NA)
    ]
    return(object)
  }
  layer_query <- layer %||% "counts"
  if (
    !inherits(object, "Assay5") ||
      length(SeuratObject::Layers(object, search = layer_query)) == 0L
  ) {
    log_message("FindVariableFeatures.StdAssay requires an Assay5 object with a counts layer.", message_type = "error")
  }
  counts_layers <- SeuratObject::Layers(object, search = layer_query)
  if (length(counts_layers) == 0L) {
    log_message("FindVariableFeatures.StdAssay requires an Assay5 object with a counts layer.", message_type = "error")
  }
  layers <- methods::slot(object, "layers")
  if (length(counts_layers) == 1L) {
    data <- layers[[counts_layers]]
    if (!inherits(data, "dgCMatrix")) {
      data <- tryCatch(
        methods::as(data, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(data)) {
        log_message("FindVariableFeatures.StdAssay requires counts convertible to dgCMatrix.", message_type = "error")
      }
    }
    hvf.info <- variable_features_vst_sparse(
      data,
      nselect = nfeatures,
      span = span,
      clip = clip,
      verbose = verbose
    )
  } else {
    data_layers <- lapply(layers[counts_layers], function(data) {
      if (inherits(data, "dgCMatrix")) {
        return(data)
      }
      out <- tryCatch(
        methods::as(data, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(out)) {
        log_message("FindVariableFeatures.StdAssay requires counts convertible to dgCMatrix.", message_type = "error")
      }
      out
    })
    feature_table <- methods::slot(object, "features")
    gene_names <- rownames(feature_table)
    layer_hvf <- vector("list", length(counts_layers))
    for (li in seq_along(counts_layers)) {
      lname <- counts_layers[[li]]
      ldata <- data_layers[[li]]
      hvf_layer <- variable_features_vst_sparse(
        ldata,
        nselect = nfeatures,
        span = span,
        clip = clip,
        verbose = verbose
      )
      layer_features <- feature_table[, lname]
      if (is.null(layer_features)) {
        layer_features <- gene_names
      } else if (is.logical(layer_features)) {
        layer_features <- gene_names[layer_features]
      }
      rownames(hvf_layer) <- layer_features
      hvf_layer <- hvf_layer[match(gene_names, layer_features), , drop = FALSE]
      rownames(hvf_layer) <- gene_names
      colnames(hvf_layer) <- paste("vf_vst", lname, colnames(hvf_layer), sep = "_")
      layer_hvf[[li]] <- hvf_layer
    }
    object[["var.features"]] <- NULL
    object[["var.features.rank"]] <- NULL
    meta_data <- methods::slot(object, "meta.data")
    new_cols <- unlist(lapply(layer_hvf, colnames), use.names = FALSE)
    hvf_columns <- do.call(c, lapply(layer_hvf, as.list))
    names(hvf_columns) <- new_cols
    meta_data <- meta_data[, setdiff(colnames(meta_data), new_cols), drop = FALSE]
    meta_data[new_cols] <- hvf_columns
    methods::slot(object, "meta.data") <- meta_data
    consensus_features <- SeuratObject::VariableFeatures(
      object,
      nfeatures = nfeatures,
      method = "vst"
    )
    SeuratObject::VariableFeatures(object) <- consensus_features
    return(object)
  }
  colnames(hvf.info) <- paste("vf_vst", counts_layers[[1L]], colnames(hvf.info), sep = "_")
  rownames(hvf.info) <- SeuratObject::Features(object, layer = counts_layers[[1L]])
  object[["var.features"]] <- NULL
  object[["var.features.rank"]] <- NULL
  object[[names(hvf.info)]] <- NULL
  object[[names(hvf.info)]] <- hvf.info
  SeuratObject::VariableFeatures(object) <- rownames(hvf.info)[
    order(hvf.info$vf_vst_counts_rank, na.last = NA)
  ]
  object
}

sct_mvp_info <- function(data_mat,
                         selection.method,
                         mean.function,
                         dispersion.function,
                         num.bin,
                         binning.method,
                         nfeatures,
                         mean.cutoff,
                         dispersion.cutoff,
                         verbose) {
  mean_fun <- mean.function %||% utils::getFromNamespace("FastExpMean", "Seurat")
  disp_fun <- dispersion.function %||% utils::getFromNamespace("FastLogVMR", "Seurat")
  feature.mean <- mean_fun(data_mat, verbose)
  feature.dispersion <- disp_fun(data_mat, verbose)
  names(feature.mean) <- names(feature.dispersion) <- rownames(data_mat)
  feature.dispersion[is.na(feature.dispersion)] <- 0
  feature.mean[is.na(feature.mean)] <- 0
  data.x.breaks <- switch(
    EXPR = binning.method,
    equal_width = num.bin,
    equal_frequency = c(-1, stats::quantile(
      x = feature.mean[feature.mean > 0],
      probs = seq.int(from = 0, to = 1, length.out = num.bin)
    )),
    stop("Unknown binning method: ", binning.method)
  )
  data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
  mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
  sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = stats::sd)
  dispersion.scaled <-
    (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
  hvf.info <- data.frame(
    mvp.mean = feature.mean,
    mvp.dispersion = feature.dispersion,
    mvp.dispersion.scaled = dispersion.scaled,
    row.names = rownames(data_mat)
  )
  hvf.info$variable <- FALSE
  hvf.info$rank <- NA_integer_
  if (identical(selection.method, "dispersion")) {
    selected <- head(
      order(hvf.info$mvp.dispersion, decreasing = TRUE),
      nfeatures
    )
  } else {
    ordered <- order(hvf.info$mvp.dispersion, decreasing = TRUE)
    keep <- hvf.info$mvp.mean[ordered] > mean.cutoff[1] &
      hvf.info$mvp.mean[ordered] < mean.cutoff[2] &
      hvf.info$mvp.dispersion.scaled[ordered] > dispersion.cutoff[1] &
      hvf.info$mvp.dispersion.scaled[ordered] < dispersion.cutoff[2]
    selected <- ordered[which(keep)]
  }
  hvf.info$variable[selected] <- TRUE
  hvf.info$rank[selected] <- seq_along(selected)
  hvf.info
}

#' @export
FindVariableFeatures.Seurat <- function(
  object,
  assay = NULL,
  selection.method = "vst",
  loess.span = 0.3,
  clip.max = "auto",
  mean.function = NULL,
  dispersion.function = NULL,
  num.bin = 20,
  binning.method = "equal_width",
  nfeatures = 2000,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf),
  verbose = TRUE,
  ...
) {
  dots <- list(...)
  unnamed_dots <- length(dots) > 0L &&
    (is.null(names(dots)) || any(!nzchar(names(dots))))
  selection.method <- switch(
    EXPR = selection.method,
    mvp = "mean.var.plot",
    disp = "dispersion",
    selection.method
  )
  if (!selection.method %in% c("vst", "dispersion", "mean.var.plot")) {
    log_message(
      "{.fn FindVariableFeatures} supports {.val vst}, {.val disp}, and {.val mean.var.plot}.",
      message_type = "error"
    )
  }
  assay <- if (is.null(assay)) {
    SeuratObject::DefaultAssay(object)
  } else {
    assay[1L]
  }
  assay_obj <- object[[assay]]
  # Custom layer/method arguments and SCTAssay residual handling are owned by
  # Seurat. Delegate these uncommon branches instead of silently dropping them.
  native_layer <- dots[["layer"]]
  remaining_dots <- dots[setdiff(names(dots), "layer")]
  if (
    unnamed_dots || length(remaining_dots) > 0L ||
      inherits(assay_obj, "SCTAssay") ||
      (!is.null(native_layer) && !identical(selection.method, "vst"))
  ) {
    fallback_args <- list(
      object = object,
      assay = assay,
      selection.method = selection.method,
      loess.span = loess.span,
      clip.max = clip.max,
      num.bin = num.bin,
      binning.method = binning.method,
      nfeatures = nfeatures,
      mean.cutoff = mean.cutoff,
      dispersion.cutoff = dispersion.cutoff,
      verbose = verbose
    )
    if (!is.null(mean.function)) {
      fallback_args$mean.function <- mean.function
    }
    if (!is.null(dispersion.function)) {
      fallback_args$dispersion.function <- dispersion.function
    }
    return(do.call(
      utils::getFromNamespace("FindVariableFeatures.Seurat", "Seurat"),
      c(fallback_args, dots)
    ))
  }
  if (identical(selection.method, "vst")) {
    if (
      length(SeuratObject::Layers(assay_obj, search = "counts")) == 0L
    ) {
      log_message("FindVariableFeatures.Seurat requires an assay with a counts layer.", message_type = "error")
    }
    assay_obj <- FindVariableFeatures.StdAssay(
      object = assay_obj,
      nfeatures = nfeatures,
      layer = native_layer,
      span = loess.span,
      clip = if (identical(clip.max, "auto")) NULL else clip.max,
      verbose = verbose
    )
    methods::slot(object, "assays")[[assay]] <- assay_obj
    return(SeuratObject::LogSeuratCommand(object))
  }
  data_layers <- SeuratObject::Layers(assay_obj, search = "data")
  data_layers <- data_layers[grepl("^data(\\.|$)", data_layers)]
  if (length(data_layers) == 0L) {
    log_message("FindVariableFeatures.Seurat requires an assay with a data layer.", message_type = "error")
  }

  if (inherits(assay_obj, "StdAssay") && length(data_layers) > 1L) {
    gene_names <- rownames(assay_obj)
    key <- if (identical(selection.method, "dispersion")) "disp" else "mvp"
    layer_hvf <- lapply(data_layers, function(data_layer) {
      data_mat <- SeuratObject::LayerData(assay_obj, layer = data_layer)
      hvf.info <- sct_mvp_info(
        data_mat = data_mat,
        selection.method = selection.method,
        mean.function = mean.function,
        dispersion.function = dispersion.function,
        num.bin = num.bin,
        binning.method = binning.method,
        nfeatures = nfeatures,
        mean.cutoff = mean.cutoff,
        dispersion.cutoff = dispersion.cutoff,
        verbose = verbose
      )
      layer_features <- rownames(data_mat)
      rownames(hvf.info) <- layer_features
      hvf.info <- hvf.info[match(gene_names, layer_features), , drop = FALSE]
      rownames(hvf.info) <- gene_names
      colnames(hvf.info) <- paste(
        "vf",
        key,
        data_layer,
        colnames(hvf.info),
        sep = "_"
      )
      hvf.info
    })
    meta_data <- methods::slot(assay_obj, "meta.data")
    new_cols <- unlist(lapply(layer_hvf, colnames), use.names = FALSE)
    hvf_columns <- do.call(c, lapply(layer_hvf, as.list))
    names(hvf_columns) <- new_cols
    meta_data <- meta_data[, setdiff(colnames(meta_data), new_cols), drop = FALSE]
    meta_data[new_cols] <- hvf_columns
    methods::slot(assay_obj, "meta.data") <- meta_data
    consensus_features <- SeuratObject::VariableFeatures(
      assay_obj,
      nfeatures = nfeatures,
      method = key
    )
    SeuratObject::VariableFeatures(assay_obj) <- consensus_features
    methods::slot(object, "assays")[[assay]] <- assay_obj
    return(SeuratObject::LogSeuratCommand(object))
  }

  data_mat <- tryCatch(
    SeuratObject::GetAssayData(assay_obj, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(data_mat)) {
    log_message("FindVariableFeatures.Seurat requires an assay with a data layer.", message_type = "error")
  }
  hvf.info <- sct_mvp_info(
    data_mat = data_mat,
    selection.method = selection.method,
    mean.function = mean.function,
    dispersion.function = dispersion.function,
    num.bin = num.bin,
    binning.method = binning.method,
    nfeatures = nfeatures,
    mean.cutoff = mean.cutoff,
    dispersion.cutoff = dispersion.cutoff,
    verbose = verbose
  )
  top.features <- rownames(hvf.info)[order(hvf.info$rank, na.last = NA)]
  SeuratObject::VariableFeatures(assay_obj) <- top.features
  vf.name <- "mvp.variable"
  hvf.info <- hvf.info[, 1:3, drop = FALSE]
  assay_obj[[colnames(hvf.info)]] <- hvf.info
  assay_obj[[vf.name]] <- rownames(assay_obj[[]]) %in% top.features
  methods::slot(object, "assays")[[assay]] <- assay_obj
  SeuratObject::LogSeuratCommand(object)
}

#' Find variable features
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return The input object with variable features recorded.
#' @export
FindVariableFeatures <- function(object, ...) {
  UseMethod("FindVariableFeatures")
}

#' @export
FindVariableFeatures.default <- function(
  object,
  method = NULL,
  nfeatures = 2000L,
  verbose = TRUE,
  selection.method = selection.method,
  ...
) {
  vst <- utils::getFromNamespace("VST", "Seurat")
  method_use <- method %||% vst
  get("FindVariableFeatures.default", asNamespace("Seurat"))(
    object = object,
    method = method_use,
    nfeatures = nfeatures,
    verbose = verbose,
    ...
  )
}
