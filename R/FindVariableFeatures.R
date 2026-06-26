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
    stop("FindVariableFeatures.StdAssay supports selection.method = 'vst'.", call. = FALSE)
  }
  if (
    !inherits(object, "Assay5") ||
      length(SeuratObject::Layers(object, search = "counts")) == 0L
  ) {
    stop("FindVariableFeatures.StdAssay requires an Assay5 object with a counts layer.", call. = FALSE)
  }
  counts_layers <- SeuratObject::Layers(object, search = "counts")
  counts_layers <- counts_layers[grepl("^counts(\\.|$)", counts_layers)]
  if (length(counts_layers) == 0L) {
    stop("FindVariableFeatures.StdAssay requires an Assay5 object with a counts layer.", call. = FALSE)
  }
  layers <- methods::slot(object, "layers")
  data <- if (length(counts_layers) == 1L) {
    layers[[counts_layers]]
  } else {
    do.call(cbind, unname(layers[counts_layers]))
  }
  if (!inherits(data, "dgCMatrix")) {
    data <- tryCatch(
      methods::as(data, "dgCMatrix"),
      error = function(e) NULL
    )
    if (is.null(data)) {
      stop("FindVariableFeatures.StdAssay requires counts convertible to dgCMatrix.", call. = FALSE)
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
  rownames(hvf.info) <- SeuratObject::Features(object, layer = "counts")
  object[["var.features"]] <- NULL
  object[["var.features.rank"]] <- NULL
  object[[names(hvf.info)]] <- NULL
  object[[names(hvf.info)]] <- hvf.info
  SeuratObject::VariableFeatures(object) <- rownames(hvf.info)[
    order(hvf.info$vf_vst_counts_rank, na.last = NA)
  ]
  object
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
  if (!identical(selection.method, "vst")) {
    stop("FindVariableFeatures.Seurat supports selection.method = 'vst'.", call. = FALSE)
  }
  assay <- if (is.null(assay)) {
    SeuratObject::DefaultAssay(object)
  } else {
    assay[1L]
  }
  assay_obj <- object[[assay]]
  if (
    !inherits(assay_obj, "Assay5") ||
      length(SeuratObject::Layers(assay_obj, search = "counts")) == 0L
  ) {
    stop("FindVariableFeatures.Seurat requires an Assay5 object with a counts layer.", call. = FALSE)
  }
  assay_obj <- FindVariableFeatures.StdAssay(
    object = assay_obj,
    nfeatures = nfeatures,
    span = loess.span,
    clip = if (identical(clip.max, "auto")) NULL else clip.max,
    verbose = verbose
  )
  methods::slot(object, "assays")[[assay]] <- assay_obj
  object
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
FindVariableFeatures.default <- function(object, ...) {
  stop("FindVariableFeatures supports Seurat and StdAssay objects.", call. = FALSE)
}
