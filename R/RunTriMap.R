#' @title Run TriMap
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunUMAP2
#' @inheritParams RunDM
#' @param n_components A number of TriMap components.
#' @param n_inliers A number of nearest neighbors for forming the nearest neighbor triplets.
#' @param n_outliers A number of outliers for forming the nearest neighbor triplets.
#' @param n_random A number of random triplets per point.
#' @param distance_method Distance metric for TriMap.
#' Options are: `"euclidean"`, `"manhattan"`, `"angular"`, `"cosine"`, `"hamming"`.
#' @param lr The learning rate for TriMap.
#' @param n_iters A number of iterations for TriMap.
#' @param apply_pca Whether to apply PCA before the nearest-neighbor calculation.
#' @param opt_method Optimization method for TriMap.
#' Options are: `"dbd"`, `"sd"`, `"momentum"`.
#' @param reduction.name Name of the reduction to be stored in the Seurat object.
#' @param reduction.key Prefix for the column names of the TriMap embeddings.
#' @param backend TriMap backend. `"cpp"` uses a compiled triplet sampler
#' and optimizer; `"python"` retains the official trimap package.
#' @param ... Passed to the trimap.TRIMAP function.
#'
#' @rdname RunTriMap
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunTriMap(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "trimap"
#' )
#' }
RunTriMap <- function(object, ...) {
  UseMethod(generic = "RunTriMap", object = object)
}

#' @rdname RunTriMap
#' @method RunTriMap Seurat
#' @export
RunTriMap.Seurat <- function(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  n_components = 2,
  n_inliers = 12,
  n_outliers = 4,
  n_random = 3,
  distance_method = "euclidean",
  lr = 0.1,
  n_iters = 400,
  apply_pca = TRUE,
  opt_method = "dbd",
  reduction.name = "trimap",
  reduction.key = "TriMap_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
) {
  if (sum(c(is.null(dims), is.null(features))) == 2) {
    log_message(
      "Please specify only one of the following arguments: dims, features",
      message_type = "error"
    )
  }
  if (!is.null(features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(
      Matrix::t(
        GetAssayData5(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(data.use) < n_components) {
      log_message(
        "Please provide as many or more features than n_components: ",
        length(features),
        " features provided, ",
        n_components,
        " TriMap components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(dims)) {
    reduction <- DefaultReduction(object, pattern = reduction)
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(dims) < n_components) {
      log_message(
        "Please provide as many or more dims than n_components: ",
        length(dims),
        " dims provided, ",
        n_components,
        " TriMap components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of dims, features",
      message_type = "error"
    )
  }
  object[[reduction.name]] <- RunTriMap(
    object = data.use,
    assay = assay,
    n_components = n_components,
    n_inliers = n_inliers,
    n_outliers = n_outliers,
    n_random = n_random,
    distance_method = distance_method,
    lr = lr,
    n_iters = n_iters,
    apply_pca = apply_pca,
    opt_method = opt_method,
    reduction.key = reduction.key,
    verbose = verbose,
    seed.use = seed.use,
    backend = backend,
    ...
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunTriMap
#' @method RunTriMap default
#' @export
RunTriMap.default <- function(
  object,
  assay = NULL,
  n_components = 2,
  n_inliers = 12,
  n_outliers = 4,
  n_random = 3,
  distance_method = "euclidean",
  lr = 0.1,
  n_iters = 400,
  apply_pca = TRUE,
  opt_method = "dbd",
  reduction.key = "TriMap_",
  verbose = TRUE,
  seed.use = 11L,
  backend = c("cpp", "python"),
  ...
) {
  backend <- match.arg(backend)
  set.seed(seed = seed.use)

  if (identical(backend, "python")) {
    PrepareEnv(modules = "trimap")
    check_python("trimap", verbose = verbose)
    trimap <- reticulate::import("trimap")

    operator <- trimap$TRIMAP(
      n_dims = as.integer(n_components),
      n_inliers = as.integer(n_inliers),
      n_outliers = as.integer(n_outliers),
      n_random = as.integer(n_random),
      distance = distance_method,
      lr = lr,
      n_iters = as.integer(n_iters),
      apply_pca = apply_pca,
      opt_method = opt_method,
      verbose = verbose,
      ...
    )
    embedding <- operator$fit_transform(object)
  } else {
    dots <- list(...)
    if (length(dots) > 0L) {
      cli::cli_warn(
        "Additional TriMap arguments are used only by {.arg backend = 'python'}."
      )
    }
    data <- native_manifold_prepare_data(object, apply_pca = apply_pca)
    n_inliers <- min(as.integer(n_inliers), nrow(data) - 1L)
    metric_id <- native_manifold_metric_id(distance_method)
    knn <- manifold_exact_knn_cpp(
      data,
      k = n_inliers + 1L,
      metric = metric_id
    )
    initial <- native_manifold_initialization(
      data,
      n_components = as.integer(n_components),
      init = "pca",
      seed = seed.use %||% 0L
    )
    embedding <- trimap_optimize_cpp(
      data = data,
      initial = initial,
      knn_index = knn[["idx"]],
      n_outliers = as.integer(n_outliers),
      n_random = as.integer(n_random),
      learning_rate = as.numeric(lr),
      iterations = as.integer(n_iters),
      optimizer = match(
        match.arg(opt_method, c("dbd", "sd", "momentum")),
        c("dbd", "sd", "momentum")
      ),
      seed = as.integer(seed.use %||% 0L),
      metric = metric_id
    )
  }
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
  if (inherits(x = object, what = "dist")) {
    rownames(x = embedding) <- attr(object, "Labels")
  } else {
    rownames(x = embedding) <- rownames(object)
  }
  reduction <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(reduction)
}
