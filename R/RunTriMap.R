#' Run TriMap (Large-scale Dimensionality Reduction Using Triplets)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer to be used. Default is "data".
#' @param n_components An integer specifying the number of TriMap components. Default is 2.
#' @param n_inliers An integer specifying the number of nearest neighbors for forming the nearest neighbor triplets. Default is 12.
#' @param n_outliers An integer specifying the number of outliers for forming the nearest neighbor triplets. Default is 4.
#' @param n_random An integer specifying the number of random triplets per point. Default is 3.
#' @param distance_method A character string specifying the distance metric for TriMap. Options are: "euclidean", "manhattan", "angular", "cosine", "hamming". Default is "euclidean".
#' @param lr A numeric value specifying the learning rate for TriMap. Default is 0.1.
#' @param n_iters An integer specifying the number of iterations for TriMap. Default is 400.
#' @param apply_pca A logical value indicating whether to apply PCA before the nearest-neighbor calculation. Default is TRUE.
#' @param opt_method A character string specifying the optimization method for TriMap. Options are: "dbd", "sd", "momentum". Default is "dbd".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "trimap".
#' @param reduction.key A character string specifying the prefix for the column names of the TriMap embeddings. Default is "TriMap_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the trimap.TRIMAP function.
#'
#' @examples
#' \dontrun{
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
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
#' @rdname RunTriMap
#' @export
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
  ...
) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 2) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- Matrix::as.matrix(
      Matrix::t(
        GetAssayData5(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " TriMap components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " TriMap components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims, features")
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
    seed.use = seed.use
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
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_python("trimap")
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
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
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
