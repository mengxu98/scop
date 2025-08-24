#' Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer to be used. Default is "data".
#' @param n_components An integer specifying the number of PHATE components. Default is 2.
#' @param knn An integer specifying the number of nearest neighbors on which to build kernel. Default is 5.
#' @param decay An integer specifying the sets decay rate of kernel tails. Default is 40.
#' @param n_landmark An integer specifying the number of landmarks to use in fast PHATE. Default is 2000.
#' @param t A character string specifying the power to which the diffusion operator is powered. This sets the level of diffusion. If ‘auto’, t is selected according to the knee point in the Von Neumann Entropy of the diffusion operator. Default is "auto".
#' @param gamma A numeric value specifying the informational distance constant between -1 and 1. gamma=1 gives the PHATE log potential, gamma=0 gives a square root potential. Default is 1.
#' @param n_pca An integer specifying the number of principal components to use for calculating neighborhoods. For extremely large datasets, using n_pca < 20 allows neighborhoods to be calculated in roughly log(n_samples) time. Default is 100.
#' @param knn_dist A character string specifying the distance metric for k-nearest neighbors. Recommended values: "euclidean, "cosine, "precomputed". Default is "euclidean".
#' @param knn_max An integer specifying the maximum number of neighbors for which alpha decaying kernel is computed for each point. For very large datasets, setting knn_max to a small multiple of knn can speed up computation significantly. Default is NULL.
#' @param t_max An integer specifying the maximum \code{t} to test. Default is 100.
#' @param do_cluster A logical value indicating whether to perform clustering on the PHATE embeddings. Default is FALSE.
#' @param n_clusters An integer specifying the number of clusters to be identified. Default is "auto".
#' @param max_clusters An integer specifying the maximum number of clusters to test. Default is 100.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "phate".
#' @param reduction.key A character string specifying the prefix for the column names of the PHATE embeddings. Default is "PHATE_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the phate.PHATE function.
#'
#' @rdname RunPHATE
#' @export
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunPHATE(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "phate"
#' )
RunPHATE <- function(object, ...) {
  UseMethod(generic = "RunPHATE", object = object)
}

#' @rdname RunPHATE
#' @method RunPHATE Seurat
#' @export
RunPHATE.Seurat <- function(
    object,
    reduction = "pca",
    dims = NULL,
    features = NULL,
    assay = NULL,
    layer = "data",
    n_components = 2,
    knn = 5,
    decay = 40,
    n_landmark = 2000,
    t = "auto",
    gamma = 1,
    n_pca = 100,
    knn_dist = "euclidean",
    knn_max = NULL,
    t_max = 100,
    do_cluster = FALSE,
    n_clusters = "auto",
    max_clusters = 100,
    reduction.name = "phate",
    reduction.key = "PHATE_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 2) {
    log_message(
      "Please specify only one of the following arguments: dims, features",
      message_type = "error"
    )
  }
  if (!is.null(x = features)) {
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
    if (ncol(x = data.use) < n_components) {
      log_message(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " PHATE components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      log_message(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " PHATE components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of dims or features",
      message_type = "error"
    )
  }
  object[[reduction.name]] <- RunPHATE(
    object = data.use,
    assay = assay,
    n_components = n_components,
    knn = knn,
    decay = decay,
    n_landmark = n_landmark,
    t = t,
    gamma = gamma,
    n_pca = n_pca,
    knn_dist = knn_dist,
    knn_max = knn_max,
    t_max = t_max,
    do_cluster = do_cluster,
    n_clusters = n_clusters,
    max_clusters = max_clusters,
    reduction.key = reduction.key,
    verbose = verbose,
    seed.use = seed.use
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunPHATE
#' @method RunPHATE default
#' @export
RunPHATE.default <- function(
    object,
    assay = NULL,
    n_components = 2,
    knn = 5,
    decay = 40,
    n_landmark = 2000,
    t = "auto",
    gamma = 1,
    n_pca = 100,
    knn_dist = "euclidean",
    knn_max = NULL,
    t_max = 100,
    do_cluster = FALSE,
    n_clusters = "auto",
    max_clusters = 100,
    reduction.key = "PHATE_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  check_python("phate")
  phate <- reticulate::import("phate")

  if (is.numeric(knn_max) && length(knn_max) > 0) {
    knn_max <- as.integer(knn_max)
  } else {
    knn_max <- NULL
  }
  operator <- phate$PHATE(
    n_components = as.integer(n_components),
    knn = as.integer(knn),
    decay = as.integer(decay),
    n_landmark = as.integer(n_landmark),
    t = as.character(t),
    gamma = as.numeric(gamma),
    n_pca = as.integer(n_pca),
    knn_dist = as.character(knn_dist),
    knn_max = knn_max,
    random_state = as.integer(seed.use),
    verbose = as.integer(verbose),
    ...
  )
  embedding <- operator$fit_transform(object, t_max = as.integer(t_max))
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  if (isTRUE(do_cluster)) {
    if (is.numeric(n_clusters)) {
      n_clusters <- as.integer(n_clusters)
    }
    if (is.numeric(max_clusters)) {
      max_clusters <- as.integer(max_clusters)
    }
    clusters <- phate$cluster$kmeans(
      operator,
      n_clusters = n_clusters,
      max_clusters = max_clusters,
      random_state = as.integer(seed.use)
    )
    clusters <- clusters + 1
    clusters <- factor(clusters, levels = sort(unique(clusters)))
    names(clusters) <- rownames(embedding)
    SeuratObject::Misc(reduction, slot = "clusters") <- clusters
  }
  return(reduction)
}
