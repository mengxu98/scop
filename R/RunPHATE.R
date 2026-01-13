#' @title Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunUMAP2
#' @inheritParams RunDM
#' @param n_components The number of PHATE components.
#' Default is `2`.
#' @param knn A number of nearest neighbors on which to build kernel.
#' Default is `5`.
#' @param decay The sets decay rate of kernel tails.
#' Default is `40`.
#' @param n_landmark A number of landmarks to use in fast PHATE.
#' Default is `2000`.
#' @param t The power to which the diffusion operator is powered.
#' This sets the level of diffusion.
#' If `"auto"`, `t` is selected according to the knee point in the Von Neumann Entropy of the diffusion operator.
#' Default is `"auto"`.
#' @param gamma The informational distance constant between `-1` and `1`.
#' `gamma=1` gives the PHATE log potential, `gamma=0` gives a square root potential.
#' Default is `1`.
#' @param n_pca A number of principal components to use for calculating neighborhoods.
#' For extremely large datasets, using `n_pca < 20` allows neighborhoods to be calculated in roughly `log(n_samples)` time.
#' Default is `100`.
#' @param knn_dist The distance metric for k-nearest neighbors.
#' Recommended values: `"euclidean"`, `"cosine"`, `"precomputed"`.
#' Default is `"euclidean"`.
#' @param knn_max The maximum number of neighbors for which alpha decaying kernel is computed for each point.
#' For very large datasets, setting `knn_max` to a small multiple of `knn` can speed up computation significantly.
#' Default is `NULL`.
#' @param t_max The maximum `t` to test.
#' Default is `100`.
#' @param do_cluster Whether to perform clustering on the PHATE embeddings.
#' Default is `FALSE`.
#' @param n_clusters A number of clusters to be identified.
#' Default is `"auto"`.
#' @param max_clusters The maximum number of clusters to test.
#' Default is `100`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"phate"`.
#' @param reduction.key The prefix for the column names of the PHATE embeddings.
#' Default is `"PHATE_"`.
#' @param ... Additional arguments to be passed to phate.PHATE.
#'
#' @rdname RunPHATE
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunPHATE(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "phate"
#' )
#' }
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
    seed.use = 11,
    ...) {
  if (sum(c(is.null(dims), is.null(features))) == 2) {
    log_message(
      "Please specify only one of the following arguments: dims, features",
      message_type = "error"
    )
  }
  if (!is.null(features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data_use <- as_matrix(
      Matrix::t(
        GetAssayData5(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(data_use) < n_components) {
      log_message(
        "Please provide as many or more features than n_components: ",
        length(features),
        " features provided, ",
        n_components,
        " PHATE components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(dims)) {
    data_use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(dims) < n_components) {
      log_message(
        "Please provide as many or more dims than n_components: ",
        length(dims),
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
    object = data_use,
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
    seed.use = 11,
    ...) {
  set.seed(seed = seed.use)
  PrepareEnv()
  check_python("phate", verbose = verbose)
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
  colnames(embedding) <- paste0(reduction.key, seq_len(ncol(embedding)))
  rownames(embedding) <- rownames(object)

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
