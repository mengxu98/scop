#' Run LargeVis (Dimensionality Reduction with a LargeVis-like method)
#'
#' @md
#' @inheritParams uwot::lvish
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used.
#' Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used.
#' Default is NULL.
#' @param features A character vector specifying the features to be used.
#' Default is NULL.
#' @param assay A character string specifying the assay to be used.
#' Default is NULL.
#' @param layer A character string specifying the layer to be used.
#' Default is "data".
#' @param n_components An integer specifying the number of LargeVis components.
#' Default is 2.
#' @param pca_method Method to carry out any PCA dimensionality reduction when the pca parameter is specified.
#' Allowed values are: "irlba", "rsvd", "bigstatsr", "svd", "auto"(the default.
#' Uses "irlba", unless more than 50 case "svd" is used.)
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object.
#' Default is "largevis".
#' @param reduction.key A character string specifying the prefix for the column names of the LargeVis embeddings.
#' Default is "LargeVis_".
#' @param verbose A logical value indicating whether to print verbose output.
#' Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the [uwot::lvish] function.
#'
#' @rdname RunLargeVis
#' @export
#'
#' @examples
#' \dontrun{
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunLargeVis(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "largevis"
#' )
#' }
RunLargeVis <- function(object, ...) {
  UseMethod(generic = "RunLargeVis", object = object)
}

#' @rdname RunLargeVis
#' @method RunLargeVis Seurat
#' @export
RunLargeVis.Seurat <- function(
    object,
    reduction = "pca",
    dims = NULL,
    features = NULL,
    assay = NULL,
    layer = "data",
    perplexity = 50,
    n_neighbors = perplexity * 3,
    n_components = 2,
    metric = "euclidean",
    n_epochs = -1,
    learning_rate = 1,
    scale = "maxabs",
    init = "lvrandom",
    init_sdev = NULL,
    repulsion_strength = 7,
    negative_sample_rate = 5,
    nn_method = NULL,
    n_trees = 50,
    search_k = 2 * n_neighbors * n_trees,
    n_threads = NULL,
    n_sgd_threads = 0,
    grain_size = 1,
    kernel = "gauss",
    pca = NULL,
    pca_center = TRUE,
    pcg_rand = TRUE,
    fast_sgd = FALSE,
    batch = FALSE,
    opt_args = NULL,
    epoch_callback = NULL,
    pca_method = NULL,
    reduction.name = "largevis",
    reduction.key = "LargeVis_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 3) {
    log_message(
      "Please specify only one of the following arguments: dims, features",
      message_type = "error"
    )
  }
  if (!is.null(x = features)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object = object)
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
      log_message(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " LargeVis components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Seurat::Embeddings(
      object[[reduction]]
    )[, dims]
    assay <- SeuratObject::DefaultAssay(
      object = object[[reduction]]
    )
    if (length(x = dims) < n_components) {
      log_message(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " LargeVis components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of dims, features",
      message_type = "error"
    )
  }
  object[[reduction.name]] <- RunLargeVis(
    object = data.use,
    assay = assay,
    perplexity = perplexity,
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    scale = scale,
    init = init,
    init_sdev = init_sdev,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    nn_method = nn_method,
    n_trees = n_trees,
    search_k = search_k,
    n_threads = n_threads,
    n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    kernel = kernel,
    pca = pca,
    pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    batch = batch,
    opt_args = opt_args,
    epoch_callback = epoch_callback,
    pca_method = pca_method,
    reduction.key = reduction.key,
    verbose = verbose,
    seed.use = seed.use
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunLargeVis
#' @method RunLargeVis default
#' @export
RunLargeVis.default <- function(
    object,
    assay = NULL,
    perplexity = 50,
    n_neighbors = perplexity * 3,
    n_components = 2,
    metric = "euclidean",
    n_epochs = -1,
    learning_rate = 1,
    scale = "maxabs",
    init = "lvrandom",
    init_sdev = NULL,
    repulsion_strength = 7,
    negative_sample_rate = 5,
    nn_method = NULL,
    n_trees = 50,
    search_k = 2 * n_neighbors * n_trees,
    n_threads = NULL,
    n_sgd_threads = 0,
    grain_size = 1,
    kernel = "gauss",
    pca = NULL,
    pca_center = TRUE,
    pcg_rand = TRUE,
    fast_sgd = FALSE,
    batch = FALSE,
    opt_args = NULL,
    epoch_callback = NULL,
    pca_method = NULL,
    reduction.key = "LargeVis_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  embedding <- uwot::lvish(
    X = object,
    perplexity = perplexity,
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    scale = scale,
    init = init,
    init_sdev = init_sdev,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    nn_method = nn_method,
    n_trees = n_trees,
    search_k = search_k,
    n_threads = n_threads,
    n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    kernel = kernel,
    pca = pca,
    pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    verbose = verbose,
    batch = batch,
    opt_args = opt_args,
    epoch_callback = epoch_callback,
    pca_method = pca_method,
    ...
  )
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
