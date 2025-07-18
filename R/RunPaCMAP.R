#' Run PaCMAP (Pairwise Controlled Manifold Approximation)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer to be used. Default is "data".
#' @param n_components An integer specifying the number of PaCMAP components. Default is 2.
#' @param n.neighbors An integer specifying the number of neighbors considered in the k-Nearest Neighbor graph. Default to 10 for dataset whose sample size is smaller than 10000. For large dataset whose sample size (n) is larger than 10000, the default value is: 10 + 15 * (log10(n) - 4).
#' @param MN_ratio A numeric value specifying the ratio of the ratio of the number of mid-near pairs to the number of neighbors. Default is 0.5.
#' @param FP_ratio A numeric value specifying the ratio of the ratio of the number of further pairs to the number of neighbors. Default is 2.
#' @param distance_method A character string specifying the distance metric to be used. Default is "euclidean".
#' @param lr A numeric value specifying the learning rate of the AdaGrad optimizer. Default is 1.
#' @param num_iters An integer specifying the number of iterations for PaCMAP optimization. Default is 450.
#' @param apply_pca A logical value indicating whether pacmap should apply PCA to the data before constructing the k-Nearest Neighbor graph. Using PCA to preprocess the data can largely accelerate the DR process without losing too much accuracy. Notice that this option does not affect the initialization of the optimization process. Default is TRUE.
#' @param init A character string specifying the initialization of the lower dimensional embedding. One of "pca" or "random". Default is "random".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "pacmap".
#' @param reduction.key A character string specifying the prefix for the column names of the PaCMAP embeddings. Default is "PaCMAP_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the pacmap.PaCMAP function.
#'
#' @rdname RunPaCMAP
#' @export
#'
#' @examples
#' \dontrun{
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunPaCMAP(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "pacmap"
#' )
#' }
RunPaCMAP <- function(object, ...) {
  UseMethod(generic = "RunPaCMAP", object = object)
}

#' @rdname RunPaCMAP
#' @method RunPaCMAP Seurat
#' @export
RunPaCMAP.Seurat <- function(
    object,
    reduction = "pca",
    dims = NULL,
    features = NULL,
    assay = NULL,
    layer = "data",
    n_components = 2,
    n.neighbors = NULL,
    MN_ratio = 0.5,
    FP_ratio = 2,
    distance_method = "euclidean",
    lr = 1,
    num_iters = 450L,
    apply_pca = TRUE,
    init = "random",
    reduction.name = "pacmap",
    reduction.key = "PaCMAP_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) < 1) {
    log_message(
      "Please specify only one of the following arguments: dims, features, or graph",
      message_type = "error"
    )
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
      log_message(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " PaCMAP components requested",
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
        " PaCMAP components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of dims or features",
      message_type = "error"
    )
  }
  object[[reduction.name]] <- RunPaCMAP(
    object = data.use,
    assay = assay,
    n_components = n_components,
    n.neighbors = n.neighbors,
    MN_ratio = MN_ratio,
    FP_ratio = FP_ratio,
    distance_method = distance_method,
    lr = lr,
    num_iters = num_iters,
    apply_pca = apply_pca,
    init = init,
    reduction.key = reduction.key,
    verbose = verbose,
    seed.use = seed.use
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunPaCMAP
#' @method RunPaCMAP default
#' @export
RunPaCMAP.default <- function(
    object,
    assay = NULL,
    n_components = 2,
    n.neighbors = NULL,
    MN_ratio = 0.5,
    FP_ratio = 2,
    distance_method = "euclidean",
    lr = 1,
    num_iters = 450L,
    apply_pca = TRUE,
    init = "random",
    reduction.key = "PaCMAP_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_python("pacmap")
  pacmap <- reticulate::import("pacmap")

  operator <- pacmap$PaCMAP(
    n_components = as.integer(n_components),
    n_neighbors = n.neighbors,
    MN_ratio = MN_ratio,
    FP_ratio = FP_ratio,
    distance = distance_method,
    lr = lr,
    num_iters = num_iters,
    apply_pca = apply_pca,
    verbose = verbose,
    random_state = as.integer(seed.use)
  )
  embedding <- operator$fit_transform(object, init = init)

  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(reduction)
}
