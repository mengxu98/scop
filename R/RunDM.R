#' @title Run DM (diffusion map)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction The reduction to be used.
#' Default is `"pca"`.
#' @param dims The dimensions to be used.
#' Default is `1:30`.
#' @param features The features to be used.
#' Default is `NULL`.
#' @param assay The assay to be used.
#' Default is `NULL`.
#' @param layer The layer to be used.
#' Default is `"data"`.
#' @param ndcs A number of diffusion components (dimensions) to be computed.
#' Default is `2`.
#' @param sigma The diffusion scale parameter of the Gaussian kernel.
#' Currently supported values are `"local"` (default) and `"global"`.
#' @param k A number of nearest neighbors to be used for the construction of the graph.
#' Default is `30`.
#' @param dist.method The distance metric to be used for the construction of the knn graph.
#' Currently supported values are `"euclidean"` and `"cosine"`.
#' Default is `"euclidean"`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"dm"`.
#' @param reduction.key The prefix for the column names of the basis vectors.
#' Default is `"DM_"`.
#' @param seed.use An integer specifying the random seed to be used.
#' Default is `11`.
#' @param ... Additional arguments to be passed to [destiny::DiffusionMap].
#'
#' @rdname RunDM
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDM(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "dm"
#' )
RunDM <- function(object, ...) {
  UseMethod(generic = "RunDM", object = object)
}

#' @rdname RunDM
#' @method RunDM Seurat
#' @export
RunDM.Seurat <- function(
    object,
    reduction = "pca",
    dims = 1:30,
    features = NULL,
    assay = NULL,
    layer = "data",
    ndcs = 2,
    sigma = "local",
    k = 30,
    dist.method = "euclidean",
    reduction.name = "dm",
    reduction.key = "DM_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  log_message(
    "Running {.pkg destiny::DiffusionMap}",
    message_type = "running",
    verbose = verbose
  )
  if (!is.null(features)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object = object)
    data_use <- as_matrix(
      Matrix::t(
        GetAssayData5(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(data_use) < ndcs) {
      log_message(
        "Please provide as many or more features than ndcs: {.val {length(features)}}",
        " features provided, {.val {ndcs}} Diffusion components requested",
        message_type = "error"
      )
    }
  } else if (!is.null(dims)) {
    reduction <- DefaultReduction(
      object,
      pattern = reduction
    )
    data_use <- SeuratObject::Embeddings(object[[reduction]])[, dims]
    assay <- SeuratObject::DefaultAssay(object = object[[reduction]])
    if (length(dims) < ndcs) {
      log_message(
        "Please provide as many or more dims than ndcs: {.val {length(dims)}}",
        " dims provided, {.val {ndcs}} DiffusionMap components requested",
        message_type = "error"
      )
    }
  } else {
    log_message(
      "Please specify one of {.arg dims} and {.arg features}",
      message_type = "error"
    )
  }

  object[[reduction.name]] <- RunDM(
    object = data_use,
    assay = assay,
    layer = layer,
    ndcs = ndcs,
    sigma = sigma,
    k = k,
    dist.method = dist.method,
    reduction.key = reduction.key,
    seed.use = seed.use,
    verbose = verbose,
    ...
  )
  object <- Seurat::LogSeuratCommand(object = object)
  log_message(
    "{.pkg destiny::DiffusionMap} computed successfully and set {.val {reduction.name}} reduction",
    message_type = "success",
    verbose = verbose
  )

  return(object)
}

#' @rdname RunDM
#' @method RunDM default
#' @export
RunDM.default <- function(
    object,
    assay = NULL,
    layer = "data",
    ndcs = 2,
    sigma = "local",
    k = 30,
    dist.method = "euclidean",
    reduction.key = "DM_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  check_r("destiny", verbose = FALSE)

  set.seed(seed = seed.use)
  dm_results <- destiny::DiffusionMap(
    data = object,
    n_eigs = ndcs,
    sigma = sigma,
    k = k,
    distance = dist.method,
    verbose = FALSE,
    ...
  )

  embeddings <- dm_results@eigenvectors
  rownames(embeddings) <- rownames(object)
  colnames(embeddings) <- paste0(reduction.key, 1:ndcs)
  reduction <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(
      slot = layer,
      dm.results = dm_results
    )
  )
  return(reduction)
}
