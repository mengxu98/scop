#' Run DM (diffusion map)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used.
#' @param dims An integer vector specifying the dimensions to be used. Default is 1:30.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer to be used. Default is "data".
#' @param ndcs An integer specifying the number of diffusion components (dimensions) to be computed. Default is 2.
#' @param sigma A character string specifying the diffusion scale parameter of the Gaussian kernel. Currently supported values are "local" (default) and "global".
#' @param k An integer specifying the number of nearest neighbors to be used for the construction of the graph. Default is 30.
#' @param dist.method A character string specifying the distance metric to be used for the construction of the knn graph. Currently supported values are "euclidean" and "cosine". Default is "euclidean".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "dm".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "DM_".
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param ... Additional arguments to be passed to the \link[destiny]{DiffusionMap} function.
#'
#' @rdname RunDM
#' @export
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunDM(
#'   object = pancreas_sub,
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
#' )
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
  ...
) {
  if (!is.null(x = features)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object = object)
    data.use <- Matrix::as.matrix(
      Matrix::t(
        SeuratObject::GetAssayData(
          object = object,
          layer = layer,
          assay = assay
        )[features, ]
      )
    )
    if (ncol(x = data.use) < ndcs) {
      stop(
        "Please provide as many or more features than ndcs: ",
        length(x = features),
        " features provided, ",
        ndcs,
        " Diffusion components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- SeuratObject::Embeddings(object[[reduction]])[, dims]
    assay <- SeuratObject::DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < ndcs) {
      stop(
        "Please provide as many or more dims than ndcs: ",
        length(x = dims),
        " dims provided, ",
        ndcs,
        " DiffusionMap components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims, features")
  }
  reduction.data <- RunDM(
    object = data.use,
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
  object[[reduction.name]] <- reduction.data
  object <- Seurat::LogSeuratCommand(object = object)
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
  ...
) {
  check_r("destiny")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  dm.results <- destiny::DiffusionMap(
    data = Matrix::as.matrix(object),
    n_eigs = ndcs,
    sigma = sigma,
    k = k,
    distance = dist.method,
    verbose = verbose,
    ...
  )

  cell.embeddings <- dm.results@eigenvectors
  rownames(x = cell.embeddings) <- rownames(object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ndcs)
  reduction <- SeuratObject::CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = layer, dm.results = dm.results)
  )
  return(reduction)
}
