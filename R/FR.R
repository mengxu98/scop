#' Run Force-Directed Layout (Fruchterman-Reingold algorithm)
#'
#' @md
#' @param object An object. This can be a Seurat object, a Neighbor object, or a Graph object.
#' @param reduction A character string specifying the reduction to be used. Default is NULL.
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param layer A character string specifying the layer to be used. Default is "data".
#' @param graph A character string specifying the name of the Graph object to be used. Default is NULL.
#' @param neighbor A character string specifying the name of the Neighbor object to be used. Default is NULL.
#' @param k.param An integer specifying the number of nearest neighbors to consider. Default is 20.
#' @param ndim An integer specifying the number of dimensions for the force-directed layout. Default is 2.
#' @param niter An integer specifying the number of iterations for the force-directed layout. Default is 500.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "fr".
#' @param reduction.key A character string specifying the prefix for the column names of the force-directed layout embeddings. Default is "FR_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the [igraph::layout_with_fr] function.
#'
#' @export
#'
#' @rdname RunFR
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunFR(
#'   object = pancreas_sub,
#'   features = Seurat::VariableFeatures(pancreas_sub)
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   reduction = "fr"
#' )
RunFR <- function(object, ...) {
  UseMethod(generic = "RunFR", object = object)
}

#' @rdname RunFR
#' @method RunFR Seurat
#' @export
RunFR.Seurat <- function(
  object,
  reduction = NULL,
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  graph = NULL,
  neighbor = NULL,
  k.param = 20,
  ndim = 2,
  niter = 500,
  reduction.name = "FR",
  reduction.key = "FR_",
  verbose = TRUE,
  seed.use = 11L,
  ...
) {
  if (
    sum(c(
      is.null(x = dims),
      is.null(x = features),
      is.null(neighbor),
      is.null(x = graph)
    )) ==
      4
  ) {
    stop(
      "Please specify only one of the following arguments: dims, features, neighbor or graph"
    )
  }
  if (!is.null(x = graph)) {
    if (!inherits(x = object[[graph]], what = "Graph")) {
      stop(
        "Please specify a Graph object name, ",
        "instead of the name of a ",
        class(object[[graph]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[graph]]
  } else if (!is.null(x = neighbor)) {
    if (!inherits(x = object[[neighbor]], what = "Neighbor")) {
      stop(
        "Please specify a Neighbor object name, ",
        "instead of the name of a ",
        class(object[[neighbor]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[neighbor]]
  } else if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- Matrix::t(
      Seurat::GetAssayData(
        object = object,
        layer = layer,
        assay = assay
      )[features, ]
    )
    data.use <- Seurat::FindNeighbors(
      data.use,
      k.param = k.param
    )[["snn"]]
  } else if (!is.null(x = dims)) {
    data.use <- Seurat::Embeddings(
      object = object[[reduction]]
    )
    if (max(dims) > ncol(x = data.use)) {
      stop("More dimensions specified in dims than have been computed")
    }
    data.use <- data.use[, dims]
    data.use <- Seurat::FindNeighbors(
      data.use,
      k.param = k.param
    )[["snn"]]
  } else {
    stop("Please specify one of dims, features, neighbor, or graph")
  }
  object[[reduction.name]] <- RunFR(
    object = data.use,
    ndim = ndim,
    niter = niter,
    seed.use = seed.use,
    verbose = verbose,
    reduction.key = reduction.key,
    ...
  )
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunFR
#' @method RunFR default
#' @export
RunFR.default <- function(
  object,
  ndim = 2,
  niter = 500,
  reduction.key = "FR_",
  verbose = TRUE,
  seed.use = 11L,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (inherits(object, "Neighbor")) {
    object <- Seurat::as.Graph(object)
  }
  g <- igraph::graph_from_adjacency_matrix(object, weighted = TRUE)
  embedding <- igraph::layout_with_fr(
    graph = g,
    dim = ndim,
    niter = niter,
    ...
  )
  colnames(x = embedding) <- paste0(
    reduction.key,
    seq_len(ncol(x = embedding))
  )
  rownames(x = embedding) <- rownames(object)
  reduction <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = object@assay.used,
    global = TRUE
  )

  return(reduction)
}
