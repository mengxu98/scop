#' Run Force-Directed Layout (Fruchterman-Reingold algorithm)
#'
#' @md
#' @param object An object. This can be a Seurat object, a Neighbor object, or a Graph object.
#' @param reduction A character string specifying the reduction to be used.
#' Default is NULL.
#' @param dims An integer vector specifying the dimensions to be used.
#' Default is NULL.
#' @param features A character vector specifying the features to be used.
#' Default is NULL.
#' @param assay A character string specifying the assay to be used.
#' Default is NULL.
#' @param layer A character string specifying the layer to be used.
#' Default is "data".
#' @param graph A character string specifying the name of the Graph object to be used.
#' Default is NULL.
#' @param neighbor A character string specifying the name of the Neighbor object to be used.
#' Default is NULL.
#' @param k.param An integer specifying the number of nearest neighbors to consider.
#' Default is 20.
#' @param ndim An integer specifying the number of dimensions for the force-directed layout.
#' Default is 2.
#' @param niter An integer specifying the number of iterations for the force-directed layout.
#' Default is 500.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object.
#' Default is "fr".
#' @param reduction.key A character string specifying the prefix for the column names of the force-directed layout embeddings.
#' Default is "FR_".
#' @param verbose A logical value indicating whether to print verbose output.
#' Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used.
#' Default is 11.
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
#'   features = SeuratObject::VariableFeatures(pancreas_sub)
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
    ...) {
  log_message(
    "Running force-directed layout",
    message_type = "info"
  )
  if (sum(c(is.null(dims), is.null(features), is.null(neighbor), is.null(graph))) == 4) {
    log_message(
      "Please specify only one of the following arguments: dims, features, neighbor or graph",
      message_type = "error"
    )
  }
  if (!is.null(graph)) {
    if (!inherits(object[[graph]], what = "Graph")) {
      log_message(
        "Please specify a Graph object name, ",
        "instead of the name of a ",
        class(object[[graph]]),
        " object",
        message_type = "error"
      )
    }
    data_use <- object[[graph]]
  } else if (!is.null(neighbor)) {
    if (!inherits(object[[neighbor]], what = "Neighbor")) {
      log_message(
        "Please specify a Neighbor object name, ",
        "instead of the name of a ",
        class(object[[neighbor]]),
        " object",
        message_type = "error"
      )
    }
    data_use <- object[[neighbor]]
  } else if (!is.null(features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data_use <- GetAssayData5(
      object = object,
      layer = layer,
      assay = assay
    )
    data_use <- Matrix::t(data_use[features, ])
    log_message(
      "Computing nearest neighbor graph and SNN",
      message_type = "info"
    )
    data_use <- Seurat::FindNeighbors(
      data_use,
      k.param = k.param,
      verbose = FALSE
    )[["snn"]]
  } else if (!is.null(dims)) {
    data_use <- Seurat::Embeddings(
      object = object[[reduction]]
    )
    if (max(dims) > ncol(x = data_use)) {
      log_message(
        "More dimensions specified in dims than have been computed",
        message_type = "error"
      )
    }
    data_use <- data_use[, dims]
    data_use <- Seurat::FindNeighbors(
      data_use,
      k.param = k.param
    )[["snn"]]
  } else {
    log_message(
      "Please specify one of dims, features, neighbor, or graph",
      message_type = "error"
    )
  }

  object[[reduction.name]] <- RunFR(
    object = data_use,
    assay = assay,
    ndim = ndim,
    niter = niter,
    seed.use = seed.use,
    verbose = verbose,
    reduction.key = reduction.key,
    ...
  )
  object <- Seurat::LogSeuratCommand(object = object)

  log_message(
    "Force-directed layout computed",
    message_type = "info"
  )
  return(object)
}

#' @rdname RunFR
#' @method RunFR default
#' @export
RunFR.default <- function(
    object,
    assay = NULL,
    ndim = 2,
    niter = 500,
    reduction.key = "FR_",
    verbose = TRUE,
    seed.use = 11L,
    ...) {
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  if (inherits(object, "Neighbor")) {
    object <- Seurat::as.Graph(object)
  }
  g <- igraph::graph_from_adjacency_matrix(
    object,
    weighted = TRUE
  )
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
    assay = assay,
    global = TRUE
  )

  return(reduction)
}
