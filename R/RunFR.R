#' @title Run Force-Directed Layout (Fruchterman-Reingold algorithm)
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param object An object. This can be a Seurat object, a Neighbor object, or a Graph object.
#' @param reduction The reduction to be used.
#' Default is `NULL`.
#' @param dims The dimensions to be used.
#' Default is `NULL`.
#' @param features The features to be used.
#' Default is `NULL`.
#' @param assay The assay to be used.
#' Default is `NULL`.
#' @param layer The layer to be used.
#' Default is `"data"`.
#' @param graph The name of the Graph object to be used.
#' Default is `NULL`.
#' @param neighbor The name of the Neighbor object to be used.
#' Default is `NULL`.
#' @param k.param The number of nearest neighbors to consider.
#' Default is `20`.
#' @param ndim The number of dimensions for the force-directed layout.
#' Default is `2`.
#' @param niter The number of iterations for the force-directed layout.
#' Default is `500`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"fr"`.
#' @param reduction.key The prefix for the column names of the force-directed layout embeddings.
#' Default is `"FR_"`.
#' @param seed.use The random seed to be used.
#' Default is `11`.
#' @param ... Additional arguments to be passed to [igraph::layout_with_fr].
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
    verbose = verbose
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
      assay = assay,
      verbose = FALSE
    )
    data_use <- Matrix::t(data_use[features, ])
    log_message(
      "Computing nearest neighbor graph and SNN",
      verbose = verbose
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
    verbose = verbose
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
