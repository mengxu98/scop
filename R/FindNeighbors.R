#' @export
FindNeighbors.Seurat <- function(
  object,
  reduction = "pca",
  dims = 1:10,
  assay = NULL,
  features = NULL,
  k.param = 20,
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1 / 15,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  do.plot = FALSE,
  graph.name = NULL,
  l2.norm = FALSE,
  cache.index = FALSE,
  ...
) {
  seurat_find_neighbors <- get("FindNeighbors.Seurat", envir = asNamespace("Seurat"))
  seurat_find_neighbors(
    object = object,
    reduction = reduction,
    dims = dims,
    assay = assay,
    features = features,
    k.param = k.param,
    return.neighbor = return.neighbor,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.method = nn.method,
    n.trees = n.trees,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    verbose = verbose,
    do.plot = do.plot,
    graph.name = graph.name,
    l2.norm = l2.norm,
    cache.index = cache.index,
    ...
  )
}

#' Find nearest neighbors
#'
#' @param object Object containing reduced-dimensional data.
#' @param ... Passed to methods.
#'
#' @return The input object with neighbor information added.
#' @export
FindNeighbors <- function(object, ...) {
  UseMethod("FindNeighbors")
}

#' @export
FindNeighbors.default <- function(object, ...) {
  stop("FindNeighbors supports Seurat objects.", call. = FALSE)
}
