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
  extra <- list(...)
  seurat_find_neighbors <- get("FindNeighbors.Seurat", envir = asNamespace("Seurat"))
  run_seurat <- function() {
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
  supported <- length(extra) == 0L &&
    !is.null(reduction) &&
    length(reduction) == 1L &&
    reduction %in% SeuratObject::Reductions(object) &&
    !is.null(dims) &&
    is.null(features) &&
    isFALSE(return.neighbor) &&
    isTRUE(compute.SNN) &&
    identical(nn.method, "annoy") &&
    identical(annoy.metric, "euclidean") &&
    isTRUE(all.equal(as.numeric(nn.eps), 0)) &&
    isFALSE(do.plot) &&
    isFALSE(l2.norm) &&
    isFALSE(cache.index)
  if (!isTRUE(supported)) {
    return(run_seurat())
  }

  emb <- SeuratObject::Embeddings(object[[reduction]])
  dims <- as.integer(dims)
  dims <- dims[is.finite(dims) & dims >= 1L & dims <= ncol(emb)]
  if (!length(dims)) {
    return(run_seurat())
  }
  data.use <- emb[, dims, drop = FALSE]
  if (!is.double(data.use)) {
    storage.mode(data.use) <- "double"
  }
  n_cells <- nrow(data.use)
  if (n_cells < 2L) {
    return(run_seurat())
  }
  k.use <- min(as.integer(k.param), n_cells - 1L)
  if (!is.finite(k.use) || k.use < 1L) {
    return(run_seurat())
  }
  n_trees <- max(1L, as.integer(n.trees))
  cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  cores <- max(1L, min(8L, as.integer(cores)))
  raw_idx <- tryCatch(
    annoy_build_search(
      data = data.use,
      k = min(k.use + 1L, n_cells),
      n_trees = n_trees,
      cores = cores
    ),
    error = function(e) NULL
  )
  if (is.null(raw_idx)) {
    return(run_seurat())
  }
  nn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = k.use)
  for (i in seq_len(n_cells)) {
    hits <- raw_idx[i, ]
    hits <- hits[!is.na(hits) & hits != i]
    if (length(hits) < k.use) {
      fill <- setdiff(seq_len(n_cells), c(i, hits))
      hits <- c(hits, fill[seq_len(min(length(fill), k.use - length(hits)))])
    }
    nn_idx[i, ] <- hits[seq_len(k.use)]
  }

  cell_names <- rownames(data.use)
  if (is.null(cell_names)) {
    cell_names <- colnames(object)
  }
  assay.use <- assay %||% object[[reduction]]@assay.used
  if (length(assay.use) != 1L || is.na(assay.use) || !nzchar(assay.use)) {
    assay.use <- SeuratObject::DefaultAssay(object)
  }
  graph_names <- graph.name
  if (is.null(graph_names)) {
    graph_names <- paste0(assay.use, c("_nn", "_snn"))
  }
  if (length(graph_names) == 1L) {
    graph_names <- c(graph_names, paste0(graph_names, "_snn"))
  }

  nn_i <- rep(seq_len(n_cells), each = k.use)
  nn_j <- as.vector(t(nn_idx))
  nn_graph <- Matrix::sparseMatrix(
    i = nn_i,
    j = nn_j,
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(cell_names, cell_names)
  )
  nn_graph <- methods::as(nn_graph, "Graph")
  nn_graph@assay.used <- assay.use

  membership <- Matrix::sparseMatrix(
    i = nn_i,
    j = nn_j,
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(cell_names, cell_names)
  )
  snn <- Matrix::tcrossprod(membership)
  snn@x <- snn@x / (2 * k.use - snn@x)
  snn@x[snn@x <= prune.SNN] <- 0
  snn <- Matrix::drop0(snn)
  Matrix::diag(snn) <- 0
  snn <- Matrix::drop0(snn)
  snn <- methods::as(snn, "generalMatrix")
  snn <- methods::as(snn, "Graph")
  snn@assay.used <- assay.use

  object@graphs[[graph_names[[1L]]]] <- nn_graph
  object@graphs[[graph_names[[2L]]]] <- snn
  return(SeuratObject::LogSeuratCommand(object = object))
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
