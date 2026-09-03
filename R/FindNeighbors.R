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
  cores = 1L,
  ...
) {
  extra <- list(...)
  run_seurat <- function() {
    find_neighbors_seurat_seurat(
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
  if (length(extra) || isTRUE(do.plot)) {
    return(run_seurat())
  }

  if (!is.null(dims)) {
    if (
      length(reduction) != 1L ||
        !reduction %in% SeuratObject::Reductions(object)
    ) {
      return(run_seurat())
    }
    assay.use <- SeuratObject::DefaultAssay(object[[reduction]])
    data.use <- SeuratObject::Embeddings(object[[reduction]])
    if (!length(dims) || any(!is.finite(dims)) || any(dims < 1L)) {
      return(run_seurat())
    }
    if (max(dims) > ncol(data.use)) {
      stop("More dimensions specified in dims than have been computed")
    }
    data.use <- data.use[, dims, drop = FALSE]
    neighbor.graphs <- FindNeighbors.default(
      object = data.use,
      k.param = k.param,
      return.neighbor = return.neighbor,
      compute.SNN = compute.SNN,
      prune.SNN = prune.SNN,
      nn.method = nn.method,
      n.trees = n.trees,
      annoy.metric = annoy.metric,
      nn.eps = nn.eps,
      verbose = verbose,
      l2.norm = l2.norm,
      cache.index = cache.index,
      cores = cores
    )
  } else {
    assay.use <- assay %||% SeuratObject::DefaultAssay(object)
    neighbor.graphs <- FindNeighbors.Assay(
      object = object[[assay.use]],
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
      l2.norm = l2.norm,
      cache.index = cache.index,
      cores = cores
    )
  }

  if (length(neighbor.graphs) == 1L && !identical(names(neighbor.graphs), "nn")) {
    neighbor.graphs <- list(nn = neighbor.graphs)
  }
  graph.name <- graph.name %||% if (return.neighbor) {
    paste0(assay.use, ".", names(neighbor.graphs))
  } else {
    paste0(assay.use, "_", names(neighbor.graphs))
  }
  if (length(graph.name) == 1L && isTRUE(verbose)) {
    message("Only one graph name supplied, storing nearest-neighbor graph only")
  }
  for (ii in seq_along(graph.name)) {
    value <- neighbor.graphs[[ii]]
    if (inherits(value, "Graph")) {
      SeuratObject::DefaultAssay(value) <- assay.use
    }
    object[[graph.name[[ii]]]] <- value
  }
  SeuratObject::LogSeuratCommand(object = object)
}

# Convert a nearest-neighbour index matrix to reference cell names without an
# element-wise `apply()` call.  Matrix indexing is column-major just like the
# index matrices returned by Seurat, so this preserves both values and
# dimnames while avoiding one R closure invocation per neighbour.
knn_indices_to_names <- function(indices, reference_names) {
  out <- matrix(
    reference_names[as.integer(indices)],
    nrow = nrow(indices),
    ncol = ncol(indices),
    dimnames = dimnames(indices)
  )
  out
}

# `knn_cross_topk_native()` has already rejected non-finite values before this
# helper is reached.  `matrixStats::rowRanks()` therefore has the same ranking
# semantics as `rank()` here (including average tie ranks), while avoiding an
# R-level `apply()` call for every reference/query vector.
knn_rank_rows <- function(x) {
  check_r("matrixStats", verbose = FALSE)
  matrixStats::rowRanks(as.matrix(x), ties.method = "average")
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

find_neighbors_seurat_seurat <- function(object, ...) {
  method <- get("FindNeighbors.Seurat", envir = asNamespace("Seurat"))
  method(object = object, ...)
}

find_neighbors_default_seurat <- function(object, ...) {
  method <- get("FindNeighbors.default", envir = asNamespace("Seurat"))
  method(object = object, ...)
}

find_neighbors_native_distance <- function(object, k) {
  if (
    !is.matrix(object) || !is.numeric(object) ||
      nrow(object) != ncol(object) || any(!is.finite(object))
  ) {
    return(NULL)
  }
  result <- tryCatch(
    thisutils::run_dense_topk(
      object,
      k = k,
      by = "row",
      decreasing = FALSE
    ),
    error = function(e) NULL
  )
  if (
    is.null(result) || !is.matrix(result[["idx"]]) ||
      !identical(dim(result[["idx"]]), c(nrow(object), k))
  ) {
    return(NULL)
  }
  list(idx = result[["idx"]], distance = result[["value"]])
}

find_neighbors_l2_rows <- function(mat) {
  norms <- sqrt(rowSums(mat * mat))
  if (any(!is.finite(norms)) || any(norms <= 0)) {
    return(NULL)
  }
  mat / norms
}

# Map a finite embedding matrix to k nearest neighbors. Native code covers the
# exact RANN path only; Annoy always returns NULL so callers use Seurat.
# Returns NULL when the request is outside the validated native contract.
find_neighbors_native_embedding <- function(
  object,
  query = NULL,
  k,
  nn.method,
  n.trees,
  annoy.metric,
  nn.eps,
  cores,
  need_distance
) {
  if (
    !is.matrix(object) || !is.numeric(object) || any(!is.finite(object)) ||
      is.null(rownames(object)) ||
      !is.numeric(cores) || length(cores) != 1L || is.na(cores) ||
      cores < 1L || !isTRUE(all.equal(as.numeric(cores), as.integer(cores)))
  ) {
    return(NULL)
  }
  cores <- as.integer(cores)
  query_mat <- if (is.null(query)) {
    object
  } else {
    if (
      !is.matrix(query) || !is.numeric(query) || any(!is.finite(query)) ||
        ncol(query) != ncol(object) || is.null(rownames(query))
    ) {
      return(NULL)
    }
    query
  }

  if (identical(nn.method, "rann")) {
    if (
      !identical(annoy.metric, "euclidean") ||
        !is.numeric(nn.eps) || length(nn.eps) != 1L || is.na(nn.eps) ||
        nn.eps != 0
    ) {
      return(NULL)
    }
    if (isTRUE(need_distance) || !is.null(query)) {
      result <- tryCatch(
        cross_knn_f32(
          reference = object,
          query = query_mat,
          k = k,
          metric = "euclidean",
          cores = cores
        ),
        error = function(e) NULL
      )
      if (
        is.null(result) || !is.matrix(result[["idx"]]) ||
          !identical(dim(result[["idx"]]), c(nrow(query_mat), k))
      ) {
        return(NULL)
      }
      return(list(
        idx = result[["idx"]],
        distance = result[["distance"]],
        cell.names = rownames(query_mat),
        metric = "euclidean",
        ndim = ncol(object)
      ))
    }
    idx <- tryCatch(
      exact_knn_f32(data = object, k = k, cores = cores),
      error = function(e) NULL
    )
    if (is.null(idx) || !identical(dim(idx), c(nrow(object), k))) {
      return(NULL)
    }
    return(list(
      idx = idx,
      distance = NULL,
      cell.names = rownames(object),
      metric = "euclidean",
      ndim = ncol(object)
    ))
  }

  NULL
}

find_neighbors_as_neighbor <- function(idx, distance, cell_names, metric, ndim) {
  idx_mat <- as.matrix(idx)
  storage.mode(idx_mat) <- "numeric"
  dist_mat <- as.matrix(distance)
  storage.mode(dist_mat) <- "numeric"
  methods::new(
    Class = "Neighbor",
    nn.idx = idx_mat,
    nn.dist = dist_mat,
    alg.info = list(metric = metric, ndim = as.integer(ndim)),
    cell.names = cell_names
  )
}

find_neighbors_as_graphs <- function(indices, cell_names, compute.SNN, prune.SNN) {
  n_cells <- length(cell_names)
  k <- ncol(indices)
  if (!n_cells || nrow(indices) != n_cells || anyNA(indices)) {
    return(NULL)
  }
  i <- rep(seq_len(n_cells), each = k)
  j <- as.vector(t(indices))
  nn <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(cell_names, cell_names)
  )
  nn <- methods::as(nn, "Graph")
  graphs <- list(nn = nn)
  if (isTRUE(compute.SNN)) {
    membership <- Matrix::sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(n_cells, n_cells),
      dimnames = list(cell_names, cell_names)
    )
    snn <- Matrix::tcrossprod(membership)
    snn@x <- snn@x / (2 * k - snn@x)
    snn@x[snn@x < prune.SNN] <- 0
    snn <- Matrix::drop0(snn)
    snn <- methods::as(snn, "generalMatrix")
    snn <- methods::as(snn, "Graph")
    graphs[["snn"]] <- snn
  }
  graphs
}

#' @export
FindNeighbors.default <- function(
  object,
  query = NULL,
  distance.matrix = FALSE,
  k.param = 20,
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1 / 15,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  l2.norm = FALSE,
  cache.index = FALSE,
  index = NULL,
  cores = 1L,
  ...
) {
  extra <- list(...)
  fallback <- function() {
    find_neighbors_default_seurat(
      object = object,
      query = query,
      distance.matrix = distance.matrix,
      k.param = k.param,
      return.neighbor = return.neighbor,
      compute.SNN = compute.SNN,
      prune.SNN = prune.SNN,
      nn.method = nn.method,
      n.trees = n.trees,
      annoy.metric = annoy.metric,
      nn.eps = nn.eps,
      verbose = verbose,
      l2.norm = l2.norm,
      cache.index = cache.index,
      index = index,
      ...
    )
  }

  if (
    length(extra) || !is.matrix(object) || !is.numeric(object) ||
      !is.logical(distance.matrix) || length(distance.matrix) != 1L ||
      is.na(distance.matrix) ||
      !is.logical(return.neighbor) || length(return.neighbor) != 1L ||
      is.na(return.neighbor) ||
      !is.logical(compute.SNN) || length(compute.SNN) != 1L ||
      is.na(compute.SNN) ||
      !is.logical(l2.norm) || length(l2.norm) != 1L || is.na(l2.norm) ||
      !is.logical(cache.index) || length(cache.index) != 1L ||
      is.na(cache.index) || cache.index || !is.null(index) ||
      !is.numeric(prune.SNN) || length(prune.SNN) != 1L ||
      is.na(prune.SNN) ||
      length(nn.method) != 1L || is.na(nn.method) ||
      length(annoy.metric) != 1L || is.na(annoy.metric)
  ) {
    return(fallback())
  }

  n_reference <- nrow(object)
  k <- suppressWarnings(as.integer(k.param))
  if (
    length(k) != 1L || is.na(k) || k < 1L || n_reference < 2L ||
      !isTRUE(all.equal(as.numeric(k.param), as.numeric(k))) ||
      is.null(rownames(object))
  ) {
    return(fallback())
  }
  k_adjusted <- n_reference < k
  if (n_reference < k) {
    k <- n_reference - 1L
  }

  if (isTRUE(distance.matrix)) {
    if (!is.null(query) || isTRUE(l2.norm) || isTRUE(return.neighbor)) {
      return(fallback())
    }
    result <- find_neighbors_native_distance(object = object, k = k)
    if (is.null(result)) {
      return(fallback())
    }
    if (isTRUE(k_adjusted)) {
      warning(
        "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
        call. = FALSE
      )
    }
    graphs <- find_neighbors_as_graphs(
      indices = result[["idx"]],
      cell_names = rownames(object),
      compute.SNN = compute.SNN,
      prune.SNN = prune.SNN
    )
    return(graphs %||% fallback())
  }

  data_use <- object
  query_use <- query
  if (isTRUE(l2.norm)) {
    data_use <- find_neighbors_l2_rows(data_use)
    if (is.null(data_use)) {
      return(fallback())
    }
    if (!is.null(query_use)) {
      query_use <- find_neighbors_l2_rows(query_use)
      if (is.null(query_use)) {
        return(fallback())
      }
    }
  }

  need_distance <- isTRUE(return.neighbor) || !is.null(query)
  result <- find_neighbors_native_embedding(
    object = data_use,
    query = query_use,
    k = k,
    nn.method = nn.method,
    n.trees = n.trees,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    cores = cores,
    need_distance = need_distance
  )
  if (is.null(result)) {
    return(fallback())
  }
  if (isTRUE(k_adjusted)) {
    warning(
      "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
      call. = FALSE
    )
  }

  if (isTRUE(need_distance)) {
    return(find_neighbors_as_neighbor(
      idx = result[["idx"]],
      distance = result[["distance"]],
      cell_names = result[["cell.names"]],
      metric = result[["metric"]],
      ndim = result[["ndim"]]
    ))
  }

  graphs <- find_neighbors_as_graphs(
    indices = result[["idx"]],
    cell_names = result[["cell.names"]],
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN
  )
  graphs %||% fallback()
}

#' @export
FindNeighbors.dist <- function(
  object,
  k.param = 20,
  return.neighbor = FALSE,
  compute.SNN = !return.neighbor,
  prune.SNN = 1 / 15,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  l2.norm = FALSE,
  cache.index = FALSE,
  cores = 1L,
  ...
) {
  FindNeighbors.default(
    object = as.matrix(object),
    distance.matrix = TRUE,
    k.param = k.param,
    return.neighbor = return.neighbor,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.method = nn.method,
    n.trees = n.trees,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    verbose = verbose,
    l2.norm = l2.norm,
    cache.index = cache.index,
    cores = cores,
    ...
  )
}

#' @export
FindNeighbors.Assay <- function(
  object,
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
  l2.norm = FALSE,
  cache.index = FALSE,
  cores = 1L,
  ...
) {
  features <- features %||% SeuratObject::VariableFeatures(object)
  data.use <- Matrix::t(SeuratObject::GetAssayData(
    object = object,
    layer = "data"
  )[features, , drop = FALSE])
  FindNeighbors.default(
    object = as.matrix(data.use),
    k.param = k.param,
    return.neighbor = return.neighbor,
    compute.SNN = compute.SNN,
    prune.SNN = prune.SNN,
    nn.method = nn.method,
    n.trees = n.trees,
    annoy.metric = annoy.metric,
    nn.eps = nn.eps,
    verbose = verbose,
    l2.norm = l2.norm,
    cache.index = cache.index,
    cores = cores,
    ...
  )
}

# SeuratObject's v5 assay class is not an S3 subclass of `Assay`, so it needs
# an explicit method even though both assay generations share the same data
# extraction contract here.
#' @export
FindNeighbors.Assay5 <- FindNeighbors.Assay

knn_cross_topk_native <- function(
  reference,
  query,
  k,
  distance_metric,
  cores = 1L
) {
  if (
    !is.matrix(reference) || !is.matrix(query) ||
      !is.numeric(reference) || !is.numeric(query) ||
      ncol(reference) != ncol(query) || nrow(reference) == 0L ||
      nrow(query) == 0L || k < 1L || k > nrow(reference) ||
      !distance_metric %in% c("cosine", "euclidean", "pearson", "spearman") ||
      !is.numeric(cores) || length(cores) != 1L || is.na(cores) || cores < 1L
  ) {
    return(NULL)
  }
  if (distance_metric %in% c("pearson", "spearman")) {
    if (any(!is.finite(reference)) || any(!is.finite(query))) {
      return(NULL)
    }
    if (identical(distance_metric, "spearman")) {
      reference <- knn_rank_rows(reference)
      query <- knn_rank_rows(query)
    }
    reference <- reference - rowMeans(reference)
    query <- query - rowMeans(query)
    if (any(rowSums(reference * reference) <= 0) ||
      any(rowSums(query * query) <= 0)) {
      return(NULL)
    }
    distance_metric <- "cosine"
  }
  tryCatch(
    cross_knn_f32(
      reference = reference,
      query = query,
      k = as.integer(k),
      metric = distance_metric,
      cores = as.integer(cores)
    ),
    error = function(e) NULL
  )
}
