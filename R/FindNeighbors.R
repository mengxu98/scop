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
  backend = c("exact", "annoy"),
  cores = 1L,
  ...
) {
  if (
    !inherits(object, "Seurat") ||
      !identical(nn.method, "annoy") ||
      !identical(annoy.metric, "euclidean") ||
      !isFALSE(return.neighbor) ||
      !isFALSE(l2.norm) ||
      !is.null(features) ||
      is.null(dims) ||
      !isFALSE(cache.index) ||
      !isFALSE(do.plot) ||
      !isTRUE(nn.eps == 0)
  ) {
    stop("FindNeighbors.Seurat received unsupported arguments for the scop implementation.", call. = FALSE)
  }
  assay <- SeuratObject::DefaultAssay(object[[reduction]])
  data.use <- SeuratObject::Embeddings(object[[reduction]])[
    ,
    dims,
    drop = FALSE
  ]
  cell.names <- rownames(data.use)
  n.cells <- nrow(data.use)
  cores <- suppressWarnings(as.integer(cores[[1L]]))
  if (is.na(cores) || cores < 1L) {
    stop("cores must be a positive integer.", call. = FALSE)
  }
  backend <- match.arg(backend)
  if (identical(backend, "exact")) {
    nn.idx <- exact_knn_f32(data.use, k.param, cores)
  } else {
    nn.idx <- annoy_build_search(data.use, k.param, n.trees, cores)
  }

  j <- as.numeric(t(nn.idx))
  i <- rep(seq_len(n.cells), each = k.param)
  nn.matrix <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(n.cells, n.cells),
    dimnames = list(cell.names, cell.names)
  )
  nn.matrix <- methods::as(nn.matrix, "Graph")
  SeuratObject::DefaultAssay(nn.matrix) <- assay

  if (compute.SNN) {
    snn.matrix <- Seurat:::ComputeSNN(
      nn_ranked = nn.idx,
      prune = prune.SNN
    )
    rownames(snn.matrix) <- cell.names
    colnames(snn.matrix) <- cell.names
    snn.matrix <- SeuratObject::as.Graph(snn.matrix)
    SeuratObject::DefaultAssay(snn.matrix) <- assay
  }
  graph.name <- if (is.null(graph.name)) {
    paste0(assay, "_", c("nn", "snn"))
  } else {
    graph.name
  }
  object[[graph.name[1]]] <- nn.matrix
  if (compute.SNN && length(graph.name) >= 2) {
    object[[graph.name[2]]] <- snn.matrix
  }
  object
}

#' @export
FindNeighbors <- function(object, ...) {
  UseMethod("FindNeighbors")
}

#' @export
FindNeighbors.default <- function(object, ...) {
  stop("FindNeighbors supports Seurat objects.", call. = FALSE)
}
