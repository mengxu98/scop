make_native_neighbor_matrix <- function() {
  x <- rbind(
    a = c(0, 0),
    b = c(0.2, 0.1),
    c = c(4, 4),
    d = c(4.2, 4.1),
    e = c(8, 0)
  )
  storage.mode(x) <- "double"
  x
}

make_native_seurat_object <- function(seed = 71) {
  set.seed(seed)
  counts <- Matrix::rsparsematrix(80, 40, density = 0.12)
  counts@x <- abs(counts@x * 10)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 40, verbose = FALSE)
  obj <- Seurat::ScaleData(
    obj,
    features = SeuratObject::VariableFeatures(obj),
    verbose = FALSE
  )
  Seurat::RunPCA(
    obj,
    features = SeuratObject::VariableFeatures(obj),
    npcs = 10,
    verbose = FALSE
  )
}

test_that("FindNeighbors default delegates slower Annoy matrix queries", {
  reference <- make_native_neighbor_matrix()
  query <- rbind(q1 = c(0.1, 0.1), q2 = c(4.1, 4.1))
  calls <- 0L
  testthat::local_mocked_bindings(
    annoy_cross_knn = function(...) {
      calls <<- calls + 1L
    },
    .package = "scop"
  )

  out <- FindNeighbors(
    reference,
    query = query,
    k.param = 2,
    return.neighbor = TRUE,
    verbose = FALSE
  )

  expect_identical(calls, 0L)
  expect_s4_class(out, "Neighbor")
  expect_identical(dim(out@nn.idx), c(2L, 2L))
  expect_identical(out@cell.names, rownames(query))
  expect_identical(out@alg.info$metric, "euclidean")
  expect_true(all(is.finite(out@nn.dist)))
})

test_that("delegated cosine Neighbor distances match Seurat semantics", {
  reference <- rbind(
    a = c(1, 0),
    b = c(0.8, 0.2),
    c = c(0, 1),
    d = c(-1, 0)
  )
  query <- rbind(q1 = c(1, 0.1), q2 = c(0.1, 1))

  actual <- FindNeighbors(
    reference,
    query = query,
    k.param = 2,
    annoy.metric = "cosine",
    return.neighbor = TRUE,
    verbose = FALSE
  )
  expected <- seurat_reference_method(
    "FindNeighbors",
    "default",
    reference,
    query = query,
    k.param = 2,
    annoy.metric = "cosine",
    return.neighbor = TRUE,
    verbose = FALSE
  )

  expect_identical(actual@nn.idx, expected@nn.idx)
  expect_equal(actual@nn.dist, expected@nn.dist, tolerance = 1e-6)
})

test_that("delegated matrix flow builds Seurat-compatible NN and SNN graphs", {
  x <- make_native_neighbor_matrix()
  actual <- FindNeighbors(x, k.param = 3, verbose = FALSE)
  expected <- seurat_reference_method(
    "FindNeighbors",
    "default",
    x,
    k.param = 3,
    verbose = FALSE
  )

  expect_named(actual, c("nn", "snn"))
  expect_s4_class(actual$nn, "Graph")
  expect_s4_class(actual$snn, "Graph")
  expect_identical(as.matrix(actual$nn), as.matrix(expected$nn))
  expect_equal(as.matrix(actual$snn), as.matrix(expected$snn), tolerance = 1e-12)
})

test_that("precomputed distance matrices use the native graph flow", {
  x <- make_native_neighbor_matrix()
  distance <- as.matrix(stats::dist(x))
  testthat::local_mocked_bindings(
    find_neighbors_default_seurat = function(...) {
      stop("distance matrix unexpectedly fell back to Seurat")
    },
    .package = "scop"
  )

  actual <- FindNeighbors(
    distance,
    distance.matrix = TRUE,
    k.param = 3,
    verbose = FALSE
  )
  expected <- seurat_reference_method(
    "FindNeighbors",
    "default",
    distance,
    distance.matrix = TRUE,
    k.param = 3,
    verbose = FALSE
  )

  expect_identical(as.matrix(actual$nn), as.matrix(expected$nn))
  expect_equal(as.matrix(actual$snn), as.matrix(expected$snn), tolerance = 1e-12)
})

test_that("unsupported index and RANN requests fall back to Seurat", {
  x <- make_native_neighbor_matrix()
  fallback_calls <- 0L
  fallback <- get("find_neighbors_default_seurat", asNamespace("scop"))
  testthat::local_mocked_bindings(
    find_neighbors_default_seurat = function(...) {
      fallback_calls <<- fallback_calls + 1L
      fallback(...)
    },
    .package = "scop"
  )

  expect_s4_class(
    FindNeighbors(
      x,
      k.param = 2,
      nn.method = "rann",
      return.neighbor = TRUE,
      verbose = FALSE
    ),
    "Neighbor"
  )
  expect_s4_class(
    FindNeighbors(
      x,
      k.param = 2,
      cache.index = TRUE,
      return.neighbor = TRUE,
      verbose = FALSE
    ),
    "Neighbor"
  )
  expect_identical(fallback_calls, 2L)
})

test_that("Seurat method stores a delegated Neighbor result", {
  obj <- make_native_seurat_object()
  native_calls <- 0L
  testthat::local_mocked_bindings(
    annoy_cross_knn = function(...) {
      native_calls <<- native_calls + 1L
    },
    .package = "scop"
  )

  out <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:10,
    k.param = 10,
    return.neighbor = TRUE,
    verbose = FALSE
  )

  expect_identical(native_calls, 0L)
  expect_true("RNA.nn" %in% names(out@neighbors))
  expect_s4_class(out@neighbors$RNA.nn, "Neighbor")
})

test_that("dist and Assay inputs retain the complete Seurat workflow", {
  x <- make_native_neighbor_matrix()
  dist_out <- FindNeighbors(stats::dist(x), k.param = 3, verbose = FALSE)
  expect_named(dist_out, c("nn", "snn"))

  obj <- make_native_seurat_object()
  assay_out <- FindNeighbors(
    obj[["RNA"]],
    features = SeuratObject::VariableFeatures(obj),
    k.param = 5,
    return.neighbor = TRUE,
    verbose = FALSE
  )
  expect_s4_class(assay_out, "Neighbor")
})
