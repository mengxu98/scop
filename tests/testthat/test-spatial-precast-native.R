test_that("real PRECAST receives selected counts and a nonempty distance graph", {
  skip_if_not_installed("PRECAST")
  check_r("feiyoung/PRECAST", verbose = FALSE)
  set.seed(42)
  counts <- matrix(rpois(40 * 120, 6) + 1, 40,
    dimnames = list(paste0("g", 1:40), paste0("s", 1:120)))
  counts[1:20, rep(1:30, 2) + rep(c(0, 60), each = 30)] <-
    counts[1:20, rep(1:30, 2) + rep(c(0, 60), each = 30)] + 10
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), assay = "RNA")
  object[["Spatial"]] <- SeuratObject::CreateAssayObject(counts = Matrix::Matrix(counts * 2, sparse = TRUE))
  grid <- expand.grid(x = seq_len(10) * 55, y = seq_len(6) * 55)
  object$x <- rep(grid$x, 2)
  object$y <- rep(grid$y, 2)
  object$sample <- rep(c("A", "B"), each = 60)
  before <- object
  out <- RunSpatialIntegration(object, assay = "Spatial", sample.by = "sample",
    par_params = list(maxIter = 3, maxIter_ICM = 2, coreNum = 1, coreNum_int = 1,
      init.nstart = 1, int.model = "kmeans", beta_grid = c(.2, .5), verbose = FALSE),
    run_params = list(K = 2, q = 3), verbose = FALSE)
  bundle <- out@tools$SpatialIntegration$methods$PRECAST
  expect_true(all(is.finite(bundle$embedding)))
  expect_false(anyNA(bundle$domains))
  expect_setequal(names(bundle$domains), colnames(object))
  expect_identical(bundle$parameters$backend_parameters$adj_params$type, "fixed_number")
  expect_true(all(bundle$parameters$adjacency_summary > 0))
  raw <- bundle$raw_result
  for (sample in seq_along(raw@seulist)) {
    native <- raw@seulist[[sample]]
    expect_identical(SeuratObject::DefaultAssay(native), "Spatial")
    expect_equal(as.matrix(GetAssayData5(native, layer = "counts")),
      counts[, colnames(native), drop = FALSE] * 2)
    graph <- raw@AdjList[[sample]]
    expect_true(all(Matrix::rowSums(graph) > 0))
    # Independent physical assertion: with a 55-pixel lattice the default
    # nearest-neighbor graph must include first-neighbor connections.
    ij <- which(as.matrix(graph) != 0, arr.ind = TRUE)
    xy <- cbind(native$col, native$row)
    distance <- sqrt(rowSums((xy[ij[, 1], , drop = FALSE] - xy[ij[, 2], , drop = FALSE])^2))
    expect_true(any(abs(distance - 55) < 1e-8))
  }
  expect_identical(object, before)
  expect_error(spatial_integration_validate_neighbor_count(raw,
    list(type = "fixed_number", number = 60)), "fewer neighbors")
  broken <- raw
  broken@AdjList[[1]] <- raw@AdjList[[1]] * 0
  expect_error(spatial_integration_validate_adjacency(broken), "no edges")
  broken@AdjList[[1]] <- raw@AdjList[[1]] * NA_real_
  expect_error(spatial_integration_validate_adjacency(broken), "invalid adjacency")
})
