test_that("real PRECAST receives selected counts and a nonempty distance graph", {
  skip_if_not_installed("PRECAST")
  skip_if_not_installed("BiocNeighbors")
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
  expect_identical(bundle$parameters$adjacency_builder, "scop::spatial_graph_compute")
  raw <- bundle$raw_result
  for (sample in seq_along(raw@seulist)) {
    native <- raw@seulist[[sample]]
    expect_identical(SeuratObject::DefaultAssay(native), "Spatial")
    expect_equal(as.matrix(GetAssayData5(native, layer = "counts")),
      counts[, colnames(native), drop = FALSE] * 2)
    graph <- raw@AdjList[[sample]]
    expect_equal(as.numeric(Matrix::colSums(graph)), rep(6, ncol(native)))
    # Independently verify every column against Euclidean nearest distances.
    xy <- cbind(native$col, native$row)
    correct <- vapply(seq_len(nrow(xy)), function(i) {
      d2 <- rowSums(sweep(xy, 2, xy[i, ], "-")^2)
      d2[i] <- Inf
      neighbors <- which(as.numeric(graph[, i]) != 0)
      all(d2[neighbors] <= sort(d2, partial = 6)[6] + 1e-8)
    }, logical(1))
    expect_true(all(correct))
  }
  expect_identical(object, before)
  # Long rows defeat PRECAST's native single-axis candidate restriction.
  # The exact graph must also be invariant to swapping axes and translation.
  long_input <- list(coords_list = lapply(raw@seulist, function(sample) {
    data.frame(cell_id = colnames(sample), x = seq_len(ncol(sample)) * 55,
      y = 0, row.names = colnames(sample))
  }))
  long_graph <- spatial_integration_set_precast_adjacency(raw, long_input,
    list(type = "fixed_number", number = 6))
  rotated <- long_input
  rotated$coords_list <- lapply(rotated$coords_list, function(coords) {
    coords$y <- coords$x + 10000
    coords$x <- 20000
    coords
  })
  rotated_graph <- spatial_integration_set_precast_adjacency(raw, rotated,
    list(type = "fixed_number", number = 6))
  expect_identical(long_graph@AdjList, rotated_graph@AdjList)
  for (graph in long_graph@AdjList) {
    correct <- vapply(seq_len(ncol(graph)), function(i) {
      d <- abs(seq_len(ncol(graph)) - i)
      d[i] <- Inf
      neighbors <- which(as.numeric(graph[, i]) != 0)
      length(neighbors) == 6 && all(d[neighbors] <= sort(d, partial = 6)[6])
    }, logical(1))
    expect_true(all(correct))
  }
  expect_error(spatial_integration_validate_neighbor_count(raw,
    list(type = "fixed_number", number = 60)), "fewer neighbors")
  broken <- raw
  broken@AdjList[[1]] <- raw@AdjList[[1]] * 0
  expect_error(spatial_integration_validate_adjacency(broken), "no edges")
  broken@AdjList[[1]] <- raw@AdjList[[1]] * NA_real_
  expect_error(spatial_integration_validate_adjacency(broken), "invalid adjacency")
})
