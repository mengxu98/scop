make_network_sample_guard_object <- function() {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(
    matrix(1, 2, 4, dimnames = list(c("g1", "g2"), letters[1:4])), sparse = TRUE))
  object$col <- c(0, 1, 0, 1)
  object$row <- c(0, 0, 0, 0)
  object$sample <- c("S1", "S1", "S2", "S2")
  object
}

test_that("spatial network sample guard rejects mixed samples before graph construction", {
  object <- make_network_sample_guard_object()
  calls <- 0L
  testthat::local_mocked_bindings(spatial_graph_compute = function(...) {
    calls <<- calls + 1L
    stop("graph reached")
  }, .package = "scop")
  run <- function(srt = object, sample.by = "sample", ...) {
    RunSpatialNetwork(srt, k = 1, sample.by = sample.by, verbose = FALSE, ...)
  }
  expect_error(run(), "requires one sample")
  for (column in list("missing", "", NA_character_, c("sample", "other"), 1L)) {
    expect_error(run(sample.by = column), "sample.by")
  }
  for (missing_id in c(NA_character_, "")) {
    bad <- object
    bad$sample[2] <- missing_id
    expect_error(run(bad), "non-missing, non-empty")
  }
  expect_identical(calls, 0L)
  expect_error(run(sample.by = NULL), "graph reached")
  single <- object[, 1:2]
  single$sample <- factor(c("S1", "S1"), levels = c("S1", "unused"))
  expect_error(run(single), "graph reached")
  expect_identical(calls, 2L)
})

test_that("spatial network sample guard uses selected image nodes and stores its field", {
  object <- make_network_sample_guard_object()
  object[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 0), row.names = c("a", "b")),
    type = "centroids", assay = SeuratObject::DefaultAssay(object), key = "guard_")
  testthat::local_mocked_bindings(spatial_graph_compute = function(coords, method, k, ...) {
    expect_identical(coords$cell_id, c("a", "b"))
    list(nodes = coords,
      edges = data.frame(from = 1L, to = 2L, distance = 1, weight = 1),
      parameters = list(method = method, k = k))
  }, .package = "scop")
  out <- RunSpatialNetwork(object, image = "slice1", k = 1,
    sample.by = "sample", verbose = FALSE)
  graph <- out@tools$SpatialNetwork$graphs[[out@tools$SpatialNetwork$active_graph]]
  expect_identical(graph$nodes$cell_id, c("a", "b"))
  expect_identical(graph$edges$from, "a")
  expect_identical(graph$edges$to, "b")
  expect_identical(graph$parameters$sample.by, "sample")
  expect_identical(out@tools$SpatialNetwork$parameters$sample.by, "sample")
})
