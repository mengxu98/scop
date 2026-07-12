test_that("spatial graph core handles sparse edge cases", {
  skip_if_not_installed("BiocNeighbors")
  collinear <- data.frame(cell_id = letters[1:5], x = 1:5, y = 0)
  graph <- scop:::spatial_graph_compute(collinear, k = 2)
  expect_true(all(graph$edges$from < graph$edges$to))
  expect_equal(anyDuplicated(paste(graph$edges$from, graph$edges$to)), 0L)
  expect_true(all(graph$edges$distance > 0))

  duplicate <- data.frame(cell_id = letters[1:4], x = c(0, 0, 1, 2), y = 0)
  duplicate_graph <- scop:::spatial_graph_compute(duplicate, k = 1)
  expect_true(any(duplicate_graph$edges$distance == 0))

  empty <- scop:::spatial_graph_compute(collinear, method = "radius", radius = 0.1)
  expect_equal(nrow(empty$nodes), 5L)
  expect_equal(nrow(empty$edges), 0L)
  expect_error(scop:::spatial_graph_compute(collinear[1, ], k = 1), "At least two")
  expect_error(scop:::spatial_graph_compute(collinear, k = 5), "number of nodes minus one")
})

test_that("spatial graph weights preserve raw distance", {
  skip_if_not_installed("BiocNeighbors")
  coords <- data.frame(cell_id = letters[1:3], x = c(0, 3, 0), y = c(0, 0, 4))
  binary <- scop:::spatial_graph_compute(coords, k = 1, directed = TRUE)
  inverse <- scop:::spatial_graph_compute(coords, k = 1, directed = TRUE, weight = "inverse_distance")
  gaussian <- scop:::spatial_graph_compute(coords, k = 1, directed = TRUE, weight = "gaussian", sigma = 2)
  expect_equal(binary$edges$distance, inverse$edges$distance)
  expect_equal(inverse$edges$weight, 1 / (1 + inverse$edges$distance))
  expect_equal(gaussian$edges$weight, exp(-(gaussian$edges$distance^2) / 8))
})

test_that("raw and display coordinate transforms round trip", {
  raw <- data.frame(cell_id = c("a", "b"), x = c(2, 8), y = c(4, 10))
  transform <- list(scale = 0.5, y_flip = TRUE, image_height = 20)
  display <- scop:::spatial_coords_to_display(raw, transform)
  expect_equal(display$x, c(1, 4))
  expect_equal(display$y, c(18, 15))
  expect_equal(scop:::spatial_coords_to_raw(display, transform), raw)
})

test_that("RunSpatialNetwork uses deterministic slots and collision protection", {
  skip_if_not_installed("BiocNeighbors")
  counts <- matrix(
    seq_len(12), nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(0, 1, 0, 1)
  srt$row <- c(0, 0, 1, 1)
  out <- RunSpatialNetwork(srt, k = 2, verbose = FALSE)
  expect_identical(out@tools$SpatialNetwork$active_graph, "knn_k2")
  expect_identical(out@tools$SpatialNetwork$graphs$knn_k2$source$coordinate_space, "raw")
  expect_error(RunSpatialNetwork(out, k = 2, verbose = FALSE), "already exists")
  out <- RunSpatialNetwork(out, k = 2, overwrite = TRUE, verbose = FALSE)
  out <- RunSpatialNetwork(out, method = "radius", radius = 1.5, verbose = FALSE)
  expect_identical(out@tools$SpatialNetwork$active_graph, "radius_r1p5")
  expect_true(all(c("knn_k2", "radius_r1p5") %in% names(out@tools$SpatialNetwork$graphs)))
})

test_that("dedicated spatial result plots accept result-only input", {
  misty <- list(results = list(
    improvements = data.frame(target = c("A", "B"), measure = "gain.R2", value = c(0.1, 0.4)),
    contributions = data.frame(target = "A", view = "paraview.3", value = 0.3)
  ))
  expect_s3_class(MistyRPlot(res = misty), "ggplot")
  expect_s3_class(MistyRPlot(res = misty, type = "contributions"), "ggplot")

  statial <- list(table = data.frame(
    imageID = rep(c("s1", "s2"), each = 2), test = rep(c("A:B", "B:A"), 2),
    r = rep(c(10, 20), 2), kontextual = c(-1, 0.5, 0.2, 1)
  ))
  expect_s3_class(StatialKontextualPlot(res = statial), "ggplot")
})
