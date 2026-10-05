test_that("neighborhood significance is recomputed for the requested threshold", {
  df <- data.frame(from = c("MG", "OL", "MG", "OL"), to = c("OL", "MG", "MG", "OL"),
    estimate = c(1, -1, 1, 1), FDR = c(.04, .04, NA, NA),
    direction = c("enriched", "depleted", "enriched", "observed"))
  out <- scop:::spatial_neighborhood_prepare_plot_table(df, "estimate", .01)
  expect_identical(out$direction, c("ns", "ns", "ns", "observed"))
  expect_identical(scop:::spatial_neighborhood_prepare_plot_table(df, "estimate", .05)$direction,
    c("enriched", "depleted", "ns", "observed"))
})

test_that("multiple sample neighborhoods require explicit selection without aggregation", {
  object <- SeuratObject::CreateSeuratObject(matrix(1, 2, 4,
    dimnames = list(c("g1", "g2"), letters[1:4])))
  object$col <- c(0, 1, 0, 1); object$row <- 0
  object$label <- c("MG", "OL", "MG", "OL")
  object$sample <- c("S1", "S1", "S2", "S2")
  object <- RunSpatialNeighborhood(object, "label", sample.by = "sample",
    k = 1, backend = "r", verbose = FALSE)
  for (type in c("heatmap", "stat", "network")) {
    expect_error(SpatialNeighborhoodPlot(object, plot_type = type), "sample")
    p <- SpatialNeighborhoodPlot(object, plot_type = type, sample = "S1")
    expect_s3_class(p, "ggplot")
  }
  p <- SpatialNeighborhoodPlot(object, sample = "S1")
  expect_identical(unique(p$data$sample), "S1")
  expect_equal(p$data$fraction, c(.5, .5))
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2)
  expect_error(SpatialNeighborhoodPlot(object, sample = "typo"), "sample")
  expect_error(SpatialNeighborhoodPlot(object, sample = c("S1", "S2")), "sample")
  single <- RunSpatialNeighborhood(object[, 1:2], "label", sample.by = "sample",
    k = 1, backend = "r", verbose = FALSE)
  expect_s3_class(SpatialNeighborhoodPlot(single), "ggplot")
  expect_error(SpatialNeighborhoodPlot(object, plot_type = "spatial", sample = "S1"), "sample")
})

test_that("comparison selection prevents repeated pairs in differential plots", {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(1, 2, 2,
    dimnames = list(c("g1", "g2"), c("a", "b"))), sparse = TRUE))
  table <- data.frame(from = c("MG", "MG"), to = c("OL", "OL"),
    comparison = c("A", "B"), condition = c("A", "B"), sample = NA_character_,
    estimate = c(1, -1), FDR = c(.04, .04), direction = c("enriched", "depleted"))
  bundle <- scop:::spatial_tag_coordinate_contract(list(pair_table = table))
  object@tools$SpatialNeighborhood <- list(active_method = "spicyR", methods = list(spicyR = bundle))
  expect_error(SpatialNeighborhoodPlot(object), "comparison")
  p <- SpatialNeighborhoodPlot(object, plot_type = "stat", comparison = "A", FDR_threshold = .01)
  expect_identical(p$data$direction, "ns")
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$y, 1)
  expect_s3_class(SpatialNeighborhoodPlot(object, comparison = "absent"), "ggplot")
})
