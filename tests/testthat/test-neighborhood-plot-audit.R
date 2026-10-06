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
  make_saved_object <- function(samples = c("S1", "S2")) {
    n_cells <- 2L * length(samples)
    object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(1, 2, n_cells,
      dimnames = list(c("g1", "g2"), letters[seq_len(n_cells)])), sparse = TRUE))
    object$col <- rep(c(0, 1), length(samples)); object$row <- 0
    object$label <- rep(c("MG", "OL"), length(samples))
    object$sample <- rep(samples, each = 2L)

    # Each two-cell sample has MG -> OL and OL -> MG edges at k = 1.
    # Use the standardized saved pair-table contract so these plotting tests
    # do not require BiocNeighbors to recompute the neighborhood graph.
    pair_table <- data.frame(
      method = "observed", comparison = "all", condition = "all",
      from = rep(c("MG", "OL"), length(samples)),
      to = rep(c("OL", "MG"), length(samples)), estimate = .5,
      statistic = NA_real_, pval = NA_real_, FDR = NA_real_,
      direction = "observed", sample = rep(samples, each = 2L),
      subject = rep(samples, each = 2L), count = 1L, total = 2L, fraction = .5
    )
    bundle <- scop:::spatial_tag_coordinate_contract(list(
      method = "observed", pair_table = pair_table,
      parameters = list(method = "observed", coordinate_space = "raw",
        group.by = "label", sample.by = "sample", k = 1L)
    ))
    object@tools$SpatialNeighborhood <- scop:::spatial_tag_coordinate_contract(list(
      method = "SpatialNeighborhood", active_method = "observed",
      methods = list(observed = bundle), pair_table = pair_table,
      parameters = bundle$parameters
    ))
    object
  }

  object <- make_saved_object()
  single <- make_saved_object("S1")
  for (type in c("heatmap", "stat", "network")) {
    expect_error(SpatialNeighborhoodPlot(object, plot_type = type), "sample")
    for (sample in c("S1", "S2")) {
      p <- SpatialNeighborhoodPlot(object, plot_type = type, sample = sample)
      expect_s3_class(p, "ggplot")
      expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2)
      if (type != "network") {
        expect_identical(unique(p$data$sample), sample)
        expect_equal(p$data$fraction, c(.5, .5))
      } else {
        expect_equal(p$layers[[1]]$data$weight, c(.5, .5))
      }
    }
    p_single <- SpatialNeighborhoodPlot(single, plot_type = type)
    expect_s3_class(p_single, "ggplot")
    expect_equal(nrow(ggplot2::ggplot_build(p_single)$data[[1]]), 2)
  }
  p <- SpatialNeighborhoodPlot(object, sample = "S1")
  expect_identical(unique(p$data$sample), "S1")
  expect_equal(p$data$fraction, c(.5, .5))
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2)
  expect_error(SpatialNeighborhoodPlot(object, sample = "typo"), "sample")
  expect_error(SpatialNeighborhoodPlot(object, sample = c("S1", "S2")), "sample")
  for (sample in list(NA_character_, character(), "", 1)) {
    expect_error(SpatialNeighborhoodPlot(object, sample = sample), "sample")
  }
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
