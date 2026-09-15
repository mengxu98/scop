plot_contract_object <- function() {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(
    matrix(1:24, 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:8))), sparse = TRUE
  ))
  object$x <- rep(c(0, 10, 0, 10), 2)
  object$y <- rep(c(0, 0, 1, 1), 2)
  object$type <- factor(rep(c("B", "A"), 4), levels = c("A", "B"))
  object[["sample label"]] <- rep(c("S1", "S2"), each = 4)
  object
}

test_that("spatial panel geometry follows coordinates and can be overridden", {
  boundaries <- data.frame(cell_id = "cell", x = c(0, 10, 10, 0), y = c(0, 0, 1, 1))
  p <- SpatialCellPlot(boundaries = boundaries)
  gt <- ggplot2::ggplotGrob(p)
  panel <- gt$layout[gt$layout$name == "panel", ]
  expect_equal(as.numeric(gt$heights[panel$t]) / as.numeric(gt$widths[panel$l]), 0.1)
  expect_equal(SpatialCellPlot(boundaries = boundaries, theme_args = list(aspect.ratio = 1))$theme$aspect.ratio, 1)
})

test_that("spatial categories retain factor order and named colors across subsets", {
  object <- plot_contract_object()
  colors <- c(A = "#112233", B = "#445566")
  for (cells in list(colnames(object), colnames(object)[c(4, 3, 1)], colnames(object)[1])) {
    p <- SpatialSpotPlot(object, group.by = "type", cells = cells, palcolor = colors)
    expect_identical(levels(p$data$.value), c("A", "B"))
    scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")
    expect_identical(unname(scale$map(c("B", "A"))), unname(colors[c("B", "A")]))
  }
  expect_identical(spatial_palette_colors(c("B", "A"), palcolor = colors), colors[c("B", "A")])
  expect_error(SpatialSpotPlot(object, group.by = "type", palcolor = c(A = "red")), "missing colors")
  values <- setNames(factor(c("A", NA, "NA", rep("A", 5)), levels = c("A", "NA")), colnames(object))
  p <- SpatialSpotPlot(object, values = values, show_na = TRUE, palcolor = c(A = "red", `NA` = "blue"))
  expect_true(is.na(p$data$.value[2]))
  expect_identical(as.character(p$data$.value[3]), "NA")
  expect_identical(ggplot2::ggplot_build(p)$data[[1]]$fill[2:3], c("grey80", "blue"))
})

test_that("spatial facets accept metadata names with spaces in point long and pie plots", {
  object <- plot_contract_object()
  p <- SpatialSpotPlot(object, group.by = "type", split.by = "sample label")
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2)
  long <- data.frame(spot = colnames(object), label = object$type)
  p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "label",
    split.by = "sample label", legend.position = "bottom")
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2)
  expect_identical(p$theme$legend.position, "bottom")
  skip_if_not_installed("scatterpie")
  values <- matrix(rep(c(.2, .8), each = 8), 8, dimnames = list(colnames(object), c("A", "B")))
  p <- SpatialSpotPlot(object, values = values, plot_type = "pie", split.by = "sample label", legend.position = "bottom")
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2)
  expect_identical(p$theme$legend.position, "bottom")
})

test_that("constant spatial features respect explicit limits and ambiguous value IDs fail", {
  object <- plot_contract_object()
  values <- setNames(rep(0, ncol(object)), colnames(object))
  p <- SpatialSpotPlot(object, values = values, lower_cutoff = 0, upper_cutoff = 1)
  expect_equal(p$scales$get_scales("colour")$limits, c(0, 1))
  expect_error(SpatialSpotPlot(object, values = values, lower_cutoff = 2, upper_cutoff = 1), "lower continuous limit")
  expect_error(SpatialSpotPlot(object, values = values, upper_quantile = 2), "Quantiles")
  names(values)[2] <- names(values)[1]
  expect_error(SpatialSpotPlot(object, values = values), "unique non-empty spot names")
  mat <- matrix(1, ncol(object), 2, dimnames = list(rep("c1", ncol(object)), c("A", "B")))
  expect_error(SpatialSpotPlot(object, values = mat), "unique non-empty spot")
})

test_that("deconvolution empty panels and constant maps remain composable", {
  object <- plot_contract_object()
  object@tools$Custom <- list(proportions = matrix(c(rep(NA_real_, 8), rep(0, 8), seq(0, 1, length.out = 8)),
    8, dimnames = list(colnames(object), c("Missing", "Absent", "Measured"))))
  p <- SpatialDeconvolutionPlot(object, "Custom", combine = FALSE)
  expect_named(p, c("Missing", "Absent", "Measured"))
  expect_match(p$Missing$data$label, "No values")
  expect_equal(p$Absent$scales$get_scales("colour")$limits, c(0, 1))
  expect_s3_class(p$Absent$scales$get_scales("colour")$breaks, "waiver")
  expect_s3_class(SpatialDeconvolutionPlot(object, "Custom"), "patchwork")
  expect_no_error(ggplot2::ggplotGrob(p$Measured))
  object@tools$Custom$proportions[1, 3] <- Inf
  expect_error(SpatialDeconvolutionPlot(object, "Custom"), "infinite")
})

test_that("cell polygons support named lists layout and explicit legend controls", {
  boundaries <- data.frame(cell_id = rep(c("a", "b"), each = 4),
    x = c(0, 1, 1, 0, 2, 3, 3, 2), y = rep(c(0, 0, 1, 1), 2),
    gene1 = rep(c(0, 1), each = 4), gene2 = 0)
  p <- SpatialCellPlot(boundaries = boundaries, features = c("gene1", "gene2"),
    combine = FALSE, legend.position = "bottom", legend.title = "Expression")
  expect_named(p, c("gene1", "gene2"))
  expect_identical(p$gene1$theme$legend.position, "bottom")
  expect_identical(p$gene1$labels$fill, "Expression")
  expect_no_error(ggplot2::ggplotGrob(p$gene2))
  expect_s3_class(SpatialCellPlot(boundaries = boundaries, features = c("gene1", "gene2"), ncol = 1), "patchwork")
})

test_that("unassigned deconvolution spots are excluded from dominant summaries", {
  weights <- matrix(c(0, 0, .2, .8), 2, byrow = TRUE,
    dimnames = list(c("unassigned", "assigned"), c("A", "B")))
  finalized <- spatial_finalize_weights(weights, rownames(weights))
  summary <- spatial_weight_summary(finalized$weights)
  expect_identical(unname(finalized$dominant), c(NA_character_, "B"))
  expect_identical(summary$dominant_counts$type, "B")
  expect_equal(summary$dominant_counts$count, 1)
  expect_equal(c(summary$n_assigned, summary$n_unassigned), c(1, 1))
  empty <- spatial_weight_summary(weights[1, , drop = FALSE])
  expect_equal(nrow(empty$dominant_counts), 0)
  expect_true(all(is.na(empty$max_prop)))
  for (bad in c(NA_real_, Inf, -1)) {
    invalid <- weights
    invalid[2, 1] <- bad
    expect_error(spatial_finalize_weights(invalid, rownames(invalid)), "finite|non-negative")
  }
  weights[1, 1] <- -.Machine$double.eps
  expect_equal(spatial_normalize_weights(weights)[1, ], c(A = 0, B = 0))
})

test_that("invalid RCTD backend weights fail before committing object changes", {
  object <- plot_contract_object()
  before <- object
  invalid <- NA_real_
  testthat::local_mocked_bindings(rctd_run_spacexr = function(st_counts, ...) {
    weights <- matrix(.5, ncol(st_counts), 2,
      dimnames = list(colnames(st_counts), c("A", "B")))
    weights[1, 1] <- invalid
    list(weights = weights, metadata = list(), api = "test", object = list())
  })
  for (value in c(NA_real_, Inf, -1)) {
    invalid <- value
    expect_error(RunRCTD(object, reference = object, reference_label = "type", min_cells = 1, verbose = FALSE),
      "finite|non-negative")
    expect_identical(object, before)
  }
})

integration_plot_contract_object <- function() {
  object <- plot_contract_object()
  for (i in 1:2) {
    cells <- colnames(object)[object[["sample label", drop = TRUE]] == paste0("S", i)]
    object[[paste0("slice", i)]] <- SeuratObject::CreateFOV(
      data.frame(x = c(10, 20, 30, 40), y = c(30, 40, 60, 50), row.names = cells),
      type = "centroids", assay = "RNA", key = paste0("v", i, "_"))
  }
  object$aligned_x <- seq_len(8) + 100
  object$aligned_y <- seq_len(8) / 10
  parameters <- list(sample.by = "sample label", cluster_colname = "type",
    coord.cols = c("x", "y"), aligned_coord_cols = c("aligned_x", "aligned_y"), coordinate_space = "raw")
  object@tools$SpatialIntegration <- list(active_method = "PRECAST", parameters = parameters,
    methods = list(PRECAST = spatial_tag_coordinate_contract(list(parameters = parameters))))
  object
}

test_that("integration image maps retain all samples and consistent named colors", {
  object <- integration_plot_contract_object()
  before <- object
  colors <- c(A = "#112233", B = "#445566")
  plots <- SpatialIntegrationPlot(object, combine = FALSE, palcolor = colors)
  expect_named(plots, c("S1", "S2"))
  expect_setequal(unlist(lapply(plots, function(p) rownames(p$data))), colnames(object))
  for (p in plots) {
    expect_equal(p$data$x, c(10, 20, 30, 40))
    expect_identical(ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")$map(c("B", "A")), unname(colors[c("B", "A")]))
  }
  expect_s3_class(SpatialIntegrationPlot(object, ncol = 1), "patchwork")
  numeric_plots <- SpatialIntegrationPlot(object, group.by = "aligned_x", combine = FALSE)
  expect_true(all(vapply(numeric_plots, function(p) is.numeric(p$data$.value), logical(1))))
  expect_equal(numeric_plots$S1$scales$get_scales("colour")$limits, numeric_plots$S2$scales$get_scales("colour")$limits)
  expect_error(SpatialIntegrationPlot(object, image = "slice1"), "covering image")
  expect_error(SpatialIntegrationPlot(object, image = c(S1 = "slice1")), "cover every sample")
  expect_identical(object, before)
})

test_that("aligned integration views bypass raw images without changing the object", {
  object <- integration_plot_contract_object()
  before <- object
  p <- SpatialIntegrationPlot(object, use_aligned = TRUE)
  expect_equal(unname(p$data$x), unname(object$aligned_x))
  expect_equal(unname(p$data$y), unname(object$aligned_y))
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2)
  expect_error(SpatialIntegrationPlot(object, use_aligned = TRUE, overlay_image = TRUE), "raw image overlay")
  alignment <- SpatialIntegrationPlot(object, plot_type = "alignment")
  expect_equal(alignment$data$x[alignment$data$coordinate == "Raw"], rep(c(10, 20, 30, 40), 2))
  expect_identical(object, before)
})
