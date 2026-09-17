followup_spatial_object <- function() {
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(
    matrix(1:12, 3, dimnames = list(paste0("g", 1:3), letters[1:4])), sparse = TRUE
  ))
  object$x <- c(0, 1, 0, 1)
  object$y <- c(0, 0, 1, 1)
  object$type <- factor(c("B", "A", "B", "A"), levels = c("A", "B"))
  object[["sample label"]] <- c("S1", "S1", "S2", "S2")
  object$Cell2location_dominant_type <- object$type
  object@tools$Cell2location <- list(abundance = matrix(
    c(1, 2, 40, 100, 2, 3, 4, 5), 4,
    dimnames = list(colnames(object), c("A", "B"))
  ))
  object
}

test_that("matrix point plots honor explicit scales in both return forms", {
  object <- followup_spatial_object()
  for (combine in c(FALSE, TRUE)) {
    p <- Cell2locationPlot(object, plot_type = "abundance", upper_cutoff = 10, combine = combine)
    expect_equal(p[[1]]$scales$get_scales("colour")$limits, c(1, 10))
    expect_equal(p[[2]]$scales$get_scales("colour")$limits, c(2, 10))
    q <- Cell2locationPlot(object, plot_type = "abundance", lower_cutoff = 0,
      upper_cutoff = 10, combine = combine)
    expect_equal(q[[1]]$scales$get_scales("colour")$limits, c(0, 10))
    expect_equal(q[[2]]$scales$get_scales("colour")$limits, c(0, 10))
  }
  q <- Cell2locationPlot(object, plot_type = "abundance", upper_quantile = .5, combine = FALSE)
  expect_equal(q[[1]]$scales$get_scales("colour")$limits, c(1, 21))
  p <- Cell2locationPlot(object, plot_type = "abundance", combine = FALSE)
  expect_equal(p[[1]]$scales$get_scales("colour")$limits, c(0, 100))
  expect_equal(p[[2]]$scales$get_scales("colour")$limits, c(0, 100))
})

test_that("dominant maps retain the named-list return contract", {
  object <- followup_spatial_object()
  p <- Cell2locationPlot(object, plot_type = "dominant", combine = FALSE, split.by = "sample label")
  expect_type(p, "list")
  expect_named(p, "Cell2location_dominant_type")
  expect_equal(nrow(ggplot2::ggplot_build(p[[1]])$layout$layout), 2)
  expect_s3_class(Cell2locationPlot(object, plot_type = "dominant"), "ggplot")
})

test_that("long points share geometry and named colors with ordinary points", {
  object <- followup_spatial_object()
  colors <- c(A = "#4477AA", B = "#EE6677")
  points <- SpatialSpotPlot(object, group.by = "type", palcolor = colors)
  long <- data.frame(spot = colnames(object), label = object$type)
  p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "label", palcolor = colors)
  point_data <- ggplot2::ggplot_build(points)$data[[1]]
  long_data <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(long_data[, c("x", "y", "fill")], point_data[, c("x", "y", "fill")])
  expect_null(p$labels$title)
  long <- long[c(1, 1, 3, 2), ]
  p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "label",
    split.by = "sample label", palcolor = colors, cells = c("a", "c"))
  expect_equal(unname(p$data$x), unname(object$x[c(1, 1, 3)]))
  expect_equal(nrow(p$data), 3)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2)
})

test_that("numeric long plots retain jitter and original color values", {
  object <- followup_spatial_object()
  long <- data.frame(spot = colnames(object), x = c(10, 20, 30, 40))
  p <- SpatialSpotPlot(object, plot.data = long, spot.by = "spot", color.by = "x",
    geom = "jitter", jitter_width = .2, jitter_height = .1)
  expect_equal(p$data$.value, long$x)
  set.seed(1)
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_true(all(abs(built$x - object$x) <= .2))
  expect_true(all(abs(built$y - object$y) <= .1))
  expect_equal(p$scales$get_scales("colour")$limits, c(10, 40))
})

test_that("constant topic panels share default breaks and keep empty panels", {
  object <- followup_spatial_object()
  object@tools$STdeconvolve <- list(theta = matrix(c(rep(0, 4), rep(.5, 4), rep(NA_real_, 4)), 4,
    dimnames = list(colnames(object), c("Absent", "Constant", "Missing"))),
    parameters = list(prefix = "STdeconvolve"))
  p <- STdeconvolvePlot(object, combine = FALSE)
  expect_equal(p[[1]]$scales$get_scales("colour")$limits, c(0, .5))
  expect_equal(p[[2]]$scales$get_scales("colour")$limits, c(0, .5))
  expect_s3_class(p[[1]]$scales$get_scales("colour")$breaks, "waiver")
  expect_s3_class(p[[2]]$scales$get_scales("colour")$breaks, "waiver")
  expect_match(p[[3]]$data$label, "No values")
})

test_that("documented bridge counts and merged image auto-resolution work", {
  skip_if_not_installed("SpatialExperiment")
  data(visium_human_pancreas_sub, package = "scop")
  object <- suppressWarnings(visium_human_pancreas_sub[, 1:16])
  spe <- srt_to_spe(object, assay = "Spatial", layer = "counts", image = "slice1")
  restored <- spe_to_srt(spe, assay = "Spatial", layer = "scop_input")
  p <- SpatialSpotPlot(restored, features = rownames(restored)[1], assay = "Spatial", layer = "counts")
  expect_equal(p$data$.value, as.numeric(GetAssayData5(object, assay = "Spatial", layer = "counts")[1, ]))
  testthat::local_mocked_bindings(spatial_integration_run_backend = function(method, input, ...) {
    list(domains = stats::setNames(rep("D1", length(input$cells)), input$cells), raw_result = list())
  })
  integrated <- suppressWarnings(RunSpatialIntegration(
    list(S1 = object[, 1:8], S2 = object[, 9:16]), assay = "Spatial", layer = "counts",
    image = c(S1 = "slice1", S2 = "slice1"), store_object = FALSE, verbose = FALSE
  ))
  plots <- SpatialIntegrationPlot(integrated, combine = FALSE)
  expect_named(plots, c("S1", "S2"))
  expect_setequal(unlist(lapply(plots, function(p) rownames(p$data))), colnames(integrated))
})
