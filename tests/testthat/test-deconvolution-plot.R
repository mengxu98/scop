test_that("DeconvolutionPlot applies its major-grid helper without recursion", {
  res <- data.frame(
    sample = rep(c("Sample1", "Sample2"), each = 2L),
    cell_type = rep(c("Alpha", "Beta"), times = 2L),
    proportion = c(0.7, 0.3, 0.4, 0.6),
    method = "CIBERSORT",
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    DeconvolutionPlot(res = res, plot_type = "bar", verbose = FALSE),
    "ggplot"
  )
  expect_s3_class(
    DeconvolutionPlot(res = res, plot_type = "box", verbose = FALSE),
    "ggplot"
  )
})

test_that("DeconvolutionPlot keeps stored values and requested sample order", {
  res <- data.frame(
    sample = rep(c("Sample1", "Sample2"), each = 2L),
    cell_type = rep(c("Alpha", "Beta"), times = 2L),
    proportion = c(1.4, 0.6, 0.8, 1.2),
    method = "CIBERSORT"
  )
  bar <- DeconvolutionPlot(
    res = res, plot_type = "bar",
    sample_order = c("Sample2", "Sample1"),
    cell_type_order = c("Beta", "Alpha")
  )
  expect_equal(sort(bar$data$abundance), sort(res$proportion))
  expect_identical(unique(as.character(bar$data$sample)), c("Sample2", "Sample1"))
  expect_identical(unique(as.character(bar$data$cell_type)), c("Beta", "Alpha"))
  boxes <- DeconvolutionPlot(res = res, plot_type = "box")
  expect_equal(sort(boxes$data$value), sort(res$proportion))
  expect_no_error(ggplot2::ggplotGrob(boxes))
})

test_that("tAgePlot retains predictions, ordering and facet membership", {
  predictions <- data.frame(
    condition = rep(c("Young", "Aged"), each = 4),
    batch = rep(c("A", "B"), times = 4),
    Mouse_tAge = c(3.1, 4.2, 5.0, 4.6, 18.2, 20.1, 21.3, 19.6)
  )
  plot <- tAgePlot(predictions, group.by = "condition", split.by = "batch")
  expect_equal(sort(plot$data$value), sort(predictions$Mouse_tAge))
  expect_identical(levels(plot$data$group.by), c("Young", "Aged"))
  expect_identical(
    plot$data$batch,
    predictions$batch[match(plot$data$cells, rownames(predictions))]
  )
  expect_no_error(ggplot2::ggplotGrob(plot))
})

test_that("immune box plots retain abundance values and comparison groups", {
  abundance <- matrix(
    c(0.7, 0.3, 0.6, 0.4, 0.4, 0.6, 0.3, 0.7),
    ncol = 2, byrow = TRUE,
    dimnames = list(paste0("Sample", 1:4), c("T cells", "B cells"))
  )
  groups <- stats::setNames(c("Control", "Control", "Treated", "Treated"), rownames(abundance))
  plot <- ImmuneAbundancePlot(
    immune.data = abundance, group.data = groups,
    group.by = "condition", plot_type = "box"
  )
  expect_equal(sort(plot$data$value), sort(as.numeric(abundance)))
  expect_identical(levels(plot$data$split.by), c("Control", "Treated"))
  expect_no_error(ggplot2::ggplotGrob(plot))
})
