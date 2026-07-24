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
