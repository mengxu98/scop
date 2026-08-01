test_that("ScissorPlot maps the scop theme for thisplot statistical panels", {
  data(pancreas_sub, package = "scop")
  n <- ncol(pancreas_sub)
  pancreas_sub$Scissor_coef <- seq(-1, 1, length.out = n)
  pancreas_sub$Scissor_status <- rep(
    c("Scissor+", "Scissor-", "Background"),
    length.out = n
  )

  p <- ScissorPlot(
    pancreas_sub,
    plot_type = "bar",
    group.by = "CellType"
  )

  expect_true(inherits(p, c("ggplot", "patchwork")))
})
