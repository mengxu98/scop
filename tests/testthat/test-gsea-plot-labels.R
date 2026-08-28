test_that("GSEAPlot uses the correctly spelled negative correlation label", {
  gsea_plot_body <- paste(deparse(body(GSEAPlot)), collapse = "\n")

  expect_match(gsea_plot_body, "Negatively correlated", fixed = TRUE)
  expect_false(grepl("Negtively correlated", gsea_plot_body, fixed = TRUE))
})
