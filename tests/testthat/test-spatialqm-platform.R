test_that("RunSpatialQM forwards the selected platform to backend metrics", {
  skip_if_not_installed("SpatialQM")

  data(visium_human_pancreas_sub, package = "scop")
  out <- suppressWarnings(RunSpatialQM(
    visium_human_pancreas_sub,
    assay = "Spatial",
    layer = "counts",
    metrics = "n_cells",
    platform = "Visium",
    verbose = FALSE
  ))

  summary <- out@tools$SpatialQM$summary
  expect_equal(summary$metric, "n_cells")
  expect_equal(summary$status, "ok")
})
