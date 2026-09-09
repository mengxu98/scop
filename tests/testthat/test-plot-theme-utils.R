test_that("apply_plot_theme accepts evaluated theme objects without crashing", {
  th <- ggplot2::theme_minimal()
  expect_identical(apply_plot_theme(th, list(base_size = 14)), th)
})

test_that("apply_plot_theme resolves theme names and functions", {
  expect_s3_class(apply_plot_theme("theme_scop"), "theme")
  expect_s3_class(apply_plot_theme("theme_this"), "theme")
  expect_s3_class(apply_plot_theme("theme_spatial"), "theme")
  expect_s3_class(apply_plot_theme(ggplot2::theme_bw), "theme")

  th <- apply_plot_theme("theme_minimal", list(base_size = 13, not_a_theme_arg = 1))
  expect_s3_class(th, "theme")
  expect_equal(th$text$size, 13)
})

test_that("apply_plot_theme handles NULL and unknown names", {
  expect_null(apply_plot_theme(NULL, allow_null = TRUE))
  expect_s3_class(apply_plot_theme(NULL), "theme")
  expect_s3_class(apply_plot_theme("not_a_theme_name"), "theme")
})

test_that("deprecated enrlichmap_nlabel stays as a NULL fallback", {
  for (fun in c("EnrichmentPlot", "GSEAPlot", "GSVAPlot")) {
    fmls <- formals(get(fun, envir = asNamespace("scop")))
    expect_identical(fmls$enrichmap_nlabel, 4, info = fun)
    expect_null(fmls$enrlichmap_nlabel, info = fun)
  }
})

test_that("FeatureStatPlot proceeds without prompts in non-interactive sessions", {
  data("pancreas_sub", package = "scop")
  cells <- colnames(pancreas_sub)[1:50]
  stat_by <- head(rownames(pancreas_sub), 55)

  expect_false(interactive())
  p <- suppressWarnings(FeatureStatPlot(
    pancreas_sub,
    stat.by = stat_by,
    group.by = "CellType",
    cells = cells,
    layer = "counts"
  ))
  expect_false(is.null(p))
})

test_that("CoverageTrackPlot preserves factor level order for palette_colors", {
  data("pbmcmultiome_sub", package = "scop")
  called_cols <- NULL
  testthat::local_mocked_bindings(
    CoveragePlot = function(..., cols) {
      called_cols <<- cols
      NULL
    },
    .package = "Signac"
  )
  rev_levels <- rev(as.character(unique(pbmcmultiome_sub$CellType)))
  pbmcmultiome_sub$test_factor <- factor(pbmcmultiome_sub$CellType, levels = rev_levels)
  CoverageTrackPlot(
    pbmcmultiome_sub,
    region = rownames(pbmcmultiome_sub[["peaks"]])[1],
    assay = "peaks",
    group.by = "test_factor",
    verbose = FALSE
  )
  expect_identical(names(called_cols), rev_levels)
})

test_that("spatial_dim_continuous_scale respects upper_quantile and squish", {
  vals <- 1:100
  sc99 <- scop:::spatial_dim_continuous_scale(vals, colors = c("blue", "red"), upper_quantile = 0.99)
  sc100 <- scop:::spatial_dim_continuous_scale(vals, colors = c("blue", "red"), upper_quantile = 1)
  expect_equal(unname(sc99$limits), c(1, 99.01))
  expect_equal(unname(sc100$limits), c(1, 100))
  expect_equal(sc99$oob(c(-10, 50, 150), range = c(0, 100)), c(0, 50, 100))
})
