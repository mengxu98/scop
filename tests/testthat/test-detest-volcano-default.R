test_that("DE plots accept result tables without an object", {
  skip_if_not_installed("geomtextpath")

  res <- data.frame(
    gene = paste0("G", seq_len(8)),
    group1 = rep(c("A", "B"), each = 4),
    avg_log2FC = rep(c(-2, -1, 1, 2), 2),
    p_val = 0.001,
    p_val_adj = 0.008,
    pct.1 = rep(c(0.1, 0.2, 0.8, 0.9), 2),
    pct.2 = 0.5
  )
  calls <- list(
    list(fun = "DEtestPlot", args = list(plot_type = "volcano")),
    list(fun = "DEtestPlot", args = list(plot_type = "manhattan")),
    list(fun = "DEtestPlot", args = list(plot_type = "ring")),
    list(fun = "VolcanoPlot", args = list()),
    list(fun = "DEtestManhattanPlot", args = list()),
    list(fun = "DEtestRingPlot", args = list())
  )
  for (call in calls) {
    args <- c(list(res = res, nlabel = 0), call$args)
    actual <- do.call(call$fun, args)
    expected <- do.call(call$fun, c(list(object = NULL), args))
    expect_s3_class(actual, "ggplot")
    expect_equal(
      ggplot2::ggplot_build(actual)$data,
      ggplot2::ggplot_build(expected)$data
    )
  }
})

test_that("volcano defaults use signed layouts on real pancreas data", {
  skip_if_not_installed("Seurat")

  data("pancreas_sub", package = "scop")
  groups <- c("Ductal", "Endocrine", "Pre-endocrine")
  cells <- unlist(lapply(groups, function(group) {
    head(colnames(pancreas_sub)[as.character(pancreas_sub$CellType) == group], 60L)
  }), use.names = FALSE)
  srt <- Seurat::NormalizeData(
    suppressWarnings(pancreas_sub[, cells]),
    verbose = FALSE
  )
  out <- RunDEtest(
    srt,
    group.by = "CellType",
    only.pos = FALSE,
    fc.threshold = 1,
    min.pct = 0.1,
    verbose = FALSE
  )

  colored_data <- function(plots) {
    do.call(rbind, lapply(plots, function(plot) {
      do.call(rbind, lapply(plot$layers[5:6], function(layer) layer$data))
    }))
  }
  expect_s4_class(out, "Seurat")
  expect_true(any(out@tools$DEtest_CellType$AllMarkers_wilcox$avg_log2FC < 0))

  plot_default <- DEtestPlot(
    out,
    group.by = "CellType",
    plot_type = "volcano",
    nlabel = 0,
    combine = FALSE
  )
  plot_direct <- VolcanoPlot(
    out,
    group.by = "CellType",
    nlabel = 0,
    combine = FALSE
  )
  plot_positive_only <- DEtestPlot(
    out,
    group.by = "CellType",
    plot_type = "volcano",
    DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
    nlabel = 0,
    combine = FALSE
  )

  default_colored <- colored_data(plot_default)
  direct_colored <- colored_data(plot_direct)
  positive_colored <- colored_data(plot_positive_only)
  expect_true(any(default_colored$avg_log2FC_raw < 0))
  expect_true(all(default_colored$y_plot[default_colored$avg_log2FC_raw < 0] < 0))
  expect_true(any(direct_colored$avg_log2FC_raw < 0))
  expect_true(all(direct_colored$y_plot[direct_colored$avg_log2FC_raw < 0] < 0))
  expect_false(any(positive_colored$avg_log2FC_raw < 0))
})

test_that("RunDEtest normalizes a missing data layer for real pancreas data", {
  skip_if_not_installed("Seurat")

  data("pancreas_sub", package = "scop")
  expect_false("data" %in% SeuratObject::Layers(pancreas_sub[["RNA"]]))

  pancreas_sub <- RunDEtest(
    pancreas_sub,
    group.by = "CellType",
    fc.threshold = 1,
    only.pos = FALSE,
    verbose = FALSE
  )
  plots <- DEtestPlot(
    pancreas_sub,
    group.by = "CellType",
    plot_type = "volcano",
    label.size = 2,
    combine = FALSE
  )

  expect_true("data" %in% SeuratObject::Layers(pancreas_sub[["RNA"]]))
  expect_length(plots, length(unique(pancreas_sub$CellType)))
})
