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
