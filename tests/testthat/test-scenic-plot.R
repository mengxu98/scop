make_scenic_plot_mock <- function(seed = 1) {
  set.seed(seed)
  counts <- Matrix::Matrix(matrix(
    stats::rpois(8 * 12, lambda = 5),
    nrow = 8,
    ncol = 12,
    dimnames = list(
      paste0("Gene", seq_len(8)),
      paste0("Cell", seq_len(12))
    )
  ), sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$CellType <- rep(c("A", "B", "C"), each = 4)

  auc <- matrix(
    stats::runif(6 * 12),
    nrow = 6,
    ncol = 12,
    dimnames = list(
      paste0(c("Jun", "Atf3", "Fos", "Klf4", "Sox9", "Birc3"), "(+)"),
      colnames(srt)
    )
  )
  srt@tools$SCENIC <- list(scores_cells_by_regulon = t(auc))
  list(srt = srt, auc = auc)
}

test_that("SCENICPlot rss heatmap returns drawable plot object", {
  dat <- make_scenic_plot_mock()
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "rss_heatmap",
    features = rownames(dat$auc)[1:4],
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$plots[[1]], "ggplot")
  expect_false("matrix_list" %in% names(out$plot))
  expect_s3_class(out$heatmap$plot, "ggplot")
  expect_true("matrix_list" %in% names(out$heatmap))
})

test_that("SCENICPlot activity heatmap returns drawable plot object", {
  dat <- make_scenic_plot_mock(seed = 2)
  out <- SCENICPlot(
    dat$srt,
    group.by = "CellType",
    plot_type = "activity_heatmap",
    features = rownames(dat$auc)[1:4],
    verbose = FALSE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$plots[[1]], "ggplot")
  expect_false("matrix_list" %in% names(out$plot))
  expect_s3_class(out$heatmap$plot, "ggplot")
  expect_true("matrix_list" %in% names(out$heatmap))
})
