test_that("spatial coordinate resolver prefers x/y before col/row", {
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 3), sparse = TRUE)
  rownames(counts) <- c("Gene1", "Gene2")
  colnames(counts) <- paste0("Spot", 1:3)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(10, 20, 30)
  srt$y <- c(3, 4, 5)
  srt$col <- c(1, 2, 3)
  srt$row <- c(1, 1, 1)
  srt$group <- c("A", "B", "A")

  p <- SpatialSpotPlot(srt, group.by = "group", overlay_image = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, c(10, 20, 30))
  expect_equal(p$data$y, c(3, 4, 5))
})

test_that("spatial plots handle single-value and empty pie states", {
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 3), sparse = TRUE)
  rownames(counts) <- c("Gene1", "Gene2")
  colnames(counts) <- paste0("Spot", 1:3)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(1, 2, 3)
  srt$y <- c(1, 1, 1)
  values <- stats::setNames(rep(0, 3), colnames(srt))

  p_single <- SpatialSpotPlot(
    srt,
    values = values,
    overlay_image = FALSE,
    theme_use = NULL
  )
  expect_s3_class(p_single, "ggplot")

  zero_pies <- data.frame(A = rep(0, 3), B = rep(0, 3), row.names = colnames(srt))
  p_empty <- SpatialSpotPlot(
    srt,
    values = zero_pies,
    plot_type = "pie",
    overlay_image = FALSE,
    theme_use = NULL
  )
  expect_s3_class(p_empty, "ggplot")
})
