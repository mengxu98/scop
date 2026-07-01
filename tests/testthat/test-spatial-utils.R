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

  p_explicit <- SpatialSpotPlot(
    srt,
    group.by = "group",
    coord.cols = c("col", "row"),
    overlay_image = FALSE
  )
  expect_equal(p_explicit$data$x, c(1, 2, 3))
  expect_equal(p_explicit$data$y, c(1, 1, 1))
})

test_that("spatial coordinate resolver falls back to col/row and validates explicit columns", {
  counts <- Matrix::Matrix(matrix(1, nrow = 2, ncol = 3), sparse = TRUE)
  rownames(counts) <- c("Gene1", "Gene2")
  colnames(counts) <- paste0("Spot", 1:3)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$col <- c(5, 6, 7)
  srt$row <- c(2, 2, 3)
  srt$group <- c("A", "B", "A")

  p <- SpatialSpotPlot(srt, group.by = "group", overlay_image = FALSE)
  expect_equal(p$data$x, c(5, 6, 7))
  expect_equal(p$data$y, c(2, 2, 3))

  expect_error(
    SpatialSpotPlot(srt, group.by = "group", coord.cols = "col", overlay_image = FALSE),
    "at least two"
  )
  expect_error(
    SpatialSpotPlot(srt, group.by = "group", coord.cols = c("missing_x", "row"), overlay_image = FALSE),
    "Missing metadata"
  )
})

test_that("spatial palette resolver keeps shared defaults and explicit palettes", {
  expect_identical(scop_spatial_palette(NULL, c("A", "B"), type = "discrete"), "Chinese")
  expect_identical(scop_spatial_palette(NULL, c(1, 2), type = "continuous"), "Spectral")
  expect_identical(scop_spatial_palette(NULL, c(-1, 1), type = "diverging"), "RdBu")
  expect_identical(scop_spatial_palette("Set1", c("A", "B"), type = "discrete"), "Set1")
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
