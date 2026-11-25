test_that("Seurat and SpatialExperiment converters preserve coordinates", {
  testthat::skip_if_not_installed("SpatialExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")

  counts <- Matrix::Matrix(matrix(c(1, 0, 2, 3, 4, 5), nrow = 2), sparse = TRUE)
  rownames(counts) <- c("Gene1", "Gene2")
  colnames(counts) <- paste0("Spot", 1:3)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(10, 20, 30)
  srt$y <- c(1, 2, 3)
  srt$sample <- c("S1", "S1", "S2")

  spe <- srt_to_spe(srt, layer = "counts")
  expect_s4_class(spe, "SpatialExperiment")
  expect_equal(unname(as.numeric(SpatialExperiment::spatialCoords(spe)[, 1])), unname(srt$x))
  expect_equal(unname(as.character(SummarizedExperiment::colData(spe)$sample)), unname(srt$sample))

  out <- spe_to_srt(spe, assay = "Spatial")
  expect_s4_class(out, "Seurat")
  expect_equal(unname(out$x), unname(srt$x))
  expect_equal(unname(out$y), unname(srt$y))
})
