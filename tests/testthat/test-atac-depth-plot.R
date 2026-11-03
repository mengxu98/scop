make_atac_depth_plot_fixture <- function() {
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  counts <- Matrix::Matrix(matrix(c(1, 2, 2, 3, 3, 4, 4, 5), 2L,
    dimnames = list(c("chr1-1-100", "chr1-101-200"), c("A", "B", "C", "D"))), sparse = TRUE)
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 101), width = 100))
  object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
  object$group <- c("one", "one", "two", "two")
  object$verified.depth <- c(10, 20, 30, 40)
  embedding <- cbind(LSI_1 = 1:4, LSI_2 = 4:1, LSI_3 = rep(1, 4))
  rownames(embedding) <- colnames(object)
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(embedding, assay = "peaks", key = "LSI_")
  object
}

test_that("depth correlations are aligned, signed, and support one component", {
  object <- make_atac_depth_plot_fixture()
  before <- serialize(object, NULL)
  result <- ATACDepthCorPlot(object, depth.column = "verified.depth", return.data = TRUE)
  expect_equal(result$correlation, c(1, -1, NA))
  expect_equal(result$n.cells, rep(4L, 3))
  one <- ATACDepthCorPlot(object, n = 1, return.data = TRUE)
  expect_equal(nrow(one), 1L)
  expect_equal(one$correlation, 1)
  expect_equal(ATACDepthCorPlot(object, method = "spearman", return.data = TRUE)$correlation, c(1, -1, NA))
  expect_s3_class(ATACDepthCorPlot(object, n = 1), "ggplot")
  expect_equal(ATACDepthCorPlot(object, n = NULL, return.data = TRUE),
    ATACDepthCorPlot(object, n = 999, return.data = TRUE))
  expect_equal(ATACDepthCorPlot(object, cells = c("D", "A", "C"),
    depth.column = "verified.depth", return.data = TRUE)$correlation, c(1, -1, NA))
  expect_identical(serialize(object, NULL), before)
  object$nCount_peaks <- c(20, 1, 8, 3)
  expect_equal(ATACDepthCorPlot(object, n = 1, return.data = TRUE)$correlation, 1)
  expect_equal(ATACDepthCorPlot(object, n = 1, depth.column = "nCount_peaks", return.data = TRUE)$correlation,
    stats::cor(1:4, c(20, 1, 8, 3)))
})


test_that("depth diagnostics do not silently discard cells or malformed inputs", {
  object <- make_atac_depth_plot_fixture()
  expect_error(ATACDepthCorPlot(object, n = 0), "positive integer")
  expect_error(ATACDepthCorPlot(object, n = 1.5), "positive integer")
  expect_error(ATACDepthCorPlot(object, reduction = "missing"), "existing dimensional")
  expect_error(ATACDepthCorPlot(object, depth.column = "missing"), "metadata column")
  expect_error(ATACDepthCorPlot(object, cells = c("A", "A", "C")), "unique cell")
  expect_error(ATACDepthCorPlot(object, cells = c("A", "B")), "at least three")
  object$verified.depth <- 1
  expect_error(ATACDepthCorPlot(object, depth.column = "verified.depth"), "must vary")
  object$verified.depth <- c(1, 2, NA, 4)
  expect_error(ATACDepthCorPlot(object, depth.column = "verified.depth"), "finite")
  object$verified.depth <- c(1, 2, -1, 4)
  expect_error(ATACDepthCorPlot(object, depth.column = "verified.depth"), "non-negative")
  object$verified.depth <- letters[1:4]
  expect_error(ATACDepthCorPlot(object, depth.column = "verified.depth"), "finite")
  suppressWarnings(object[["lsi"]]@cell.embeddings <- object[["lsi"]]@cell.embeddings[1:3, , drop = FALSE])
  expect_error(ATACDepthCorPlot(object), "every selected assay cell")
})


test_that("Depth correlations use selected-assay cells in partial-assay objects", {
  object <- make_atac_depth_plot_fixture()
  matrix <- matrix(seq_len(20), nrow = 4,
    dimnames = list(colnames(object), as.character(-2:2)))
  object[["peaks"]]@positionEnrichment[["TSS"]] <- matrix
  rna <- Matrix::Matrix(matrix(1, 2, 5,
    dimnames = list(c("g1", "g2"), c("D", "RNA_only", "A", "B", "C"))), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  partial[["peaks"]] <- object[["peaks"]]
  partial[["lsi"]] <- object[["lsi"]]
  partial$group <- c("two", "rna", "one", "one", "two")
  partial$nCount_peaks <- c(9, NA, 3, 5, 7)
  expect_equal(ATACDepthCorPlot(partial, assay = "peaks", return.data = TRUE)$correlation,
    c(1, -1, NA))
})

test_that("depth plot input errors are actionable", {
  expect_error(ATACDepthCorPlot(list()), "Seurat")
  object <- make_atac_depth_plot_fixture()
  expect_error(ATACDepthCorPlot(object, assay = "missing"), "ChromatinAssay")
})
