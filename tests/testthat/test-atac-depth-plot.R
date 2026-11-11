test_that("depth correlations are aligned, signed, and support one component", {
  object <- make_atac_plot_test_object()
  before <- serialize(object, NULL)
  result <- ATACDepthCorPlot(object, depth.column = "verified.depth", return.data = TRUE)
  expect_equal(result$correlation, c(1, -1, NA))
  expect_equal(result$n.cells, rep(4L, 3))
  one <- ATACDepthCorPlot(object, n = 1, return.data = TRUE)
  expect_equal(nrow(one), 1L)
  expect_equal(one$correlation, 1)
  expect_equal(ATACDepthCorPlot(object, method = "spearman", return.data = TRUE)$correlation, c(1, -1, NA))
  expect_s3_class(ATACDepthCorPlot(object, n = 1), "ggplot")
  expect_equal(
    ATACDepthCorPlot(object, n = NULL, return.data = TRUE),
    ATACDepthCorPlot(object, n = 999, return.data = TRUE)
  )
  expect_equal(ATACDepthCorPlot(object,
    cells = c("D", "A", "C"),
    depth.column = "verified.depth", return.data = TRUE
  )$correlation, c(1, -1, NA))
  expect_identical(serialize(object, NULL), before)
  object$nCount_peaks <- c(20, 1, 8, 3)
  expect_equal(ATACDepthCorPlot(object, n = 1, return.data = TRUE)$correlation, 1)
  expect_equal(
    ATACDepthCorPlot(object, n = 1, depth.column = "nCount_peaks", return.data = TRUE)$correlation,
    stats::cor(1:4, c(20, 1, 8, 3))
  )
})


test_that("depth diagnostics do not silently discard cells or malformed inputs", {
  object <- make_atac_plot_test_object()
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
  object <- make_atac_plot_test_object()
  matrix <- matrix(seq_len(20),
    nrow = 4,
    dimnames = list(colnames(object), as.character(-2:2))
  )
  object[["peaks"]]@positionEnrichment[["TSS"]] <- matrix
  rna <- Matrix::Matrix(matrix(1, 2, 5,
    dimnames = list(c("g1", "g2"), c("D", "RNA_only", "A", "B", "C"))
  ), sparse = TRUE)
  partial <- SeuratObject::CreateSeuratObject(rna, assay = "RNA")
  partial[["peaks"]] <- object[["peaks"]]
  partial[["lsi"]] <- object[["lsi"]]
  partial$group <- c("two", "rna", "one", "one", "two")
  partial$nCount_peaks <- c(9, NA, 3, 5, 7)
  expect_equal(
    ATACDepthCorPlot(partial, assay = "peaks", return.data = TRUE)$correlation,
    c(1, -1, NA)
  )
})

test_that("depth plot input errors are actionable", {
  expect_error(ATACDepthCorPlot(list()), "Seurat")
  object <- make_atac_plot_test_object()
  expect_error(ATACDepthCorPlot(object, assay = "missing"), "ChromatinAssay")
})
