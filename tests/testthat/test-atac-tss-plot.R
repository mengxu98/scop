make_atac_tss_plot_fixture <- function() {
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

test_that("TSS profiles require genuine positional matrices and match Signac", {
  object <- make_atac_tss_plot_fixture()
  object$TSS.enrichment <- 5
  expect_error(ATACTSSPlot(object), "fast = FALSE")
  matrix <- matrix(seq_len(20) / 5, nrow = 4,
    dimnames = list(colnames(object), as.character(-2:2)))
  object[["peaks"]]@positionEnrichment[["TSS"]] <- rbind(matrix,
    expected = rep(1, 5), motif = c(0, 0, 1, 0, 0))
  before <- serialize(object, NULL)
  actual <- suppressWarnings(ATACTSSPlot(object, group.by = "group"))
  native <- suppressWarnings(Signac::TSSPlot(object, group.by = "group"))
  expect_equal(actual$data, native$data)
  subset_plot <- suppressWarnings(ATACTSSPlot(object, group.by = "group", idents = "one"))
  expect_identical(unique(subset_plot$data$group), "one")
  single <- suppressWarnings(ATACTSSPlot(object, group.by = "group", cells = "A"))
  expect_equal(single$data$norm.value, as.numeric(matrix["A", ]))
  expect_identical(serialize(object, NULL), before)
  expect_error(ATACTSSPlot(object, group.by = c("group", "group")), "one existing")
  expect_error(ATACTSSPlot(object, idents = "missing"), "existing non-missing")
  object[["peaks"]]@positionEnrichment[["TSS"]][1, 1] <- Inf
  expect_error(ATACTSSPlot(object), "non-finite")
})


test_that("TSS profiles use selected-assay cells in partial-assay objects", {
  object <- make_atac_tss_plot_fixture()
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
  actual <- suppressWarnings(ATACTSSPlot(partial, assay = "peaks", group.by = "group"))
  native <- suppressWarnings(ATACTSSPlot(object, group.by = "group"))
  actual_data <- actual$data[order(actual$data$group, actual$data$position), ]
  native_data <- native$data[order(native$data$group, native$data$position), ]
  rownames(actual_data) <- rownames(native_data) <- NULL
  expect_equal(actual_data, native_data)
})


test_that("tss plot input errors are actionable", {
  expect_error(ATACTSSPlot(list()), "Seurat")
  object <- make_atac_tss_plot_fixture()
  expect_error(ATACTSSPlot(object, assay = "missing"), "ChromatinAssay")
  expect_error(ATACTSSPlot(object, group.by = "missing"), "metadata column")
  object$group[1] <- NA
  expect_error(ATACTSSPlot(object, group.by = "group"), "non-missing group")
})
