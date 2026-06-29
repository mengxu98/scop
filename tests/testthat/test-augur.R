test_that("RunAugur writes per-cell AUC metadata from factor cell types", {
  counts <- matrix(
    c(
      5, 0, 4, 0,
      0, 3, 0, 2,
      2, 1, 3, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:4))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$celltype <- factor(c("alpha", "beta", "alpha", "beta"))
  srt$tech <- factor(c("A", "A", "B", "B"))

  testthat::local_mocked_bindings(
    augur_cpp = function(...) {
      list(
        AUC = data.frame(
          cell_type = c("alpha", "beta"),
          auc = c(0.8, 0.6)
        )
      )
    }
  )

  out <- RunAugur(
    srt,
    celltype.by = "celltype",
    label.by = "tech",
    n_subsamples = 1,
    subsample_size = 2,
    min_cells = 1,
    verbose = FALSE
  )

  expect_true("augur_auc" %in% colnames(out@meta.data))
  expect_true("augur_rank" %in% colnames(out@meta.data))
  expect_equal(unname(out$augur_auc), c(0.8, 0.6, 0.8, 0.6))
  expect_equal(unname(out$augur_rank), c(1, 2, 1, 2))
})
