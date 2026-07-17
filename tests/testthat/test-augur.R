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

test_that("Augur sparse subsampling matches the dense variable-gene reference", {
  mat <- methods::as(Matrix::Matrix(
    matrix(
      c(
        0, 0, 0, 0, 0, 0,
        1, 0, 1, 0, 1, 0,
        2, 2, 2, 2, 2, 2,
        0, 4, 0, 4, 0, 4,
        5, 1, 5, 1, 2, 1
      ),
      nrow = 5,
      byrow = TRUE,
      dimnames = list(paste0("Gene", 1:5), paste0("Cell", 1:6))
    ),
    sparse = TRUE
  ), "dgCMatrix")
  selected <- c(1L, 3L, 5L)
  dense <- t(as.matrix(mat[, selected, drop = FALSE]))
  expected <- dense[, apply(dense, 2L, function(x) any(x != x[[1L]])), drop = FALSE]

  actual <- augur_subsample_cpp(mat, selected)
  expect_identical(actual, expected)
})
