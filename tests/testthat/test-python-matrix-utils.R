test_that("Python matrix bridge preserves sparse storage and orientation", {
  bridge <- getFromNamespace("python_cells_by_features", "scop")
  sparse <- Matrix::rsparsematrix(120, 80, density = 0.03)
  dense <- matrix(seq_len(35), nrow = 7)

  sparse_out <- bridge(sparse)
  dense_out <- bridge(dense)

  expect_s4_class(sparse_out, "sparseMatrix")
  expect_identical(dim(sparse_out), rev(dim(sparse)))
  expect_equal(as.matrix(sparse_out), t(as.matrix(sparse)))
  expect_identical(dense_out, t(dense))
})

test_that("minimal AnnData bridge preserves names and aligns embeddings", {
  build <- getFromNamespace("python_anndata_from_matrix", "scop")
  x <- Matrix::sparseMatrix(
    i = c(1, 3, 2),
    j = c(1, 1, 2),
    x = c(2, 4, 3),
    dims = c(3, 2),
    dimnames = list(c("g1", "g2", "g3"), c("c1", "c2"))
  )
  embedding <- matrix(
    1:4,
    nrow = 2,
    dimnames = list(c("c2", "c1"), c("PC1", "PC2"))
  )
  module <- new.env(parent = emptyenv())
  module$AnnData <- function(X) {
    out <- new.env(parent = emptyenv())
    out$X <- X
    out$obsm <- list()
    out
  }

  adata <- build(x, list(pca = embedding), anndata_module = module)

  expect_s4_class(adata$X, "CsparseMatrix")
  expect_equal(as.matrix(adata$X), t(as.matrix(x)))
  expect_identical(adata$obs_names, colnames(x))
  expect_identical(adata$var_names, rownames(x))
  expect_equal(adata$obsm$pca, embedding[colnames(x), , drop = FALSE])
})
