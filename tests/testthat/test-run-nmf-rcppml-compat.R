test_that("NMF sparse feature variances match base R", {
  set.seed(71)
  x <- Matrix::rsparsematrix(80, 25, density = 0.18)
  expected <- apply(x, 1L, stats::var)

  expect_equal(nmf_feature_variances(x), expected, tolerance = 1e-12)
})

test_that("NMF dense feature variances match base R", {
  mat <- matrix(
    c(1, 2, 3, 4, 2, 3, 5, 7, 3, 4, 8, 11),
    nrow = 3,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:4))
  )
  testthat::local_mocked_bindings(check_r = function(...) TRUE, .package = "scop")

  expect_equal(
    nmf_feature_variances(mat),
    apply(mat, 1, stats::var),
    tolerance = 1e-12
  )
})

test_that("RunNMF works when RcppML has no global thread setter", {
  skip_if_not_installed("RcppML")
  library(Matrix)

  set.seed(11)
  x <- Matrix::rsparsematrix(60, 30, density = 0.2)
  x@x <- abs(x@x)
  rownames(x) <- paste0("gene_", seq_len(nrow(x)))
  colnames(x) <- paste0("cell_", seq_len(ncol(x)))

  fit <- RunNMF.default(
    x,
    nbes = 5,
    maxit = 5,
    verbose = FALSE,
    seed.use = 11
  )

  expect_equal(dim(fit@cell.embeddings), c(ncol(x), 5L))
  expect_equal(dim(fit@feature.loadings), c(nrow(x), 5L))
})
