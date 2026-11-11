test_that("RunCCA accepts sparse matrices and matches Seurat's dense result", {
  skip_if_not_installed("Seurat")

  set.seed(92)
  object1 <- Matrix::rsparsematrix(40, 18, density = 0.3)
  object2 <- Matrix::rsparsematrix(40, 16, density = 0.3)
  rownames(object1) <- rownames(object2) <- paste0("g", seq_len(40))
  colnames(object1) <- paste0("a", seq_len(ncol(object1)))
  colnames(object2) <- paste0("b", seq_len(ncol(object2)))

  actual <- RunCCA.default(
    object1,
    object2,
    num.cc = 6,
    seed.use = 92,
    verbose = FALSE
  )
  expected <- get("RunCCA.default", envir = asNamespace("Seurat"))(
    as.matrix(object1),
    as.matrix(object2),
    num.cc = 6,
    seed.use = 92,
    verbose = FALSE
  )

  expect_equal(actual$d, expected$d, tolerance = 1e-10)
  qa <- qr.Q(qr(actual$ccv))
  qb <- qr.Q(qr(expected$ccv))
  expect_equal(
    svd(crossprod(qa, qb), nu = 0, nv = 0)$d,
    rep(1, 6),
    tolerance = 1e-8
  )
})
