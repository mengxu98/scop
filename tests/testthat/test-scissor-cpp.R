test_that("Scissor network wrappers run on real-data-derived inputs", {
  data("pancreas_sub", package = "scop")
  counts <- methods::as(
    GetAssayData5(pancreas_sub, assay = "RNA", layer = "counts"),
    "dgCMatrix"
  )
  x <- t(log1p(as.matrix(counts[seq_len(30), seq_len(120), drop = FALSE])))
  x <- x[, apply(x, 2, stats::sd) > 0, drop = FALSE]
  x <- x[, seq_len(min(20, ncol(x))), drop = FALSE]
  y_gauss <- as.numeric(scale(rowMeans(x[, seq_len(5), drop = FALSE])))
  y_bin <- rep(c(0, 1), length.out = nrow(x))

  omega <- abs(stats::cor(x))
  omega[!is.finite(omega)] <- 0
  diag(omega) <- 0
  omega[upper.tri(omega) | lower.tri(omega)] <- pmax(
    omega[upper.tri(omega) | lower.tri(omega)],
    1e-3
  )
  omega <- (omega + t(omega)) / 2
  omega <- methods::as(Matrix::Matrix(omega, sparse = TRUE), "dgCMatrix")
  lambda <- c(0.2, 0.1, 0.05)

  gauss <- scissor_gaussian_net_fit_cpp(
    x = x,
    y = y_gauss,
    omega = omega,
    alpha = 0.7,
    lambda = lambda,
    nlambda = length(lambda),
    foldid = NULL,
    inzero = FALSE,
    isd = FALSE,
    thresh = 1e-6,
    maxit = 200L,
    threshP = 1e-5
  )
  expect_true(all(is.finite(gauss$Beta)))
  expect_equal(nrow(gauss$fit), length(lambda))

  binomial <- scissor_binomial_net_fit_cpp(
    x = x,
    y = y_bin,
    omega = omega,
    alpha = 0.7,
    lambda = lambda,
    nlambda = length(lambda),
    foldid = NULL,
    inzero = FALSE,
    isd = FALSE,
    thresh = 1e-6,
    maxit = 200L,
    threshP = 1e-5
  )
  expect_true(all(is.finite(binomial$Beta)))
  expect_equal(nrow(binomial$fit), length(lambda))
})
