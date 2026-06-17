test_that("Milo C++ neighborhood kernels count undirected one-hop samples", {
  knn_idx <- matrix(
    as.integer(c(
      2, 3,
      1, 3,
      1, 2,
      5, 3,
      4, 3
    )),
    nrow = 5,
    byrow = TRUE
  )
  out <- milo_nhood_counts_cpp(
    knn_idx = knn_idx,
    sampled_vertices = as.integer(c(1, 4)),
    sample_id = as.integer(c(1, 1, 2, 2, 1)),
    n_samples = 2L,
    k_dist = rep(1, 5)
  )

  expect_equal(out$counts, matrix(as.integer(c(2, 1, 1, 2)), nrow = 2, byrow = TRUE))
  expect_equal(out$nhood_size, as.integer(c(3, 3)))
  expect_equal(out$sampled_vertices, as.integer(c(1, 4)))
})

test_that("Milo weighted FDR preserves NA p-values and monotone adjustment", {
  out <- milo_weighted_fdr_cpp(
    pvalues = c(0.01, 0.20, NA, 0.03),
    weights = c(1, 2, 1, 1)
  )

  expect_equal(is.na(out), c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(out[c(1, 2, 4)], c(0.04, 0.20, 0.06), tolerance = 1e-12)
})

test_that("RunMilo exposes r and cpp backend choices", {
  expect_equal(eval(formals(RunMilo)$backend), c("r", "cpp"))
})
