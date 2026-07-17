# Tests for Palantir C++ backend
#
# Structure:
#   1. Helper: mock embedding/KNN data
#   2. palantir_compute_kernel_cpp — adaptive anisotropic kernel
#   3. palantir_normalize_kernel_cpp — diffusion map normalization
#   4. palantir_multiscale_space_cpp — eigenvalue-scaled components
#   5. palantir_numpy_random_sample_cpp — NumPy-compatible RNG
#   6. palantir_maxmin_waypoints_cpp — max-min waypoint sampling
#   7. palantir_pseudotime_cpp — pseudotime via Dijkstra + refinement
#   8. palantir_markov_chain_cpp — directed Markov chain
#   9. palantir_absorption_cpp — absorption probabilities (I-Q solve)

# ---------------------------------------------------------------------------
# 1. Helpers
# ---------------------------------------------------------------------------

make_palantir_mock <- function(
  n_cells = 30,
  n_dims = 5,
  n_neighbors = 8,
  seed = 42
) {
  set.seed(seed)
  embedding <- matrix(
    stats::rnorm(n_cells * n_dims, sd = 2),
    nrow = n_cells,
    ncol = n_dims
  )
  storage.mode(embedding) <- "double"

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  knn_dist <- matrix(NA_real_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    dists <- sqrt(colSums((t(embedding) - embedding[i, ])^2))
    dists[i] <- Inf
    ord <- order(dists)[seq_len(n_neighbors)]
    knn_idx[i, ] <- ord
    knn_dist[i, ] <- dists[ord]
  }

  list(
    embedding = embedding,
    knn_idx = knn_idx,
    knn_dist = knn_dist,
    n_cells = n_cells,
    n_dims = n_dims,
    n_neighbors = n_neighbors
  )
}

# ---------------------------------------------------------------------------
# 2. palantir_compute_kernel_cpp
# ---------------------------------------------------------------------------

test_that("palantir_compute_kernel_cpp returns correct structure", {
  dat <- make_palantir_mock()
  out <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.0
  )

  expect_type(out, "list")
  expect_true(all(c("i", "j", "x", "n", "adaptive_std") %in% names(out)))
  expect_equal(out$n, dat$n_cells)
  expect_length(out$adaptive_std, dat$n_cells)
  expect_true(all(out$adaptive_std > 0))
  # All values should be finite and positive
  expect_true(all(is.finite(out$x)))
  expect_true(all(out$x > 0))
  # i, j should be within [1, n]
  expect_true(all(out$i >= 1 & out$i <= dat$n_cells))
  expect_true(all(out$j >= 1 & out$j <= dat$n_cells))
  # Kernel should have some non-self edges
  non_self <- out$i != out$j
  expect_gt(sum(non_self), 0)
})

test_that("palantir_compute_kernel_cpp with alpha > 0 produces different kernel", {
  dat <- make_palantir_mock()
  out0 <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.0
  )
  out1 <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.5
  )
  # Values should differ when alpha > 0
  expect_false(isTRUE(all.equal(out0$x, out1$x)))
})

test_that("palantir_compute_kernel_cpp handles all-zero embedding", {
  n_cells <- 10
  n_dims <- 3
  embedding <- matrix(0, nrow = n_cells, ncol = n_dims)
  storage.mode(embedding) <- "double"
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  knn_dist <- matrix(0, nrow = n_cells, ncol = 3)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- head(candidates, 3)
    knn_dist[i, ] <- rep(0, 3)
  }

  out <- palantir_compute_kernel_cpp(
    data = embedding,
    knn_idx = knn_idx,
    knn_dist = knn_dist,
    knn = n_cells - 1,
    alpha = 0.0
  )
  # adaptive_std should be a small positive number (clamped at 1e-10)
  expect_true(all(out$adaptive_std >= 1e-10))
  # Values should be exp(0) = 1, then duplicated by W + W.T
  expect_true(all(out$x > 0))
  expect_true(all(is.finite(out$x)))
})

test_that("palantir_compute_kernel_cpp is deterministic", {
  dat <- make_palantir_mock(seed = 42)
  out1 <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.0
  )
  out2 <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.0
  )
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 3. palantir_normalize_kernel_cpp
# ---------------------------------------------------------------------------

test_that("palantir_normalize_kernel_cpp returns correct structure", {
  dat <- make_palantir_mock()
  kernel <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.0
  )
  norm <- palantir_normalize_kernel_cpp(
    kernel_i = as.integer(kernel$i),
    kernel_j = as.integer(kernel$j),
    kernel_x = as.numeric(kernel$x),
    n = as.integer(kernel$n)
  )

  expect_type(norm, "list")
  expect_true(all(c("T_i", "T_j", "T_x", "L_i", "L_j", "L_x", "n", "row_sum") %in% names(norm)))
  expect_equal(norm$n, dat$n_cells)
  expect_length(norm$row_sum, dat$n_cells)
  expect_true(all(norm$row_sum > 0))
  expect_true(all(is.finite(norm$T_x)))
  expect_true(all(norm$T_x > 0))
})

test_that("palantir_normalize_kernel_cpp produces row-stochastic T matrix", {
  dat <- make_palantir_mock(n_cells = 10)
  kernel <- palantir_compute_kernel_cpp(
    data = dat$embedding,
    knn_idx = dat$knn_idx,
    knn_dist = dat$knn_dist,
    knn = dat$n_neighbors,
    alpha = 0.0
  )
  norm <- palantir_normalize_kernel_cpp(
    kernel_i = as.integer(kernel$i),
    kernel_j = as.integer(kernel$j),
    kernel_x = as.numeric(kernel$x),
    n = as.integer(kernel$n)
  )

  # Each row of T should sum to approximately 1
  T_mat <- Matrix::sparseMatrix(
    i = norm$T_i,
    j = norm$T_j,
    x = norm$T_x,
    dims = c(norm$n, norm$n)
  )
  row_sums <- Matrix::rowSums(T_mat)
  expect_equal(row_sums, rep(1, norm$n), tolerance = 1e-10)
})

test_that("palantir_normalize_kernel_cpp handles uniform kernel", {
  n <- 5
  kernel_i <- as.integer(c(1,1,2,2,3,3,4,4,5,5))
  kernel_j <- as.integer(c(2,3,1,4,1,5,2,5,3,4))
  kernel_x <- rep(1.0, 10)

  norm <- palantir_normalize_kernel_cpp(
    kernel_i = kernel_i,
    kernel_j = kernel_j,
    kernel_x = kernel_x,
    n = n
  )

  expect_equal(norm$n, n)
  expect_true(all(is.finite(norm$T_x)))
})

# ---------------------------------------------------------------------------
# 4. palantir_multiscale_space_cpp
# ---------------------------------------------------------------------------

test_that("palantir_multiscale_space_cpp returns correct dimensions", {
  n <- 20
  n_eigs <- 8
  set.seed(1)
  eigvecs <- matrix(rnorm(n * n_eigs), nrow = n, ncol = n_eigs)
  # Normalize columns
  eigvecs <- apply(eigvecs, 2, function(x) x / sqrt(sum(x^2)))
  eigvals <- seq(0.99, 0.5, length.out = n_eigs)  # decreasing

  ms <- palantir_multiscale_space_cpp(eigvecs, eigvals)

  expect_true(is.matrix(ms))
  expect_equal(nrow(ms), n)
  # Number of components = n_use - 1 (first eigenvector skipped)
  expect_gte(ncol(ms), 2)   # at least 3 eigs used → 2 components minimum
  expect_lte(ncol(ms), n_eigs - 1)
  expect_true(all(is.finite(ms)))
})

test_that("palantir_multiscale_space_cpp scaling formula", {
  n <- 10
  n_eigs <- 5
  set.seed(2)
  eigvecs <- matrix(rnorm(n * n_eigs), nrow = n, ncol = n_eigs)
  eigvecs <- apply(eigvecs, 2, function(x) x / sqrt(sum(x^2)))
  eigvals <- c(0.95, 0.8, 0.6, 0.4, 0.2)

  ms <- palantir_multiscale_space_cpp(eigvecs, eigvals)

  # Manual check on first component (j=1 in 0-based indexing, eigenvalue 0.8)
  manual <- eigvecs[, 2] * eigvals[2] / (1.0 - eigvals[2])
  expect_equal(ms[, 1], manual, tolerance = 1e-10)
})

test_that("palantir_multiscale_space_cpp handles edge case: all equal eigenvalues", {
  n <- 15
  eigvecs <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  eigvecs <- apply(eigvecs, 2, function(x) x / sqrt(sum(x^2)))
  eigvals <- rep(0.9, 5)

  ms <- palantir_multiscale_space_cpp(eigvecs, eigvals)
  expect_true(is.matrix(ms))
  expect_true(all(is.finite(ms)))
})

# ---------------------------------------------------------------------------
# 5. palantir_numpy_random_sample_cpp
# ---------------------------------------------------------------------------

test_that("palantir_numpy_random_sample_cpp returns in [0, 1]", {
  out <- palantir_numpy_random_sample_cpp(100, seed = 42)
  expect_length(out, 100)
  expect_true(all(out >= 0))
  expect_true(all(out < 1))
  expect_type(out, "double")
})

test_that("palantir_numpy_random_sample_cpp is reproducible", {
  out1 <- palantir_numpy_random_sample_cpp(50, seed = 42)
  out2 <- palantir_numpy_random_sample_cpp(50, seed = 42)
  expect_equal(out1, out2)
})

test_that("palantir_numpy_random_sample_cpp different seeds differ", {
  out1 <- palantir_numpy_random_sample_cpp(100, seed = 1)
  out2 <- palantir_numpy_random_sample_cpp(100, seed = 2)
  expect_false(isTRUE(all.equal(out1, out2)))
})

# ---------------------------------------------------------------------------
# 6. palantir_maxmin_waypoints_cpp
# ---------------------------------------------------------------------------

test_that("palantir_maxmin_waypoints_cpp returns valid indices", {
  n <- 50
  n_dims <- 4
  set.seed(1)
  ms_data <- matrix(rnorm(n * n_dims), nrow = n, ncol = n_dims)

  wp <- palantir_maxmin_waypoints_cpp(ms_data, num_waypoints = 10, seed = 42)
  expect_true(all(wp >= 1))
  expect_true(all(wp <= n))
  expect_true(all(wp == as.integer(wp)))
})

test_that("palantir_maxmin_waypoints_cpp respects minimum waypoints", {
  n <- 20
  ms_data <- matrix(rnorm(n * 3), nrow = n, ncol = 3)

  # Request fewer than n_cols → should be bumped up
  wp <- palantir_maxmin_waypoints_cpp(ms_data, num_waypoints = 2, seed = 42)
  expect_gte(length(wp), 3)   # at least max(3, n_cols)
})

test_that("palantir_maxmin_waypoints_cpp is deterministic", {
  n <- 30
  ms_data <- matrix(rnorm(n * 4), nrow = n, ncol = 4)

  wp1 <- palantir_maxmin_waypoints_cpp(ms_data, num_waypoints = 8, seed = 42)
  wp2 <- palantir_maxmin_waypoints_cpp(ms_data, num_waypoints = 8, seed = 42)
  expect_equal(wp1, wp2)
})

# ---------------------------------------------------------------------------
# 7. palantir_pseudotime_cpp
# ---------------------------------------------------------------------------

test_that("palantir_pseudotime_cpp returns valid structure", {
  n <- 20
  d <- 3
  set.seed(1)
  ms_data <- matrix(rnorm(n * d), nrow = n, ncol = d)

  out <- palantir_pseudotime_cpp(
    ms_data = ms_data,
    start_cell = 0,          # 0-based
    waypoints = as.integer(c(1, 5, 10, 15, 20)),
    knn = 5,
    max_iterations = 10,
    n_jobs = 1
  )

  expect_type(out, "list")
  expect_true(all(c("pseudotime", "W") %in% names(out)))
  expect_length(out$pseudotime, n)
  expect_true(all(is.finite(out$pseudotime)))
  # Pseudotime should be normalized to [0, 1]
  expect_gte(min(out$pseudotime), 0)
  expect_lte(max(out$pseudotime), 1)
})

test_that("palantir_pseudotime_cpp start_cell has minimum pseudotime", {
  n <- 15
  d <- 3
  set.seed(2)
  ms_data <- matrix(rnorm(n * d), nrow = n, ncol = d)

  start_cell <- 0  # first cell
  out <- palantir_pseudotime_cpp(
    ms_data = ms_data,
    start_cell = start_cell,
    waypoints = as.integer(seq(1, n, by = 3)),
    knn = 4,
    max_iterations = 10,
    n_jobs = 1
  )

  # Start cell should have pseudotime near 0
  expect_equal(out$pseudotime[start_cell + 1], 0, tolerance = 0.05)
})

test_that("palantir_pseudotime_cpp rejects invalid start_cell", {
  ms_data <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(
    palantir_pseudotime_cpp(
      ms_data = ms_data,
      start_cell = -1,
      waypoints = as.integer(c(1, 5)),
      knn = 3,
      max_iterations = 5
    ),
    "start_cell out of range"
  )
  expect_error(
    palantir_pseudotime_cpp(
      ms_data = ms_data,
      start_cell = 100,
      waypoints = as.integer(c(1, 5)),
      knn = 3,
      max_iterations = 5
    ),
    "start_cell out of range"
  )
})

test_that("palantir_pseudotime_cpp is deterministic", {
  n <- 12
  d <- 2
  set.seed(3)
  ms_data <- matrix(rnorm(n * d), nrow = n, ncol = d)

  args <- list(
    ms_data = ms_data,
    start_cell = 0L,
    waypoints = as.integer(c(1, 4, 8, 12)),
    knn = 4,
    max_iterations = 10,
    n_jobs = 1
  )
  out1 <- do.call(palantir_pseudotime_cpp, args)
  out2 <- do.call(palantir_pseudotime_cpp, args)
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 8. palantir_markov_chain_cpp
# ---------------------------------------------------------------------------

test_that("palantir_markov_chain_cpp returns correct structure", {
  n <- 15
  d <- 3
  set.seed(1)
  wp_data <- matrix(rnorm(n * d), nrow = n, ncol = d)
  pseudotime <- seq(0, 1, length.out = n)

  out <- palantir_markov_chain_cpp(
    wp_data = wp_data,
    knn = 5,
    pseudotime = pseudotime
  )

  expect_type(out, "list")
  expect_true(all(c("T_i", "T_j", "T_x", "n") %in% names(out)))
  expect_equal(out$n, n)
  expect_true(all(is.finite(out$T_x)))
  # T should be row-stochastic or have self-loops
  T_mat <- Matrix::sparseMatrix(
    i = out$T_i,
    j = out$T_j,
    x = out$T_x,
    dims = c(n, n)
  )
  row_sums <- Matrix::rowSums(T_mat)
  expect_equal(row_sums, rep(1, n), tolerance = 1e-10)
})

test_that("palantir_markov_chain_cpp handles single cell", {
  wp_data <- matrix(rnorm(3), nrow = 1, ncol = 3)
  pseudotime <- 0.5

  out <- palantir_markov_chain_cpp(
    wp_data = wp_data,
    knn = 3,
    pseudotime = pseudotime
  )

  expect_equal(out$n, 1)
  # Single cell should get a self-loop since no valid neighbors
  expect_true(all(out$T_i == 1))
  expect_true(all(out$T_j == 1))
  expect_equal(out$T_x, 1.0)
})

# ---------------------------------------------------------------------------
# 9. palantir_absorption_cpp
# ---------------------------------------------------------------------------

test_that("palantir_absorption_cpp returns correct dimensions", {
  n <- 10
  # Build a simple Markov chain
  T_i <- as.integer(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10))
  T_j <- as.integer(c(2,3,1,4,1,5,4,5,3,6,5,7,6,8,7,9,8,10,9,10))
  T_x <- rep(0.5, 20)  # uniform transitions

  bp <- palantir_absorption_cpp(
    T_i = T_i,
    T_j = T_j,
    T_x = T_x,
    n = n,
    terminal_state_indices = as.integer(c(1, 10))
  )

  expect_true(is.matrix(bp))
  expect_equal(nrow(bp), n)
  expect_equal(ncol(bp), 2)
  # Terminal states should have probability 1 to themselves
  expect_equal(bp[1, 1], 1.0)
  expect_equal(bp[10, 2], 1.0)
  # All probabilities in [0, 1]
  expect_true(all(bp >= 0 & bp <= 1, na.rm = TRUE))
})

test_that("palantir_absorption_cpp row-sums are 1", {
  n <- 8
  T_i <- as.integer(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
  T_j <- as.integer(c(2,3,1,4,1,5,4,5,3,6,5,7,6,8,7,8))
  T_x <- rep(0.5, 16)

  bp <- palantir_absorption_cpp(
    T_i = T_i,
    T_j = T_j,
    T_x = T_x,
    n = n,
    terminal_state_indices = as.integer(c(1, 8))
  )

  rs <- rowSums(bp)
  expect_equal(rs, rep(1, n), tolerance = 1e-10)
})

test_that("palantir_absorption_cpp all cells terminal", {
  n <- 5
  T_i <- as.integer(c(1,2,3,4,5))
  T_j <- as.integer(c(2,3,4,5,1))
  T_x <- rep(1.0, 5)

  bp <- palantir_absorption_cpp(
    T_i = T_i,
    T_j = T_j,
    T_x = T_x,
    n = n,
    terminal_state_indices = as.integer(1:5)
  )

  # Identity matrix
  expect_equal(diag(bp), rep(1, 5))
})

test_that("palantir_row_entropy_cpp matches Palantir branch post-processing", {
  probabilities <- matrix(
    c(
      1, 0, 0,
      0.5, 0.5, 0,
      -0.2, 0.8, 0.4,
      NA_real_, 0.7, 0.3
    ),
    nrow = 4,
    byrow = TRUE
  )
  legacy <- apply(probabilities, 1, function(p) {
    p <- p[p > 0]
    if (length(p) < 2L) 0 else -sum(p * log(p))
  })

  expect_equal(
    scop:::palantir_row_entropy_cpp(probabilities),
    legacy,
    tolerance = 1e-12
  )
})
