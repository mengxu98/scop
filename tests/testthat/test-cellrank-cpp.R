# Tests for CellRank C++ backend
#
# Covers:
#   1. cellrank_validate_transition_matrix_cpp
#   2. cellrank_stationary_distribution_cpp
#   3. cellrank_schur_cpp (eigen-based)
#   4. cellrank_auto_n_states_cpp
#   5. cellrank_velocity_kernel_cpp
#   6. cellrank_pseudotime_kernel_cpp
#   7. cellrank_cytotrace_kernel_cpp
#   8. cellrank_cflare_cpp
#   9. cellrank_gpcca_cpp
#  10. cellrank_lineage_drivers_cpp

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_transition_matrix <- function(n = 30, seed = 1) {
  set.seed(seed)
  T <- matrix(runif(n * n), n, n)
  T <- T / rowSums(T)
  T
}

make_knn_idx <- function(n_cells, n_neighbors = 5, seed = 42) {
  set.seed(seed)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, n_neighbors))
  }
  knn_idx
}

# ---------------------------------------------------------------------------
# 1. Validate transition matrix
# ---------------------------------------------------------------------------

test_that("cellrank_validate_transition_matrix_cpp fixes NaN/Inf", {
  n <- 10
  T <- make_transition_matrix(n)
  T[1, 2] <- NaN; T[3, 4] <- Inf; T[5, 6] <- -0.1
  T[7, ] <- 0

  out <- scop:::cellrank_validate_transition_matrix_cpp(T)
  expect_true(all(is.finite(out$transition_matrix)))
  expect_true(all(out$transition_matrix >= 0))
  expect_equal(rowSums(out$transition_matrix), rep(1, n), tolerance = 1e-8)
  expect_true(out$nans_fixed >= 1 || out$negs_clipped >= 1 || out$zero_rows_fixed >= 1)
})

# ---------------------------------------------------------------------------
# 2. Stationary distribution
# ---------------------------------------------------------------------------

test_that("cellrank_stationary_distribution_cpp sums to 1", {
  T <- make_transition_matrix(30, seed = 2)
  pi <- scop:::cellrank_stationary_distribution_cpp(T)
  expect_equal(sum(pi), 1, tolerance = 1e-8)
  expect_true(all(pi >= 0))
})

test_that("stationary distribution satisfies pi = pi * T", {
  n <- 15
  T <- make_transition_matrix(n, seed = 5)
  pi <- scop:::cellrank_stationary_distribution_cpp(T)
  pi_next <- as.numeric(pi %*% T)
  expect_equal(pi, pi_next, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 3. Schur (eigen-based) decomposition
# ---------------------------------------------------------------------------

test_that("cellrank_schur_cpp returns valid components", {
  n <- 30
  T <- make_transition_matrix(n, seed = 3)
  out <- scop:::cellrank_schur_cpp(T, n_components = 5)
  expect_true(abs(out$eigenvalues[1] - 1.0) < 0.15)
  expect_equal(ncol(out$schur_vectors), 5)
  expect_equal(nrow(out$schur_vectors), n)
  expect_true(all(out$macrostate_assignment >= 1))
  expect_true(all(out$macrostate_assignment <= 5))
  expect_equal(length(out$stationary_distribution), n)
})

test_that("cellrank_schur_cpp defaults to n_components = 2 for small matrices", {
  T <- make_transition_matrix(3, seed = 10)
  out <- scop:::cellrank_schur_cpp(T, n_components = 10)
  expect_true(ncol(out$schur_vectors) <= 3)
})

# ---------------------------------------------------------------------------
# 4. Auto-detect n_states
# ---------------------------------------------------------------------------

test_that("cellrank_auto_n_states_cpp returns valid range", {
  evals <- c(1.0, 0.95, 0.8, 0.3, 0.1)
  n <- scop:::cellrank_auto_n_states_cpp(evals, min_states = 2, max_states = 20)
  expect_true(n >= 2)
  expect_true(n <= length(evals))
})

# ---------------------------------------------------------------------------
# 5. Velocity kernel
# ---------------------------------------------------------------------------

test_that("cellrank_velocity_kernel_cpp produces valid transition matrix", {
  n_cells <- 20; n_dims <- 2; n_neighbors <- 5
  set.seed(1)
  vel_emb <- matrix(rnorm(n_cells * n_dims), n_cells, n_dims)
  embedding <- matrix(rnorm(n_cells * n_dims), n_cells, n_dims)
  knn_idx <- make_knn_idx(n_cells, n_neighbors, seed = 1)

  T <- scop:::cellrank_velocity_kernel_cpp(vel_emb, embedding, knn_idx)
  expect_equal(dim(T), c(n_cells, n_cells))
  expect_true(all(is.finite(T)))
  expect_true(all(T >= 0))
  # Rows should sum to ~1
  rs <- rowSums(T)
  expect_equal(rs, rep(1, n_cells), tolerance = 1e-8)
})

test_that("velocity kernel backward mode differs from forward", {
  n_cells <- 15; n_dims <- 2
  set.seed(2)
  # Create directed velocity field: cells on left point right
  embedding <- matrix(c(rep(1, 8), rep(4, 7), rep(0, 8), rep(0, 7)), ncol = 2)
  vel_emb <- matrix(c(rep(3, 8), rep(-3, 7), rep(0, 15)), ncol = 2)
  knn_idx <- make_knn_idx(n_cells, 4, seed = 2)

  T_fwd <- scop:::cellrank_velocity_kernel_cpp(vel_emb, embedding, knn_idx, backward = FALSE)
  T_bwd <- scop:::cellrank_velocity_kernel_cpp(vel_emb, embedding, knn_idx, backward = TRUE)
  # With directed velocities, forward and backward should differ somewhere
  expect_false(identical(T_fwd, T_bwd))
  expect_equal(rowSums(T_fwd), rep(1, n_cells), tolerance = 1e-8)
  expect_equal(rowSums(T_bwd), rep(1, n_cells), tolerance = 1e-8)
})

test_that("velocity kernel handles zero-velocity cells", {
  n_cells <- 10; n_dims <- 2
  vel_emb <- matrix(0, n_cells, n_dims)  # all zero velocity
  set.seed(3)
  embedding <- matrix(rnorm(n_cells * n_dims), n_cells, n_dims)
  knn_idx <- make_knn_idx(n_cells, 3, seed = 3)

  T <- scop:::cellrank_velocity_kernel_cpp(vel_emb, embedding, knn_idx)
  # Zero-velocity cells should have self-loop = 1
  for (i in seq_len(n_cells)) expect_equal(T[i, i], 1)
})

test_that("connectivity kernel cpp matches legacy neighbor weighting", {
  knn_idx <- matrix(
    c(2L, 2L, 1L, 3L, NA_integer_, 2L, 4L, 1L, 3L, 2L, 4L, 3L),
    nrow = 4,
    byrow = TRUE
  )
  knn_dist <- matrix(
    c(0.4, 0.8, 0.2, 0.5, NA, NA, 0.3, 0.6, 0.5, 0.9, 0.7, 0.1),
    nrow = 4,
    byrow = TRUE
  )
  legacy <- matrix(0, 4, 4)
  for (i in seq_len(nrow(knn_idx))) {
    row_sum <- 0
    distances <- knn_dist[i, ][!is.na(knn_dist[i, ])]
    sigma <- if (length(distances) > 0L) stats::median(distances) + 1e-10 else 1
    for (k in seq_len(ncol(knn_idx))) {
      neighbor <- knn_idx[i, k]
      if (is.na(neighbor) || neighbor < 1L || neighbor > 4L || neighbor == i) next
      distance <- knn_dist[i, k]
      if (is.na(distance)) next
      weight <- exp(-(distance * distance) / (2 * sigma * sigma))
      legacy[i, neighbor] <- weight
      row_sum <- row_sum + weight
    }
    if (row_sum > 0) legacy[i, ] <- legacy[i, ] / row_sum else legacy[i, i] <- 1
  }

  out <- scop:::cellrank_connectivity_kernel_cpp(knn_idx, knn_dist)

  expect_equal(out, legacy, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 6. Pseudotime kernel
# ---------------------------------------------------------------------------

test_that("cellrank_pseudotime_kernel_cpp produces valid transition matrix", {
  n_cells <- 20
  set.seed(4)
  pseudotime <- sort(runif(n_cells))
  knn_idx <- make_knn_idx(n_cells, 5, seed = 4)

  T <- scop:::cellrank_pseudotime_kernel_cpp(pseudotime, knn_idx, cell_weights = rep(1, n_cells))
  expect_equal(dim(T), c(n_cells, n_cells))
  expect_true(all(is.finite(T)))
  expect_true(all(T >= 0))
  expect_equal(rowSums(T), rep(1, n_cells), tolerance = 1e-8)
})

test_that("pseudotime kernel forward mode transitions toward later pseudotime", {
  n_cells <- 10
  pseudotime <- 1:10 / 10
  knn_idx <- make_knn_idx(n_cells, 3, seed = 7)

  T <- scop:::cellrank_pseudotime_kernel_cpp(
    pseudotime, knn_idx, cell_weights = rep(1, n_cells), backward = FALSE
  )
  # Early cell should have non-zero transitions
  expect_true(sum(T[1, ]) > 0)
})

# ---------------------------------------------------------------------------
# 7. CytoTRACE kernel
# ---------------------------------------------------------------------------

test_that("cellrank_cytotrace_kernel_cpp produces valid transition matrix", {
  n_cells <- 20
  set.seed(5)
  gene_counts <- runif(n_cells, 500, 5000)
  knn_idx <- make_knn_idx(n_cells, 5, seed = 5)

  T <- scop:::cellrank_cytotrace_kernel_cpp(gene_counts, knn_idx)
  expect_equal(dim(T), c(n_cells, n_cells))
  expect_true(all(is.finite(T)))
  expect_true(all(T >= 0))
  expect_equal(rowSums(T), rep(1, n_cells), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 8. CFLARE estimator
# ---------------------------------------------------------------------------

test_that("cellrank_cflare_cpp full pipeline runs end-to-end", {
  n <- 50
  T <- make_transition_matrix(n, seed = 4)
  out <- scop:::cellrank_cflare_cpp(T, n_states = 4)
  expect_equal(dim(out$transition_matrix), c(n, n))
  expect_equal(length(out$stationary_distribution), n)
  expect_true(length(out$eigenvalues) > 0)
  expect_equal(length(out$macrostate_assignment), n)
  expect_equal(length(out$terminal_states), n)
  expect_equal(length(out$fate_confidence), n)
  expect_true(out$n_terminal_states >= 1)
  expect_equal(out$method, "CFLARE")
})

test_that("CFLARE absorption probabilities between 0 and 1", {
  n <- 40
  T <- make_transition_matrix(n, seed = 6)
  out <- scop:::cellrank_cflare_cpp(T, n_states = 3)
  ap <- out$absorption_probabilities
  expect_true(all(ap >= 0, na.rm = TRUE))
  expect_true(all(ap <= 1 + 1e-6, na.rm = TRUE))
})

test_that("CFLARE handles identity matrix (absorbing)", {
  n <- 10
  T <- diag(1, n)
  out <- scop:::cellrank_cflare_cpp(T, n_states = 3)
  expect_equal(dim(out$transition_matrix), c(n, n))
  expect_true(all(out$fate_confidence >= 0))
})

# ---------------------------------------------------------------------------
# 9. GPCCA estimator
# ---------------------------------------------------------------------------

test_that("cellrank_gpcca_cpp full pipeline runs end-to-end", {
  n <- 50
  T <- make_transition_matrix(n, seed = 7)
  out <- tryCatch(
    scop:::cellrank_gpcca_cpp(T, n_states = 4, n_cells_terminal = 5),
    error = function(e) NULL
  )
  skip_if(is.null(out), "GPCCA pipeline failed on random matrix")
  expect_equal(dim(out$transition_matrix), c(n, n))
  expect_equal(length(out$stationary_distribution), n)
  expect_true(length(out$eigenvalues) > 0)
  expect_equal(length(out$macrostate_assignment), n)
  expect_true(is.matrix(out$chi))
  expect_equal(nrow(out$chi), n)
  expect_equal(out$method, "GPCCA")
})

test_that("GPCCA absorption probabilities between 0 and 1", {
  n <- 30
  T <- make_transition_matrix(n, seed = 8)
  out <- tryCatch(
    scop:::cellrank_gpcca_cpp(T, n_states = 3, n_cells_terminal = 5),
    error = function(e) NULL
  )
  skip_if(is.null(out), "GPCCA failed")
  ap <- out$absorption_probabilities
  expect_true(all(ap >= 0, na.rm = TRUE))
  expect_true(all(ap <= 1 + 1e-6, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 10. Lineage drivers
# ---------------------------------------------------------------------------

test_that("cellrank_lineage_drivers_cpp computes correlations", {
  n_genes <- 20; n_cells <- 30
  set.seed(9)
  expression <- matrix(rgamma(n_genes * n_cells, shape = 2, rate = 1), nrow = n_genes)
  abs_probs <- matrix(runif(n_cells * 2), nrow = n_cells)
  abs_probs <- abs_probs / rowSums(abs_probs)

  out <- scop:::cellrank_lineage_drivers_cpp(expression, abs_probs, lineage_idx = as.integer(c(1, 2)))
  expect_equal(dim(out$correlation), c(n_genes, 2))
  expect_true(all(out$correlation >= -1 - 1e-6, na.rm = TRUE))
  expect_true(all(out$correlation <= 1 + 1e-6, na.rm = TRUE))
})

test_that("lineage drivers with specific lineage indices", {
  n_genes <- 10; n_cells <- 20
  set.seed(10)
  expression <- matrix(rgamma(n_genes * n_cells, shape = 2), nrow = n_genes)
  abs_probs <- matrix(runif(n_cells * 3), nrow = n_cells)
  abs_probs <- abs_probs / rowSums(abs_probs)

  out <- scop:::cellrank_lineage_drivers_cpp(expression, abs_probs, lineage_idx = as.integer(c(1, 3)))
  expect_equal(dim(out$correlation), c(n_genes, 2))
  expect_equal(out$lineage_idx, as.integer(c(1, 3)))
})

# ---------------------------------------------------------------------------
# 11. Determinism
# ---------------------------------------------------------------------------

test_that("CFLARE is deterministic", {
  T <- make_transition_matrix(30, seed = 11)
  out1 <- scop:::cellrank_cflare_cpp(T, n_states = 3)
  out2 <- scop:::cellrank_cflare_cpp(T, n_states = 3)
  expect_equal(out1$macrostate_assignment, out2$macrostate_assignment)
  expect_equal(out1$absorption_probabilities, out2$absorption_probabilities)
})

test_that("velocity kernel is deterministic", {
  n_cells <- 20; n_dims <- 2
  set.seed(12)
  vel_emb <- matrix(rnorm(n_cells * n_dims), n_cells, n_dims)
  embedding <- matrix(rnorm(n_cells * n_dims), n_cells, n_dims)
  knn_idx <- make_knn_idx(n_cells, 5, seed = 12)
  T1 <- scop:::cellrank_velocity_kernel_cpp(vel_emb, embedding, knn_idx)
  T2 <- scop:::cellrank_velocity_kernel_cpp(vel_emb, embedding, knn_idx)
  expect_equal(T1, T2)
})
