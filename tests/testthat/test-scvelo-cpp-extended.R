# Tests for scVelo C++ backend — extended coverage
#
# Existing test-scvelo-cpp.R covers: scvelo_stochastic_embedding_cpp (9 tests)
# This file adds coverage for independently-testable pipeline components:
#   1. scvelo_filter_genes_cpp — gene filtering
#   2. scvelo_normalize_log_cpp — per-cell normalization + log1p
#   3. scvelo_moments_cpp — first-order moments (KNN smoothing)
#
# NOTE: scvelo_deterministic_cpp, scvelo_velocity_confidence_cpp,
# scvelo_velocity_transition_cpp, and scvelo_velocity_graph_cpp require
# dependent intermediate computations (Ms, residual, embedding) and are
# better tested via integration the run_scvelo_cpp wrapper or benchmark.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_scvelo_data <- function(
  n_genes = 20,
  n_cells = 30,
  n_neighbors = 8,
  seed = 42
) {
  set.seed(seed)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 2, rate = 1),
    nrow = n_genes,
    ncol = n_cells
  )
  unspliced <- 0.4 * spliced + matrix(
    stats::rnorm(n_genes * n_cells, sd = 0.3),
    nrow = n_genes,
    ncol = n_cells
  )
  unspliced <- pmax(unspliced, 0)

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, n_neighbors))
  }

  list(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    n_genes = n_genes,
    n_cells = n_cells,
    n_neighbors = n_neighbors
  )
}

# ---------------------------------------------------------------------------
# 1. scvelo_filter_genes_cpp
# ---------------------------------------------------------------------------

test_that("scvelo_filter_genes_cpp returns 0/1 indicator vector", {
  dat <- make_scvelo_data()
  out <- scvelo_filter_genes_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    min_counts = 2,
    min_counts_u = 1
  )
  expect_type(out, "integer")
  expect_true(all(out %in% c(0, 1)))
  expect_length(out, dat$n_genes)
})

test_that("scvelo_filter_genes_cpp keeps all genes with high counts", {
  n_genes <- 10
  n_cells <- 20
  spliced <- matrix(100, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(50, nrow = n_genes, ncol = n_cells)

  out <- scvelo_filter_genes_cpp(
    spliced = spliced,
    unspliced = unspliced,
    min_counts = 10,
    min_counts_u = 5
  )
  expect_equal(length(out), n_genes)
  expect_true(all(out == 1))
})

test_that("scvelo_filter_genes_cpp filters all-zero spliced genes", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(c(rep(0, n_cells), rep(100, n_cells * (n_genes - 1))),
                    nrow = n_genes, ncol = n_cells, byrow = TRUE)
  unspliced <- matrix(100, nrow = n_genes, ncol = n_cells)

  out <- scvelo_filter_genes_cpp(spliced, unspliced, min_counts = 1, min_counts_u = 5)
  expect_equal(out[1], 0)  # first gene filtered (zero spliced)
  expect_true(all(out[-1] == 1))
})

test_that("scvelo_filter_genes_cpp filters all-zero unspliced genes", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(100, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(c(rep(0, n_cells), rep(100, n_cells * (n_genes - 1))),
                      nrow = n_genes, ncol = n_cells, byrow = TRUE)

  out <- scvelo_filter_genes_cpp(spliced, unspliced, min_counts = 5, min_counts_u = 1)
  expect_equal(out[1], 0)  # first gene filtered (zero unspliced)
  expect_true(all(out[-1] == 1))
})

test_that("scvelo_filter_genes_cpp is deterministic", {
  dat <- make_scvelo_data(seed = 5)
  out1 <- scvelo_filter_genes_cpp(dat$spliced, dat$unspliced, 2, 1)
  out2 <- scvelo_filter_genes_cpp(dat$spliced, dat$unspliced, 2, 1)
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 2. scvelo_normalize_log_cpp
# ---------------------------------------------------------------------------

test_that("scvelo_normalize_log_cpp returns correct structure", {
  dat <- make_scvelo_data()
  norm <- scvelo_normalize_log_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced
  )
  expect_type(norm, "list")
  expect_true(all(c("spliced_norm", "unspliced_norm") %in% names(norm)))
  expect_equal(dim(norm$spliced_norm), dim(dat$spliced))
  expect_true(all(is.finite(norm$spliced_norm)))
  expect_true(all(is.finite(norm$unspliced_norm)))
  # Values should be non-negative (log1p of non-negative)
  expect_true(all(norm$spliced_norm >= 0))
  expect_true(all(norm$unspliced_norm >= 0))
})

test_that("scvelo_normalize_log_cpp handles all-zero input", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(0, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(0, nrow = n_genes, ncol = n_cells)

  norm <- scvelo_normalize_log_cpp(spliced, unspliced)
  expect_equal(norm$spliced_norm, matrix(0, n_genes, n_cells))
  expect_equal(norm$unspliced_norm, matrix(0, n_genes, n_cells))
})

test_that("scvelo_normalize_log_cpp is deterministic", {
  dat <- make_scvelo_data(seed = 5)
  out1 <- scvelo_normalize_log_cpp(dat$spliced, dat$unspliced)
  out2 <- scvelo_normalize_log_cpp(dat$spliced, dat$unspliced)
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 3. scvelo_moments_cpp
# ---------------------------------------------------------------------------

test_that("scvelo_moments_cpp returns correct structure", {
  dat <- make_scvelo_data()
  norm <- scvelo_normalize_log_cpp(dat$spliced, dat$unspliced)

  moments <- scvelo_moments_cpp(
    spliced = norm$spliced_norm,
    unspliced = norm$unspliced_norm,
    knn_idx = dat$knn_idx
  )
  expect_type(moments, "list")
  expect_true(all(c("Ms", "Mu") %in% names(moments)))
  expect_equal(dim(moments$Ms), dim(dat$spliced))
  expect_equal(dim(moments$Mu), dim(dat$spliced))
  expect_true(all(is.finite(moments$Ms)))
  expect_true(all(is.finite(moments$Mu)))
})

test_that("scvelo_moments_cpp includes self in neighborhood average", {
  # With 1 cell and 0 neighbors, moments should equal input
  spliced <- matrix(runif(10), nrow = 10, ncol = 1)
  unspliced <- matrix(runif(10), nrow = 10, ncol = 1)
  # KNN with no valid neighbors (all NA)
  knn_idx <- matrix(NA_integer_, nrow = 1, ncol = 2)

  moments <- scvelo_moments_cpp(spliced, unspliced, knn_idx)
  # Should equal the input (self only)
  expect_equal(moments$Ms[, 1], spliced[, 1], tolerance = 1e-10)
  expect_equal(moments$Mu[, 1], unspliced[, 1], tolerance = 1e-10)
})

test_that("scvelo_moments_cpp is deterministic", {
  dat <- make_scvelo_data(seed = 7)
  norm <- scvelo_normalize_log_cpp(dat$spliced, dat$unspliced)

  m1 <- scvelo_moments_cpp(norm$spliced_norm, norm$unspliced_norm, dat$knn_idx)
  m2 <- scvelo_moments_cpp(norm$spliced_norm, norm$unspliced_norm, dat$knn_idx)
  expect_equal(m1, m2)
})

# ---------------------------------------------------------------------------
# 4. Input validation
# ---------------------------------------------------------------------------

test_that("scvelo_normalize_log_cpp rejects mismatched dimensions", {
  expect_error(
    scvelo_normalize_log_cpp(
      spliced = matrix(1, 5, 10),
      unspliced = matrix(1, 3, 10)
    ),
    "identical dimensions"
  )
})

test_that("scvelo_moments_cpp rejects mismatched knn_idx rows", {
  spliced <- matrix(1, 5, 10)
  unspliced <- matrix(1, 5, 10)
  knn_idx <- matrix(1L, nrow = 5, ncol = 3)  # only 5 rows, need 10

  expect_error(
    scvelo_moments_cpp(spliced, unspliced, knn_idx),
    "knn_idx"
  )
})
