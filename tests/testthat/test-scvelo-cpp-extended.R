# Tests for scVelo C++ backend — extended coverage
#
# Existing test-scvelo-cpp.R covers: scanpy_stochastic_embedding_cpp (9 tests)
# This file adds coverage for independently-testable pipeline components:
#   1. scanpy_filter_genes_cpp — gene filtering
#   2. scanpy_normalize_log_cpp — per-cell normalization + log1p
#   3. scanpy_moments_cpp — first-order moments (KNN smoothing)
#
# NOTE: scanpy_deterministic_cpp, scanpy_velocity_confidence_cpp,
# scanpy_velocity_transition_cpp, and scanpy_velocity_graph_cpp require
# dependent intermediate computations (Ms, residual, embedding) and are
# better tested via integration the run_scanpy_cpp wrapper or benchmark.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_scanpy_data <- function(
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
# 1. scanpy_filter_genes_cpp
# ---------------------------------------------------------------------------

test_that("scanpy_filter_genes_cpp returns 0/1 indicator vector", {
  dat <- make_scanpy_data()
  out <- scanpy_filter_genes_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    min_counts = 2,
    min_counts_u = 1
  )
  expect_type(out, "integer")
  expect_true(all(out %in% c(0, 1)))
  expect_length(out, dat$n_genes)
})

test_that("scanpy_filter_genes_cpp keeps all genes with high counts", {
  n_genes <- 10
  n_cells <- 20
  spliced <- matrix(100, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(50, nrow = n_genes, ncol = n_cells)

  out <- scanpy_filter_genes_cpp(
    spliced = spliced,
    unspliced = unspliced,
    min_counts = 10,
    min_counts_u = 5
  )
  expect_equal(length(out), n_genes)
  expect_true(all(out == 1))
})

test_that("scanpy_filter_genes_cpp filters all-zero spliced genes", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(c(rep(0, n_cells), rep(100, n_cells * (n_genes - 1))),
    nrow = n_genes, ncol = n_cells, byrow = TRUE
  )
  unspliced <- matrix(100, nrow = n_genes, ncol = n_cells)

  out <- scanpy_filter_genes_cpp(spliced, unspliced, min_counts = 1, min_counts_u = 5)
  expect_equal(out[1], 0) # first gene filtered (zero spliced)
  expect_true(all(out[-1] == 1))
})

test_that("scanpy_filter_genes_cpp filters all-zero unspliced genes", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(100, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(c(rep(0, n_cells), rep(100, n_cells * (n_genes - 1))),
    nrow = n_genes, ncol = n_cells, byrow = TRUE
  )

  out <- scanpy_filter_genes_cpp(spliced, unspliced, min_counts = 5, min_counts_u = 1)
  expect_equal(out[1], 0) # first gene filtered (zero unspliced)
  expect_true(all(out[-1] == 1))
})

test_that("scanpy_filter_genes_cpp is deterministic", {
  dat <- make_scanpy_data(seed = 5)
  out1 <- scanpy_filter_genes_cpp(dat$spliced, dat$unspliced, 2, 1)
  out2 <- scanpy_filter_genes_cpp(dat$spliced, dat$unspliced, 2, 1)
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 2. scanpy_normalize_log_cpp
# ---------------------------------------------------------------------------

test_that("scanpy_normalize_log_cpp returns correct structure", {
  dat <- make_scanpy_data()
  norm <- scanpy_normalize_log_cpp(
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

test_that("scanpy_normalize_log_cpp handles all-zero input", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(0, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(0, nrow = n_genes, ncol = n_cells)

  norm <- scanpy_normalize_log_cpp(spliced, unspliced)
  expect_equal(norm$spliced_norm, matrix(0, n_genes, n_cells))
  expect_equal(norm$unspliced_norm, matrix(0, n_genes, n_cells))
})

test_that("scanpy_normalize_log_cpp is deterministic", {
  dat <- make_scanpy_data(seed = 5)
  out1 <- scanpy_normalize_log_cpp(dat$spliced, dat$unspliced)
  out2 <- scanpy_normalize_log_cpp(dat$spliced, dat$unspliced)
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 3. scanpy_moments_cpp
# ---------------------------------------------------------------------------

test_that("scanpy_moments_cpp returns correct structure", {
  dat <- make_scanpy_data()
  norm <- scanpy_normalize_log_cpp(dat$spliced, dat$unspliced)

  moments <- scanpy_moments_cpp(
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

test_that("scanpy_moments_cpp includes self in neighborhood average", {
  # With 1 cell and 0 neighbors, moments should equal input
  spliced <- matrix(runif(10), nrow = 10, ncol = 1)
  unspliced <- matrix(runif(10), nrow = 10, ncol = 1)
  # KNN with no valid neighbors (all NA)
  knn_idx <- matrix(NA_integer_, nrow = 1, ncol = 2)

  moments <- scanpy_moments_cpp(spliced, unspliced, knn_idx)
  # Should equal the input (self only)
  expect_equal(moments$Ms[, 1], spliced[, 1], tolerance = 1e-10)
  expect_equal(moments$Mu[, 1], unspliced[, 1], tolerance = 1e-10)
})

test_that("scanpy_moments_cpp is deterministic", {
  dat <- make_scanpy_data(seed = 7)
  norm <- scanpy_normalize_log_cpp(dat$spliced, dat$unspliced)

  m1 <- scanpy_moments_cpp(norm$spliced_norm, norm$unspliced_norm, dat$knn_idx)
  m2 <- scanpy_moments_cpp(norm$spliced_norm, norm$unspliced_norm, dat$knn_idx)
  expect_equal(m1, m2)
})

test_that("connectivity moments can skip unused second-order matrices", {
  set.seed(20260810)
  spliced <- matrix(stats::rexp(8L * 12L), 8L, 12L)
  unspliced <- matrix(stats::rexp(8L * 12L), 8L, 12L)
  knn_idx <- t(vapply(seq_len(12L), function(cell) {
    sample(setdiff(seq_len(12L), cell), 4L)
  }, integer(4L)))

  full <- scanpy_moments_connectivities_cpp(
    spliced,
    unspliced,
    knn_idx,
    compute_second_order = TRUE
  )
  first_order <- scanpy_moments_connectivities_cpp(
    spliced,
    unspliced,
    knn_idx,
    compute_second_order = FALSE
  )

  expect_equal(first_order$Ms, full$Ms, tolerance = 0)
  expect_equal(first_order$Mu, full$Mu, tolerance = 0)
  expect_null(first_order$Mss)
  expect_null(first_order$Mus)
  expect_equal(full$Mss, scanpy_second_order_moments_cpp(
    spliced,
    unspliced,
    knn_idx
  )$Mss, tolerance = 0)
})

test_that("deterministic scVelo fit matches extreme-quantile semantics", {
  set.seed(20260810)
  n_genes <- 24L
  n_cells <- 80L
  Ms <- matrix(stats::rgamma(n_genes * n_cells, shape = 1.5), n_genes)
  Ms[matrix(stats::runif(length(Ms)), n_genes) < 0.55] <- 0
  Mu <- matrix(pmax(
    0,
    Ms * stats::runif(n_genes, 0.25, 0.9) +
      matrix(stats::rnorm(length(Ms), sd = 0.18), n_genes)
  ), nrow = n_genes)
  Mu[1L, ] <- 0
  Ms[2L, ] <- 0
  Mu[3:5, ] <- matrix(stats::rexp(3L * n_cells), 3L)
  knn_idx <- t(vapply(seq_len(n_cells), function(cell) {
    sample(setdiff(seq_len(n_cells), cell), 8L)
  }, integer(8L)))
  embedding <- matrix(stats::rnorm(n_cells * 2L), n_cells, 2L)

  reference_fit <- lapply(seq_len(n_genes), function(gene) {
    s <- Ms[gene, ]
    u <- Mu[gene, ]
    normalized <- s / max(max(s), 1e-3) + u / max(max(u), 1e-3)
    keep <- normalized >= stats::quantile(normalized, 0.95, type = 7)
    denominator <- sum(s[keep]^2)
    gamma <- if (denominator > 0) sum(s[keep] * u[keep]) / denominator else 0
    total <- sum((u - mean(u))^2)
    residual_sum <- sum((u - gamma * s)^2)
    r2 <- if (total > 0) 1 - residual_sum / total else 0
    c(gamma = gamma, r2 = r2)
  })
  reference_fit <- do.call(rbind, reference_fit)
  velocity_genes <- reference_fit[, "r2"] > 0.01 &
    reference_fit[, "gamma"] > 0.01 &
    rowSums(Ms > 0) > 0 & rowSums(Mu > 0) > 0
  if (sum(velocity_genes) < 2L) {
    velocity_genes <- reference_fit[, "r2"] >
      stats::quantile(reference_fit[, "r2"], 0.80, type = 7)
  }
  expect_gt(sum(velocity_genes), 1L)
  expect_true(any(!velocity_genes))

  observed <- scanpy_deterministic_cpp(
    Ms = Ms,
    Mu = Mu,
    knn_idx = knn_idx,
    embedding = embedding,
    fit_offset = FALSE,
    perc = 95
  )

  expect_equal(observed$gamma, reference_fit[, "gamma"], tolerance = 1e-12)
  expect_equal(observed$r2, reference_fit[, "r2"], tolerance = 1e-12)
  expect_identical(as.logical(observed$velocity_genes), velocity_genes)
  expect_equal(
    observed$residual,
    Mu - reference_fit[, "gamma"] * Ms,
    tolerance = 1e-12
  )

  confidence <- scanpy_velocity_confidence_cpp(
    Ms = Ms[velocity_genes, , drop = FALSE],
    residual = observed$residual[velocity_genes, , drop = FALSE],
    knn_idx = knn_idx
  )
  centered <- sweep(
    observed$residual[velocity_genes, , drop = FALSE],
    2L,
    colMeans(observed$residual[velocity_genes, , drop = FALSE]),
    "-"
  )
  expect_equal(
    confidence$velocity_length,
    round(sqrt(colSums(centered^2)), 2),
    tolerance = 1e-12
  )
})

# ---------------------------------------------------------------------------
# 4. Input validation
# ---------------------------------------------------------------------------

test_that("scanpy_normalize_log_cpp rejects mismatched dimensions", {
  expect_error(
    scanpy_normalize_log_cpp(
      spliced = matrix(1, 5, 10),
      unspliced = matrix(1, 3, 10)
    ),
    "identical dimensions"
  )
})

test_that("scanpy_moments_cpp rejects mismatched knn_idx rows", {
  spliced <- matrix(1, 5, 10)
  unspliced <- matrix(1, 5, 10)
  knn_idx <- matrix(1L, nrow = 5, ncol = 3) # only 5 rows, need 10

  expect_error(
    scanpy_moments_cpp(spliced, unspliced, knn_idx),
    "knn_idx",
    ignore.case = TRUE
  )
})
