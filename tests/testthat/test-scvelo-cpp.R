# Tests for scVelo C++ backend (scvelo_stochastic_embedding_cpp)
#
# Structure:
#   1. Helper: mock spliced/unspliced/KNN/embedding data
#   2. Basic structural tests (output shape, valid fields)
#   3. Edge case: all genes zero
#   4. Edge case: unspliced = spliced (gamma should be ~1)
#   5. Edge case: single cell
#   6. Edge case: very few neighbors
#   7. Numerical stability
#   8. Determinism
#   9. Input validation

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_scvelo_mock <- function(
  n_genes = 10,
  n_cells = 20,
  n_neighbors = 4,
  n_dims = 2,
  seed = 42
) {
  set.seed(seed)
  # Spliced expression (genes × cells)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 2, rate = 1),
    nrow = n_genes,
    ncol = n_cells
  )
  # Unspliced: correlated with spliced + noise (beta ~ 0.5)
  beta <- 0.5
  unspliced <- beta * spliced + matrix(
    stats::rnorm(n_genes * n_cells, sd = 0.2),
    nrow = n_genes,
    ncol = n_cells
  )
  unspliced <- pmax(unspliced, 0)

  # KNN index (cells × neighbors), 1-based, self excluded
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, n_neighbors))
  }

  # Embedding (cells × dims), e.g. UMAP
  embedding <- matrix(
    stats::rnorm(n_cells * n_dims, sd = 2),
    nrow = n_cells,
    ncol = n_dims
  )

  list(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding,
    n_genes = n_genes,
    n_cells = n_cells,
    n_neighbors = n_neighbors,
    n_dims = n_dims,
    beta = beta
  )
}

# ---------------------------------------------------------------------------
# 1. Basic structural tests
# ---------------------------------------------------------------------------

test_that("scvelo_stochastic_embedding_cpp returns a named list with correct fields", {
  dat <- make_scvelo_mock()
  out <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    knn_idx = dat$knn_idx,
    embedding = dat$embedding
  )

  expect_type(out, "list")
  # Must contain at least: velocity_embedding, confidence, velocity_length, gamma
  required <- c("velocity_embedding", "confidence", "velocity_length", "gamma")
  expect_true(all(required %in% names(out)))

  # velocity_embedding: cells × dims
  expect_true(is.matrix(out$velocity_embedding))
  expect_equal(dim(out$velocity_embedding), c(dat$n_cells, dat$n_dims))
  expect_true(all(is.finite(out$velocity_embedding)))

  # confidence: length = n_cells, [0, 1]
  expect_type(out$confidence, "double")
  expect_length(out$confidence, dat$n_cells)
  expect_true(all(out$confidence >= 0, na.rm = TRUE))
  expect_true(all(out$confidence <= 1, na.rm = TRUE))

  # velocity_length: length = n_cells, non-negative
  expect_type(out$velocity_length, "double")
  expect_length(out$velocity_length, dat$n_cells)
  expect_true(all(out$velocity_length >= 0, na.rm = TRUE))

  # gamma: length = n_genes, non-negative
  expect_type(out$gamma, "double")
  expect_length(out$gamma, dat$n_genes)
  expect_true(all(out$gamma >= 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 2. Edge case: all zeros
# ---------------------------------------------------------------------------

test_that("scvelo handles all-zero expression", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(0, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(0, nrow = n_genes, ncol = n_cells)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 3))
  }
  set.seed(1)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_true(all(out$velocity_embedding == 0))
  expect_true(all(out$confidence == 0))
  expect_true(all(out$velocity_length == 0))
  expect_true(all(out$gamma == 0))
})

# ---------------------------------------------------------------------------
# 3. Edge case: unspliced = spliced (gamma should converge to ~1)
# ---------------------------------------------------------------------------

test_that("scvelo estimates gamma ~ 1 when unspliced = spliced", {
  n_genes <- 20
  n_cells <- 30
  set.seed(1)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 3, rate = 1),
    nrow = n_genes, ncol = n_cells
  )
  unspliced <- spliced  # identical
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 4)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 4))
  }
  set.seed(1)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  # Gamma should be close to 1 for all genes
  expect_equal(out$gamma, rep(1, n_genes), tolerance = 0.1)
  # Residual velocity should be near zero
  expect_true(mean(out$velocity_length) < 0.1)
})

# ---------------------------------------------------------------------------
# 4. Edge case: single cell
# ---------------------------------------------------------------------------

test_that("scvelo handles single cell gracefully", {
  n_genes <- 5
  n_cells <- 1
  spliced <- matrix(runif(n_genes), nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(runif(n_genes), nrow = n_genes, ncol = n_cells)
  # KNN: no valid neighbors since exclude_self=TRUE
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  embedding <- matrix(runif(2), nrow = 1, ncol = 2)

  expect_error(
    scop:::scvelo_stochastic_embedding_cpp(
      spliced = spliced,
      unspliced = unspliced,
      knn_idx = knn_idx,
      embedding = embedding
    ),
    NA
  )
})

# ---------------------------------------------------------------------------
# 5. Edge case: few neighbors
# ---------------------------------------------------------------------------

test_that("scvelo works with 1 neighbor", {
  n_genes <- 5
  n_cells <- 10
  set.seed(1)
  spliced <- matrix(runif(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(runif(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 1)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sample(candidates, 1)
  }
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_equal(dim(out$velocity_embedding), c(10, 2))
  expect_true(all(is.finite(out$velocity_embedding)))
})

# ---------------------------------------------------------------------------
# 6. Determinism
# ---------------------------------------------------------------------------

test_that("scvelo_stochastic_embedding_cpp is deterministic", {
  dat <- make_scvelo_mock(seed = 42)
  out1 <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    knn_idx = dat$knn_idx,
    embedding = dat$embedding
  )
  out2 <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    knn_idx = dat$knn_idx,
    embedding = dat$embedding
  )
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 7. Gamma directionality: larger spliced → smaller gamma when unspliced low
# ---------------------------------------------------------------------------

test_that("scvelo gamma reflects relationship between spliced and unspliced", {
  n_genes <- 10
  n_cells <- 30
  set.seed(1)
  spliced <- matrix(runif(n_genes * n_cells, 1, 10), nrow = n_genes, ncol = n_cells)
  # Make first 5 genes have high unspliced/spliced ratio (high gamma),
  # last 5 genes have low ratio
  unspliced <- spliced
  unspliced[1:5, ] <- unspliced[1:5, ] * 2.0    # gamma approx 2
  unspliced[6:10, ] <- unspliced[6:10, ] * 0.1  # gamma approx 0.1

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 4)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 4))
  }
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  # First 5 genes should have higher gamma than last 5
  expect_true(mean(out$gamma[1:5]) > mean(out$gamma[6:10]) + 0.5)
})

# ---------------------------------------------------------------------------
# 8. Input validation
# ---------------------------------------------------------------------------

test_that("scvelo rejects mismatched spliced/unspliced dimensions", {
  dat <- make_scvelo_mock()
  expect_error(
    scop:::scvelo_stochastic_embedding_cpp(
      spliced = dat$spliced[1:5, , drop = FALSE],  # different n_genes
      unspliced = dat$unspliced,
      knn_idx = dat$knn_idx,
      embedding = dat$embedding
    ),
    "spliced and unspliced must have identical dimensions"
  )
})

test_that("scvelo rejects mismatched knn_idx rows vs cells", {
  dat <- make_scvelo_mock()
  expect_error(
    scop:::scvelo_stochastic_embedding_cpp(
      spliced = dat$spliced,
      unspliced = dat$unspliced,
      knn_idx = dat$knn_idx[1:5, , drop = FALSE],
      embedding = dat$embedding
    ),
    "knn_idx"
  )
})

test_that("scvelo rejects mismatched embedding rows vs cells", {
  dat <- make_scvelo_mock()
  expect_error(
    scop:::scvelo_stochastic_embedding_cpp(
      spliced = dat$spliced,
      unspliced = dat$unspliced,
      knn_idx = dat$knn_idx,
      embedding = dat$embedding[1:5, , drop = FALSE]
    ),
    "embedding"
  )
})

# ---------------------------------------------------------------------------
# 9. Velocity embedding direction summary
# ---------------------------------------------------------------------------

test_that("scvelo embedding output is finite and well-scaled", {
  n_genes <- 30
  n_cells <- 50
  set.seed(1)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 2, rate = 1),
    nrow = n_genes, ncol = n_cells
  )
  unspliced <- 0.3 * spliced + matrix(
    rnorm(n_genes * n_cells, sd = 0.1), nrow = n_genes, ncol = n_cells
  )
  unspliced <- pmax(unspliced, 0)

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 6)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, min(6, length(candidates))))
  }
  set.seed(99)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scop:::scvelo_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_true(all(is.finite(out$velocity_embedding)))
  expect_true(all(is.finite(out$confidence)))
  expect_true(all(is.finite(out$velocity_length)))
  expect_true(all(is.finite(out$gamma)))

  # At least some velocity signal should be present
  expect_gt(mean(out$velocity_length), 0)
  # Some confidence should be > 0
  expect_gt(max(out$confidence), 0)
})
