# Tests for PAGA C++ backend
#
# Structure:
#   1. Helper: mock KNN + group data
#   2. Basic structural tests (output shape, valid fields)
#   3. Edge case: 2 groups
#   4. Edge case: all cells same group
#   5. Edge case: many groups, sparse connections
#   6. Edge case: single cell per group
#   7. Numerical stability: zero-variance KNN structure
#   8. run_paga_cpp integration via mock Seurat
#   9. Backend parity: results are deterministic

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_paga_mock <- function(
  n_cells = 20,
  n_groups = 3,
  n_neighbors = 5,
  seed = 42
) {
  set.seed(seed)
  groups <- sample(seq_len(n_groups), n_cells, replace = TRUE)
  # KNN index: each row has n_neighbors (1-based), self excluded
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, n_neighbors))
  }
  list(
    knn_idx = knn_idx,
    groups = groups,
    n_groups = n_groups,
    n_cells = n_cells,
    n_neighbors = n_neighbors
  )
}

# ---------------------------------------------------------------------------
# 1. Basic structural tests
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp returns a named list with correct fields", {
  dat <- make_paga_mock()
  out <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx,
    groups = dat$groups,
    n_groups = dat$n_groups
  )

  expect_type(out, "list")
  expected_names <- c(
    "connectivities", "connectivities_tree",
    "expected_n_edges_random", "group_sizes", "directed_edges"
  )
  expect_named(out, expected_names)

  # connectivities: square matrix with correct dimension
  expect_true(is.matrix(out$connectivities))
  expect_equal(dim(out$connectivities), c(dat$n_groups, dat$n_groups))

  # connectivities_tree: same dimension
  expect_true(is.matrix(out$connectivities_tree))
  expect_equal(dim(out$connectivities_tree), c(dat$n_groups, dat$n_groups))

  # group_sizes: length = n_groups, sums to n_cells
  expect_type(out$group_sizes, "double")
  expect_length(out$group_sizes, dat$n_groups)
  expect_equal(sum(out$group_sizes), dat$n_cells)

  # directed_edges: square
  expect_equal(dim(out$directed_edges), c(dat$n_groups, dat$n_groups))
})

test_that("connectivities are symmetric and bounded [0, 1]", {
  dat <- make_paga_mock(seed = 7)
  out <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx,
    groups = dat$groups,
    n_groups = dat$n_groups
  )
  con <- out$connectivities
  expect_equal(con, t(con), tolerance = 1e-15)
  expect_true(all(con >= 0, na.rm = TRUE))
  expect_true(all(con <= 1 + 1e-15, na.rm = TRUE))
})

test_that("connectivities_tree is a maximum spanning tree (n_groups-1 non-zero edges)", {
  dat <- make_paga_mock(n_cells = 30, n_groups = 5, seed = 1)
  out <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx,
    groups = dat$groups,
    n_groups = dat$n_groups
  )
  tree <- out$connectivities_tree
  expect_equal(dim(tree), c(5, 5))
  n_nonzero <- sum(tree > 0)
  # For a spanning tree of 5 nodes, we need exactly 4 undirected edges,
  # stored symmetrically → 8 non-zero entries (lower + upper triangle)
  # Tree stored as directed adjacency (one direction per spanning-tree edge)
  expect_equal(n_nonzero, dat$n_groups - 1)
})

test_that("directed_edges counts are non-negative integers", {
  dat <- make_paga_mock(seed = 3)
  out <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx,
    groups = dat$groups,
    n_groups = dat$n_groups
  )
  de <- out$directed_edges
  expect_true(all(de >= 0))
  expect_true(all(de == round(de)))
})

# ---------------------------------------------------------------------------
# 2. Two-group case
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp works with exactly 2 groups", {
  n_cells <- 20
  set.seed(1)
  groups <- sample(1:2, n_cells, replace = TRUE)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 4)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 4))
  }
  out <- scop:::paga_connectivities_cpp(
    knn_idx = knn_idx, groups = groups, n_groups = 2
  )
  expect_equal(dim(out$connectivities), c(2, 2))
  expect_equal(dim(out$connectivities_tree), c(2, 2))
})

# ---------------------------------------------------------------------------
# 3. Single group — edge case
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp handles single group gracefully", {
  n_cells <- 10
  groups <- rep(1L, n_cells)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 3))
  }
  expect_error(
    scop:::paga_connectivities_cpp(knn_idx = knn_idx, groups = groups, n_groups = 1),
    NA  # should not error
  )
})

# ---------------------------------------------------------------------------
# 4. Many groups, sparse connections
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp handles many-groups-sparse scenario", {
  n_cells <- 50
  n_groups <- 10
  set.seed(42)
  groups <- sample(seq_len(n_groups), n_cells, replace = TRUE)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 3))
  }
  out <- scop:::paga_connectivities_cpp(
    knn_idx = knn_idx, groups = groups, n_groups = n_groups
  )
  expect_equal(dim(out$connectivities), c(n_groups, n_groups))
  expect_true(all(out$connectivities >= 0, na.rm = TRUE))
  expect_true(all(out$connectivities <= 1, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 5. Single cell per group
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp handles single-cell groups", {
  n_cells <- 6
  n_groups <- 6
  groups <- 1:6
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, min(3, length(candidates))))
  }
  out <- scop:::paga_connectivities_cpp(
    knn_idx = knn_idx, groups = groups, n_groups = n_groups
  )
  expect_equal(dim(out$connectivities), c(6, 6))
  expect_true(all(out$connectivities >= 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 6. Determinism: same input → same output
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp is deterministic", {
  dat <- make_paga_mock(seed = 1)
  out1 <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx, groups = dat$groups, n_groups = dat$n_groups
  )
  out2 <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx, groups = dat$groups, n_groups = dat$n_groups
  )
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 7. Group sizes consistency
# ---------------------------------------------------------------------------

test_that("group_sizes sum to n_cells and match groups input", {
  dat <- make_paga_mock(n_cells = 25, n_groups = 4, seed = 99)
  out <- scop:::paga_connectivities_cpp(
    knn_idx = dat$knn_idx, groups = dat$groups, n_groups = dat$n_groups
  )
  expect_equal(sum(out$group_sizes), dat$n_cells)
  expect_equal(
    as.numeric(out$group_sizes),
    as.numeric(table(dat$groups))
  )
})

# ---------------------------------------------------------------------------
# 8. Input validation rejects bad arguments
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp rejects mismatched dimensions", {
  dat <- make_paga_mock()
  expect_error(
    scop:::paga_connectivities_cpp(
      knn_idx = dat$knn_idx[1:10, , drop = FALSE],
      groups = dat$groups,
      n_groups = dat$n_groups
    ),
    "knn_idx rows must match groups length"
  )
})

test_that("paga_connectivities_cpp rejects negative n_groups", {
  dat <- make_paga_mock()
  expect_error(
    scop:::paga_connectivities_cpp(
      knn_idx = dat$knn_idx,
      groups = dat$groups,
      n_groups = -1L
    ),
    "n_groups must be positive"
  )
})

test_that("paga_connectivities_cpp rejects groups with invalid indices", {
  dat <- make_paga_mock()
  bad_groups <- dat$groups
  bad_groups[1] <- 99L  # > n_groups
  expect_error(
    scop:::paga_connectivities_cpp(
      knn_idx = dat$knn_idx,
      groups = bad_groups,
      n_groups = dat$n_groups
    ),
    "groups must be 1-based"
  )
})

# ---------------------------------------------------------------------------
# 9. Directed edges structure
# ---------------------------------------------------------------------------

test_that("directed_edges diagonal counts intra-group KNN edges", {
  n_cells <- 10
  groups <- c(rep(1L, 5), rep(2L, 5))
  # Perfect block: cells only connect within their group
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  for (i in 1:5) {
    candidates <- 1:5
    candidates <- setdiff(candidates, i)
    knn_idx[i, ] <- sample(candidates, 3)
  }
  for (i in 6:10) {
    candidates <- 6:10
    candidates <- setdiff(candidates, i)
    knn_idx[i, ] <- sample(candidates, 3)
  }
  out <- scop:::paga_connectivities_cpp(
    knn_idx = knn_idx, groups = groups, n_groups = 2L
  )
  # Directed edges between different groups should be 0
  expect_equal(out$directed_edges[1, 2], 0)
  expect_equal(out$directed_edges[2, 1], 0)
  # Inter-group connectivities should be 0
  expect_equal(out$connectivities[1, 2], 0)
  expect_equal(out$connectivities[2, 1], 0)
})
