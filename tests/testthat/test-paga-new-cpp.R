# Tests for PAGA new C++ functions
#
# Covers:
#   1. paga_velocity_transitions_cpp
#   2. paga_root_cell_cpp
#   3. paga_diffusion_pseudotime_cpp (updated with n_branchings, min_group_size)

make_paga_mock <- function(n_cells = 20, n_groups = 3, n_neighbors = 5, seed = 42) {
  set.seed(seed)
  groups <- sample(seq_len(n_groups), n_cells, replace = TRUE)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, n_neighbors))
  }
  n_dims <- 2
  embedding <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)
  velocity_embedding <- matrix(rnorm(n_cells * n_dims, sd = 0.5), nrow = n_cells, ncol = n_dims)
  list(knn_idx = knn_idx, groups = groups, n_groups = n_groups,
       n_cells = n_cells, embedding = embedding,
       velocity_embedding = velocity_embedding)
}

# ---------------------------------------------------------------------------
# 1. Velocity transitions
# ---------------------------------------------------------------------------

test_that("paga_velocity_transitions_cpp returns valid structure", {
  dat <- make_paga_mock()
  out <- scop:::paga_velocity_transitions_cpp(
    velocity_embedding = dat$velocity_embedding,
    knn_idx = dat$knn_idx,
    groups = dat$groups,
    n_groups = dat$n_groups
  )
  expect_type(out, "list")
  expect_true("transitions_confidence" %in% names(out))
  expect_true("transitions_confidence_tree" %in% names(out))
  expect_true("group_sizes" %in% names(out))
  tc <- out$transitions_confidence
  expect_equal(dim(tc), c(dat$n_groups, dat$n_groups))
  expect_true(all(tc >= 0))
  # Rows should sum to ~1 (or 0 if no transitions)
  rs <- rowSums(tc)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-8))
})

test_that("paga velocity transitions is deterministic", {
  dat <- make_paga_mock(seed = 11)
  out1 <- scop:::paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups
  )
  out2 <- scop:::paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups
  )
  expect_equal(out1$transitions_confidence, out2$transitions_confidence)
})

test_that("velocity transitions with softmax_scale parameter", {
  dat <- make_paga_mock()
  out1 <- scop:::paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups, softmax_scale = 4.0
  )
  out2 <- scop:::paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups, softmax_scale = 1.0
  )
  expect_equal(dim(out1$transitions_confidence), dim(out2$transitions_confidence))
  # Different softmax scales should produce different matrices (on non-trivial data)
  expect_true(!identical(out1$transitions_confidence, out2$transitions_confidence))
})

test_that("group_sizes from velocity transitions match input", {
  dat <- make_paga_mock(n_cells = 30, n_groups = 4, seed = 5)
  out <- scop:::paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups
  )
  expect_length(out$group_sizes, dat$n_groups)
  expect_equal(sum(out$group_sizes), dat$n_cells, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2. Root cell
# ---------------------------------------------------------------------------

test_that("paga_root_cell_cpp returns valid cell indices", {
  dat <- make_paga_mock()
  root_cells <- scop:::paga_root_cell_cpp(
    embedding = dat$embedding,
    groups = dat$groups,
    root_group = 1L
  )
  expect_type(root_cells, "integer")
  expect_true(length(root_cells) >= 1)
  # All returned cells should be in group 1
  for (r in root_cells) {
    expect_true(r >= 1 && r <= dat$n_cells)
    expect_equal(dat$groups[r], 1)
  }
})

test_that("root cell first entry is in specified group", {
  n_cells <- 30; n_groups <- 3
  set.seed(7)
  groups <- sample(1:n_groups, n_cells, replace = TRUE)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  for (g in 1:n_groups) {
    if (any(groups == g)) {
      root_cells <- scop:::paga_root_cell_cpp(embedding, as.integer(groups), as.integer(g))
      expect_equal(groups[root_cells[1]], g)
    }
  }
})

test_that("root cell is deterministic for same input", {
  dat <- make_paga_mock(seed = 42)
  root1 <- scop:::paga_root_cell_cpp(dat$embedding, dat$groups, 1L)
  root2 <- scop:::paga_root_cell_cpp(dat$embedding, dat$groups, 1L)
  expect_equal(root1, root2)
})

# ---------------------------------------------------------------------------
# 3. Diffusion pseudotime (updated)
# ---------------------------------------------------------------------------

test_that("paga_diffusion_pseudotime_cpp returns valid structure", {
  n_groups <- 4
  set.seed(1)
  con <- matrix(runif(n_groups * n_groups), n_groups, n_groups)
  con <- (con + t(con)) / 2
  diag(con) <- 1

  out <- scop:::paga_diffusion_pseudotime_cpp(
    connectivities = con, root_group = 1L, n_dcs = 3L,
    n_branchings = 0L, min_group_size = 0.01
  )
  expect_type(out, "list")
  expect_true("pseudotime" %in% names(out))
  expect_true("diffusion_components" %in% names(out))
  expect_true("diffusion_eigenvalues" %in% names(out))
  expect_length(out$pseudotime, n_groups)
  # Pseudotime should be non-negative
  expect_true(all(out$pseudotime >= 0))
  # Pseudotime should be in [0, 1] (normalized)
  expect_true(max(out$pseudotime) <= 1 + 1e-10)
})

test_that("DPT root group has pseudotime 0", {
  n_groups <- 4
  set.seed(2)
  con <- matrix(runif(n_groups * n_groups), n_groups, n_groups)
  con <- (con + t(con)) / 2
  diag(con) <- 1

  out <- scop:::paga_diffusion_pseudotime_cpp(con, root_group = 1L, n_dcs = 3L)
  # Root group (group 1, index 1) should have pseudotime 0
  expect_equal(out$pseudotime[out$root_group], 0, tolerance = 1e-10)
})

test_that("DPT with n_branchings returns branch count", {
  n_groups <- 5
  set.seed(3)
  con <- diag(1, n_groups)
  con[1, 2] <- con[2, 1] <- 0.8
  con[2, 3] <- con[3, 2] <- 0.6
  con[3, 4] <- con[4, 3] <- 0.4
  for (i in 1:n_groups) for (j in 1:n_groups) {
    if (con[i, j] == 0 && i != j) con[i, j] <- 0.05
  }

  out <- scop:::paga_diffusion_pseudotime_cpp(con, 1L, n_dcs = 3L, n_branchings = 2L, min_group_size = 0.01)
  expect_true(out$n_branchings_found >= 0)
})