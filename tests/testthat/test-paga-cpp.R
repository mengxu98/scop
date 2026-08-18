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
#  10. paga_velocity_transitions_cpp
#  11. paga_root_cell_cpp
#  12. paga_diffusion_pseudotime_cpp (updated with n_branchings, min_group_size)
#  13. cell_dpt_pseudotime_cpp + RunPAGA dpt_pseudotime integration

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
  n_dims <- 2
  embedding <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)
  velocity_embedding <- matrix(
    rnorm(n_cells * n_dims, sd = 0.5),
    nrow = n_cells,
    ncol = n_dims
  )
  list(
    knn_idx = knn_idx,
    groups = groups,
    n_groups = n_groups,
    n_cells = n_cells,
    n_neighbors = n_neighbors,
    embedding = embedding,
    velocity_embedding = velocity_embedding
  )
}

# ---------------------------------------------------------------------------
# 1. Basic structural tests
# ---------------------------------------------------------------------------

test_that("paga_connectivities_cpp returns a named list with correct fields", {
  dat <- make_paga_mock()
  out <- paga_connectivities_cpp(
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
  out <- paga_connectivities_cpp(
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
  out <- paga_connectivities_cpp(
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
  out <- paga_connectivities_cpp(
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
  out <- paga_connectivities_cpp(
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
    paga_connectivities_cpp(knn_idx = knn_idx, groups = groups, n_groups = 1),
    NA # should not error
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
  out <- paga_connectivities_cpp(
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
  out <- paga_connectivities_cpp(
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
  out1 <- paga_connectivities_cpp(
    knn_idx = dat$knn_idx, groups = dat$groups, n_groups = dat$n_groups
  )
  out2 <- paga_connectivities_cpp(
    knn_idx = dat$knn_idx, groups = dat$groups, n_groups = dat$n_groups
  )
  expect_equal(out1, out2)
})

# ---------------------------------------------------------------------------
# 7. Group sizes consistency
# ---------------------------------------------------------------------------

test_that("group_sizes sum to n_cells and match groups input", {
  dat <- make_paga_mock(n_cells = 25, n_groups = 4, seed = 99)
  out <- paga_connectivities_cpp(
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
    paga_connectivities_cpp(
      knn_idx = dat$knn_idx[1:10, , drop = FALSE],
      groups = dat$groups,
      n_groups = dat$n_groups
    ),
    "knn_idx rows must match groups length", ignore.case = TRUE
  )
})

test_that("paga_connectivities_cpp rejects negative n_groups", {
  dat <- make_paga_mock()
  expect_error(
    paga_connectivities_cpp(
      knn_idx = dat$knn_idx,
      groups = dat$groups,
      n_groups = -1L
    ),
    "n_groups must be positive", ignore.case = TRUE
  )
})

test_that("paga_connectivities_cpp rejects groups with invalid indices", {
  dat <- make_paga_mock()
  bad_groups <- dat$groups
  bad_groups[1] <- 99L # > n_groups
  expect_error(
    paga_connectivities_cpp(
      knn_idx = dat$knn_idx,
      groups = bad_groups,
      n_groups = dat$n_groups
    ),
    "groups must be 1-based", ignore.case = TRUE
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
  out <- paga_connectivities_cpp(
    knn_idx = knn_idx, groups = groups, n_groups = 2L
  )
  # Directed edges between different groups should be 0
  expect_equal(out$directed_edges[1, 2], 0)
  expect_equal(out$directed_edges[2, 1], 0)
  # Inter-group connectivities should be 0
  expect_equal(out$connectivities[1, 2], 0)
  expect_equal(out$connectivities[2, 1], 0)
})

# ---------------------------------------------------------------------------
# 10. Velocity transitions
# ---------------------------------------------------------------------------

test_that("paga_velocity_transitions_cpp returns valid structure", {
  dat <- make_paga_mock()
  out <- paga_velocity_transitions_cpp(
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
  out1 <- paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups
  )
  out2 <- paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups
  )
  expect_equal(out1$transitions_confidence, out2$transitions_confidence)
})

test_that("velocity transitions with softmax_scale parameter", {
  dat <- make_paga_mock()
  out1 <- paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups,
    softmax_scale = 4.0
  )
  out2 <- paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups,
    softmax_scale = 1.0
  )
  expect_equal(dim(out1$transitions_confidence), dim(out2$transitions_confidence))
  # Different softmax scales should produce different matrices (on non-trivial data)
  expect_true(!identical(out1$transitions_confidence, out2$transitions_confidence))
})

test_that("group_sizes from velocity transitions match input", {
  dat <- make_paga_mock(n_cells = 30, n_groups = 4, seed = 5)
  out <- paga_velocity_transitions_cpp(
    dat$velocity_embedding, dat$knn_idx, dat$groups, dat$n_groups
  )
  expect_length(out$group_sizes, dat$n_groups)
  expect_equal(sum(out$group_sizes), dat$n_cells, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 11. Root cell
# ---------------------------------------------------------------------------

test_that("paga_root_cell_cpp returns valid cell indices", {
  dat <- make_paga_mock()
  root_cells <- paga_root_cell_cpp(
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
  n_cells <- 30
  n_groups <- 3
  set.seed(7)
  groups <- sample(1:n_groups, n_cells, replace = TRUE)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  for (g in 1:n_groups) {
    if (any(groups == g)) {
      root_cells <- paga_root_cell_cpp(embedding, as.integer(groups), as.integer(g))
      expect_equal(groups[root_cells[1]], g)
    }
  }
})

test_that("root cell is deterministic for same input", {
  dat <- make_paga_mock(seed = 42)
  root1 <- paga_root_cell_cpp(dat$embedding, dat$groups, 1L)
  root2 <- paga_root_cell_cpp(dat$embedding, dat$groups, 1L)
  expect_equal(root1, root2)
})

# ---------------------------------------------------------------------------
# 12. Diffusion pseudotime (updated)
# ---------------------------------------------------------------------------

test_that("paga_diffusion_pseudotime_cpp returns valid structure", {
  n_groups <- 4
  set.seed(1)
  con <- matrix(runif(n_groups * n_groups), n_groups, n_groups)
  con <- (con + t(con)) / 2
  diag(con) <- 1

  out <- paga_diffusion_pseudotime_cpp(
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

  out <- paga_diffusion_pseudotime_cpp(con, root_group = 1L, n_dcs = 3L)
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
  for (i in 1:n_groups) {
    for (j in 1:n_groups) {
      if (con[i, j] == 0 && i != j) con[i, j] <- 0.05
    }
  }

  out <- paga_diffusion_pseudotime_cpp(con, 1L, n_dcs = 3L, n_branchings = 2L, min_group_size = 0.01)
  expect_true(out$n_branchings_found >= 0)
})

# ---------------------------------------------------------------------------
# 13. Cell-level DPT + RunPAGA integration
# ---------------------------------------------------------------------------

test_that("cell_dpt_pseudotime_cpp returns per-cell pseudotime", {
  dat <- make_paga_mock(n_cells = 12, n_groups = 3, n_neighbors = 4, seed = 9)
  knn_dist <- matrix(NA_real_, dat$n_cells, ncol(dat$knn_idx))
  for (i in seq_len(dat$n_cells)) {
    for (k in seq_len(ncol(dat$knn_idx))) {
      j <- dat$knn_idx[i, k]
      knn_dist[i, k] <- sqrt(sum((dat$embedding[i, ] - dat$embedding[j, ])^2))
    }
  }
  root_cells <- paga_root_cell_cpp(
    embedding = dat$embedding[, 1:2, drop = FALSE],
    groups = dat$groups,
    root_group = 1L
  )
  out <- cell_dpt_pseudotime_cpp(
    knn_idx = dat$knn_idx,
    knn_dist = knn_dist,
    root_cell = root_cells[[1L]],
    n_dcs = 3L
  )
  expect_type(out, "list")
  expect_length(out$pseudotime, dat$n_cells)
  expect_true(all(out$pseudotime >= 0))
  expect_true(max(out$pseudotime) <= 1 + 1e-10)
  expect_equal(out$pseudotime[out$root_cell], 0, tolerance = 1e-10)
})

test_that("RunPAGA cpp infer_pseudotime writes dpt_pseudotime to meta.data", {
  skip_if_not_installed("BiocNeighbors")

  set.seed(1)
  n_cells <- 24L
  cells <- paste0("cell", seq_len(n_cells))
  counts <- Matrix::Matrix(
    matrix(stats::rpois(40 * n_cells, lambda = 3), nrow = 40),
    sparse = TRUE,
    dimnames = list(paste0("gene", seq_len(40)), cells)
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("g1", "g2", "g3"), length.out = n_cells))
  pca <- matrix(
    stats::rnorm(n_cells * 6L),
    nrow = n_cells,
    dimnames = list(cells, paste0("PC_", seq_len(6L)))
  )
  umap <- pca[, 1:2, drop = FALSE]
  colnames(umap) <- paste0("UMAP_", 1:2)
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca,
    assay = "RNA",
    key = "PC_"
  )
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap,
    assay = "RNA",
    key = "UMAP_"
  )

  out <- RunPAGA(
    srt = srt,
    group.by = "group",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    n_neighbors = 5L,
    backend = "cpp",
    infer_pseudotime = TRUE,
    root_group = "g1",
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_true("dpt_pseudotime" %in% colnames(out@meta.data))
  expect_length(out$dpt_pseudotime, ncol(out))
  expect_true(all(out$dpt_pseudotime >= 0))
  expect_true(max(out$dpt_pseudotime) <= 1 + 1e-10)
  expect_true(!is.null(out@tools$PAGA$dpt_pseudotime))
  expect_equal(out@tools$PAGA$dpt_pseudotime, out$dpt_pseudotime)
})

test_that("RunPAGA skips dense cell-level DPT on large inputs", {
  skip_if_not_installed("BiocNeighbors")

  n_cells <- 4100L
  cells <- paste0("cell", seq_len(n_cells))
  counts <- Matrix::Matrix(
    matrix(1, nrow = 5, ncol = n_cells),
    sparse = TRUE,
    dimnames = list(paste0("gene", seq_len(5)), cells)
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("g1", "g2"), length.out = n_cells))
  pca <- matrix(
    stats::rnorm(n_cells * 4L),
    nrow = n_cells,
    dimnames = list(cells, paste0("PC_", seq_len(4L)))
  )
  umap <- pca[, 1:2, drop = FALSE]
  colnames(umap) <- paste0("UMAP_", 1:2)
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca,
    assay = "RNA",
    key = "PC_"
  )
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap,
    assay = "RNA",
    key = "UMAP_"
  )

  expect_warning(
    out <- RunPAGA(
      srt = srt,
      group.by = "group",
      linear_reduction = "pca",
      nonlinear_reduction = "umap",
      n_neighbors = 5L,
      backend = "cpp",
      infer_pseudotime = TRUE,
      root_group = "g1",
      show_plot = FALSE,
      verbose = TRUE
    ),
    "skipped cell-level"
  )
  expect_false("dpt_pseudotime" %in% colnames(out@meta.data))
  expect_null(out@tools$PAGA$dpt_pseudotime)
  expect_false(is.null(out@tools$PAGA$pseudotime))
})
