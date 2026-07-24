make_monocle2_branching_srt <- function(n_genes = 90L, n_per_branch = 50L, seed = 1L) {
  set.seed(seed)
  n_cells <- n_per_branch * 3L
  branch <- rep(c("root", "branch_a", "branch_b"), each = n_per_branch)
  branch_time <- rep(seq(0, 1, length.out = n_per_branch), 3L)
  mu <- matrix(1, nrow = n_genes, ncol = n_cells)
  mu[1:20, branch == "root"] <- 5 + branch_time[branch == "root"] * 2
  mu[21:45, branch == "branch_a"] <- 3 + branch_time[branch == "branch_a"] * 8
  mu[46:70, branch == "branch_b"] <- 3 + branch_time[branch == "branch_b"] * 8
  mu[71:n_genes, ] <- 2
  counts <- matrix(
    rnbinom(n_genes * n_cells, mu = as.vector(mu), size = 2),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  srt <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$branch <- branch
  SeuratObject::VariableFeatures(srt) <- rownames(srt)
  srt
}

expect_monocle2_cpp_matches_r <- function(out_r, out_cpp, tolerance = 1e-8) {
  expect_equal(
    as.numeric(out_cpp$Monocle2_Pseudotime),
    as.numeric(out_r$Monocle2_Pseudotime),
    tolerance = tolerance
  )
  expect_equal(
    as.character(out_cpp$Monocle2_State),
    as.character(out_r$Monocle2_State)
  )
}

test_that("monocle2_order_from_mst_cpp follows Monocle2 DFS pseudotime semantics", {
  distances <- as.matrix(stats::dist(matrix(c(0, 1, 2, 1, 1), ncol = 1)))
  edges <- matrix(c(1L, 2L, 2L, 3L, 2L, 4L, 2L, 5L), ncol = 2L, byrow = TRUE)

  out1 <- scop:::monocle2_order_from_mst_cpp(distances, edges, root_cell = 1L)
  out2 <- scop:::monocle2_order_from_mst_cpp(distances, edges, root_cell = 1L)

  expect_equal(out1, out2)
  expect_equal(out1$root_cell, 1L)
  expect_equal(out1$pseudotime[1], 0)
  expect_equal(out1$pseudotime[3], 2)
  expect_true(length(unique(out1$state)) >= 2L)
})

test_that("RunMonocle2 cpp backend stores Monocle2-compatible results", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(1)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- rep(c("A", "B"), each = 60)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)

  out <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    group.by = "group",
    backend = "cpp",
    n_neighbors = 10,
    root_state = "A",
    verbose = FALSE
  )

  expect_true("DDRTree" %in% names(out@reductions))
  expect_true("Monocle2_State" %in% colnames(out@meta.data))
  expect_true("Monocle2_Pseudotime" %in% colnames(out@meta.data))
  expect_equal(out@tools$Monocle2$backend, "cpp")
  expect_equal(out@tools$Monocle2$parameters$reduction, "DDRTree")
  expect_s4_class(out@tools$Monocle2$cds, "CellDataSet")
  expect_s3_class(out@tools$Monocle2$graph, "data.frame")
  expect_true(all(out$Monocle2_Pseudotime >= 0))
  expect_true(all(is.finite(out$Monocle2_Pseudotime)))
  expect_equal(out@tools$Monocle2$parameters$engine, "monocle_mst_cpp_order")
})

test_that("RunMonocle2 cpp backend matches R ordering on the same root state", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(1)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)

  out_r <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    backend = "r",
    root_state = 1,
    verbose = FALSE
  )
  out_cpp <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    backend = "cpp",
    root_state = 1,
    verbose = FALSE
  )

  expect_equal(
    as.numeric(out_cpp$Monocle2_Pseudotime),
    as.numeric(out_r$Monocle2_Pseudotime),
    tolerance = 1e-8
  )
  expect_equal(
    as.character(out_cpp$Monocle2_State),
    as.character(out_r$Monocle2_State)
  )
})

test_that("RunMonocle2 cpp backend matches R ordering with DDRTree iteration override", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(1)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)

  out_r <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    backend = "r",
    root_state = 1,
    ddrtree_maxIter = 5,
    verbose = FALSE
  )
  out_cpp <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    backend = "cpp",
    root_state = 1,
    ddrtree_maxIter = 5,
    verbose = FALSE
  )

  expect_equal(out_cpp@tools$Monocle2$parameters$ddrtree_maxIter, 5)
  expect_equal(
    as.numeric(out_cpp$Monocle2_Pseudotime),
    as.numeric(out_r$Monocle2_Pseudotime),
    tolerance = 1e-8
  )
  expect_equal(
    as.character(out_cpp$Monocle2_State),
    as.character(out_r$Monocle2_State)
  )
})

test_that("RunMonocle2 cpp backend preserves non-DDRTree Monocle2 failures", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(1)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)

  expect_error(
    RunMonocle2(
      srt,
      features = rownames(srt)[1:30],
      backend = "r",
      reduction_method = "SimplePPT",
      root_state = 1,
      verbose = FALSE
    ),
    "unrecognized dimensionality reduction method"
  )
  expect_error(
    RunMonocle2(
      srt,
      features = rownames(srt)[1:30],
      backend = "cpp",
      reduction_method = "SimplePPT",
      root_state = 1,
      verbose = FALSE
    ),
    "unrecognized dimensionality reduction method"
  )
})

test_that("RunMonocle2 cpp backend supports Disp feature selection", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(2)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)

  out <- RunMonocle2(
    srt,
    feature_type = "Disp",
    disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 0.5 * dispersion_fit",
    backend = "cpp",
    n_neighbors = 10,
    verbose = FALSE
  )

  expect_true("DDRTree" %in% names(out@reductions))
  expect_equal(out@tools$Monocle2$parameters$reduction, "DDRTree")
  expect_s3_class(out@tools$Monocle2$dispersion_table, "data.frame")
  expect_true(length(out@tools$Monocle2$features) >= 2L)
})

test_that("RunMonocle2 cpp backend matches R with Disp feature selection", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(2)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  run_args <- list(
    feature_type = "Disp",
    disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 0.5 * dispersion_fit",
    root_state = 1,
    n_neighbors = 10,
    verbose = FALSE
  )

  out_r <- do.call(RunMonocle2, c(list(srt = srt, backend = "r"), run_args))
  out_cpp <- do.call(RunMonocle2, c(list(srt = srt, backend = "cpp"), run_args))

  expect_equal(out_cpp@tools$Monocle2$features, out_r@tools$Monocle2$features)
  expect_monocle2_cpp_matches_r(out_r, out_cpp)
})

test_that("RunMonocle2 cpp backend supports root_state by trajectory state", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  set.seed(1)
  counts <- matrix(
    rnbinom(50L * 120L, mu = 5, size = 1),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("c", 1:120))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)

  out_default <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    backend = "cpp",
    n_neighbors = 10,
    verbose = FALSE
  )
  state_levels <- sort(unique(as.character(out_default$Monocle2_State)))
  skip_if(length(state_levels) < 2L)

  out_root <- RunMonocle2(
    srt,
    features = rownames(srt)[1:30],
    backend = "cpp",
    root_state = state_levels[2L],
    n_neighbors = 10,
    verbose = FALSE
  )

  expect_true(out_root@tools$Monocle2$parameters$root_state %in% state_levels[2L])
  expect_false(identical(
    as.numeric(out_default$Monocle2_Pseudotime),
    as.numeric(out_root$Monocle2_Pseudotime)
  ))
})

test_that("RunMonocle2 cpp backend matches R ordering on a multi-state trajectory", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  srt <- make_monocle2_branching_srt()
  run_args <- list(
    features = rownames(srt)[1:80],
    root_state = 1,
    expressionFamily = "uninormal",
    norm_method = "none",
    ddrtree_maxIter = 5,
    verbose = FALSE
  )

  out_r <- do.call(RunMonocle2, c(list(srt = srt, backend = "r"), run_args))
  out_cpp <- do.call(RunMonocle2, c(list(srt = srt, backend = "cpp"), run_args))

  expect_gte(length(unique(as.character(out_r$Monocle2_State))), 2)
  expect_equal(length(unique(as.character(out_cpp$Monocle2_State))), length(unique(as.character(out_r$Monocle2_State))))
  expect_monocle2_cpp_matches_r(out_r, out_cpp)
})

test_that("RunMonocle2 cpp backend covers feature selection and root parameter variants", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  srt <- make_monocle2_branching_srt()
  common_args <- list(
    root_state = 1,
    expressionFamily = "uninormal",
    norm_method = "none",
    ddrtree_maxIter = 5,
    verbose = FALSE
  )

  explicit_out <- do.call(
    RunMonocle2,
    c(list(srt = srt, backend = "cpp", features = rownames(srt)[1:80]), common_args)
  )
  hvf_out <- do.call(
    RunMonocle2,
    c(list(srt = srt, backend = "cpp", features = NULL, feature_type = "HVF"), common_args)
  )
  group_root_out <- RunMonocle2(
    srt,
    features = rownames(srt)[1:80],
    backend = "cpp",
    group.by = "branch",
    root_state = "branch_a",
    expressionFamily = "uninormal",
    norm_method = "none",
    ddrtree_maxIter = 5,
    verbose = FALSE
  )

  expect_equal(length(explicit_out@tools$Monocle2$features), 80)
  expect_equal(hvf_out@tools$Monocle2$features, SeuratObject::VariableFeatures(srt))
  expect_true(group_root_out@tools$Monocle2$parameters$root_state %in% unique(as.character(group_root_out$Monocle2_State)))
  expect_true(all(is.finite(explicit_out$Monocle2_Pseudotime)))
  expect_true(all(is.finite(hvf_out$Monocle2_Pseudotime)))
  expect_true(all(is.finite(group_root_out$Monocle2_Pseudotime)))
})

test_that("RunMonocle2 cpp backend produces plot-compatible tools and metadata", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  srt <- make_monocle2_branching_srt(n_per_branch = 30L)
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- RunMonocle2(
    srt,
    features = rownames(srt)[1:70],
    backend = "cpp",
    group.by = "branch",
    root_state = 1,
    expressionFamily = "uninormal",
    norm_method = "none",
    ddrtree_maxIter = 5,
    show_plot = TRUE,
    verbose = FALSE
  )

  expect_true("DDRTree" %in% names(out@reductions))
  expect_true(all(c("Monocle2_State", "Monocle2_Pseudotime") %in% colnames(out@meta.data)))
  expect_s4_class(out@tools$Monocle2$cds, "CellDataSet")
  expect_s3_class(out@tools$Monocle2$graph, "data.frame")
  expect_s3_class(out@tools$Monocle2$trajectory, "LayerInstance")
  expect_equal(out@tools$Monocle2$parameters$engine, "monocle_mst_cpp_order")
})

test_that("RunMonocle2 cpp backend preserves root_state error behavior", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  srt <- make_monocle2_branching_srt(n_per_branch = 30L)

  expect_error(
    RunMonocle2(
      srt,
      features = rownames(srt)[1:70],
      backend = "cpp",
      root_state = "missing_state",
      expressionFamily = "uninormal",
      norm_method = "none",
      ddrtree_maxIter = 5,
      verbose = FALSE
    ),
    "no cells for State"
  )
})
