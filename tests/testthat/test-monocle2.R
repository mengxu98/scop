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

test_that("RunMonocle2 stores Monocle2-compatible results", {
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
    root_state = 1,
    verbose = FALSE
  )

  expect_true("DDRTree" %in% names(out@reductions))
  expect_true("Monocle2_State" %in% colnames(out@meta.data))
  expect_true("Monocle2_Pseudotime" %in% colnames(out@meta.data))
  expect_s4_class(out@tools$Monocle2$cds, "CellDataSet")
  expect_true(inherits(out@tools$Monocle2$trajectory, "Layer"))
  expect_true(all(out$Monocle2_Pseudotime >= 0))
  expect_true(all(is.finite(out$Monocle2_Pseudotime)))
})

test_that("RunMonocle2 supports Disp feature selection", {
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
    verbose = FALSE
  )

  expect_true("DDRTree" %in% names(out@reductions))
  expect_true("Monocle2_Pseudotime" %in% colnames(out@meta.data))
  expect_s4_class(out@tools$Monocle2$cds, "CellDataSet")
  expect_true(all(is.finite(out$Monocle2_Pseudotime)))
})

test_that("RunMonocle2 reorders by root_state on a multi-state trajectory", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  srt <- make_monocle2_branching_srt()
  common_args <- list(
    features = rownames(srt)[1:80],
    expressionFamily = "uninormal",
    norm_method = "none",
    ddrtree_maxIter = 5,
    verbose = FALSE
  )

  out_default <- do.call(RunMonocle2, c(list(object = srt), common_args))
  expect_gte(length(unique(as.character(out_default$Monocle2_State))), 2)

  states <- sort(unique(as.character(out_default$Monocle2_State)))
  out_root <- do.call(
    RunMonocle2,
    c(list(object = srt, root_state = states[2]), common_args)
  )
  expect_false(identical(
    as.numeric(out_default$Monocle2_Pseudotime),
    as.numeric(out_root$Monocle2_Pseudotime)
  ))
})

test_that("RunMonocle2 preserves unknown reduction failures", {
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

  expect_error(
    RunMonocle2(
      srt,
      features = rownames(srt)[1:30],
      reduction_method = "SimplePPT",
      root_state = 1,
      verbose = FALSE
    ),
    "unrecognized dimensionality reduction method"
  )
})

test_that("RunMonocle2 deprecates backend and n_neighbors arguments", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("monocle")
  skip_if_not_installed("DDRTree")

  srt <- make_monocle2_branching_srt(n_per_branch = 30L)
  expect_warning(
    RunMonocle2(
      srt,
      features = rownames(srt)[1:70],
      backend = "cpp",
      n_neighbors = 10,
      expressionFamily = "uninormal",
      norm_method = "none",
      ddrtree_maxIter = 5,
      verbose = FALSE
    ),
    "deprecated"
  )
})
