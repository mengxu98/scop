make_cpp_backend_safety_object <- function(
  n_genes = 40L,
  n_cells = 24L,
  seed = 20260731L
) {
  set.seed(seed)
  cells <- paste0("cell", seq_len(n_cells))
  genes <- paste0("gene", seq_len(n_genes))
  counts <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 3), nrow = n_genes),
    sparse = TRUE,
    dimnames = list(genes, cells)
  )
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$group <- factor(rep(c("g1", "g2", "g3"), length.out = n_cells))
  srt$col <- rep(seq_len(ceiling(sqrt(n_cells))), length.out = n_cells)
  srt$row <- rep(seq_len(ceiling(n_cells / ceiling(sqrt(n_cells)))), each = ceiling(sqrt(n_cells)))[
    seq_len(n_cells)
  ]

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

  spliced <- counts
  unspliced <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 2), nrow = n_genes),
    sparse = TRUE,
    dimnames = dimnames(counts)
  )
  srt[["spliced"]] <- SeuratObject::CreateAssay5Object(counts = spliced)
  srt[["unspliced"]] <- SeuratObject::CreateAssay5Object(counts = unspliced)

  from <- seq_len(n_cells)
  to <- c(seq.int(2L, n_cells), 1L)
  graph <- Matrix::sparseMatrix(
    i = c(from, to),
    j = c(to, from),
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(cells, cells)
  )
  srt@graphs[["connectivities"]] <- SeuratObject::as.Graph(graph)
  srt
}

test_that("dense-memory budget helper is deterministic and rejects unsafe work", {
  dense_gib <- getFromNamespace("cpp_dense_gib", "scop")
  dense_budget <- getFromNamespace("assert_cpp_dense_budget", "scop")

  expect_equal(dense_gib(1000, 1000, 1), 8e6 / 1024^3)
  expect_equal(
    dense_budget(100, 100, 2, Inf, "test"),
    dense_gib(100, 100, 2)
  )
  expect_error(
    dense_budget(1000, 1000, 2, 0.001, "test backend"),
    "exceeding max_dense_gib"
  )
  expect_error(
    dense_budget(1, 1, 1, 0, "test backend"),
    "positive number"
  )
})

test_that("approximate public backends require explicit opt-in", {
  srt <- make_cpp_backend_safety_object()

  expect_identical(eval(formals(RunCellRank)$backend)[[1]], "python")
  expect_identical(eval(formals(RunPalantir)$backend)[[1]], "python")
  expect_identical(eval(formals(RunSCENICPlus)$backend)[[1]], "python")

  expect_error(
    RunCellRank(srt = srt, backend = "cpp", verbose = FALSE),
    "allow_approximate = TRUE"
  )
  expect_error(
    RunPalantir(
      srt = srt,
      backend = "cpp",
      linear_reduction = "pca",
      early_cell = colnames(srt)[[1L]],
      verbose = FALSE
    ),
    "allow_approximate = TRUE"
  )
  expect_error(
    RunSCENICPlus(srt = srt, backend = "cpp", verbose = FALSE),
    "allow_approximate = TRUE"
  )
})

test_that("unsupported Palantir and PAGA C++ arguments fail before computation", {
  srt <- make_cpp_backend_safety_object()

  expect_error(
    RunPalantir(
      srt = srt,
      backend = "cpp",
      allow_approximate = TRUE,
      linear_reduction = "pca",
      nonlinear_reduction = "umap",
      early_cell = colnames(srt)[[1L]],
      adjust_early_cell = TRUE,
      verbose = FALSE
    ),
    "adjust_early_cell"
  )
  expect_error(
    RunPAGA(
      srt = srt,
      group.by = "group",
      linear_reduction = "pca",
      nonlinear_reduction = "umap",
      backend = "cpp",
      embedded_with_PAGA = TRUE,
      verbose = FALSE
    ),
    "embedded_with_PAGA"
  )
})

test_that("CellRank blocks an unsafe dense cell-by-cell allocation", {
  srt <- make_cpp_backend_safety_object(n_cells = 30L)
  expect_error(
    RunCellRank(
      srt = srt,
      backend = "cpp",
      allow_approximate = TRUE,
      kernel_type = "cytotrace",
      max_dense_gib = 1e-8,
      show_plot = FALSE,
      verbose = FALSE
    ),
    "exceeding max_dense_gib"
  )
})

test_that("CellRank and scVelo reject ignored plotting or analysis options", {
  srt <- make_cpp_backend_safety_object()

  expect_error(
    RunCellRank(
      srt = srt,
      backend = "cpp",
      allow_approximate = TRUE,
      show_plot = TRUE,
      verbose = FALSE
    ),
    "show_plot"
  )
  expect_error(
    RunSCVELO(
      srt = srt,
      assay_y = c("spliced", "unspliced"),
      backend = "cpp",
      denoise = TRUE,
      show_plot = FALSE,
      verbose = FALSE
    ),
    "denoise"
  )
})

test_that("scVelo blocks dense expression work before matrix coercion", {
  srt <- make_cpp_backend_safety_object()
  expect_error(
    RunSCVELO(
      srt = srt,
      assay_y = c("spliced", "unspliced"),
      linear_reduction = "pca",
      nonlinear_reduction = "umap",
      backend = "cpp",
      filter_genes = FALSE,
      max_dense_gib = 1e-8,
      show_plot = FALSE,
      verbose = FALSE
    ),
    "exceeding max_dense_gib"
  )
})

test_that("CytoSPACE blocks dense inputs before expression coercion", {
  spatial <- make_cpp_backend_safety_object(n_cells = 12L)
  reference <- make_cpp_backend_safety_object(n_cells = 18L, seed = 20260732L)
  reference$cell_type <- factor(rep(c("a", "b"), length.out = ncol(reference)))

  expect_error(
    RunCytoSPACE(
      spatial,
      reference = reference,
      reference_label = "cell_type",
      features = rownames(spatial),
      mean_cell_numbers = 1,
      max_dense_gib = 1e-8,
      verbose = FALSE
    ),
    "exceeding max_dense_gib"
  )
})

test_that("CytoSPACE public workflow records fraction and assignment backends", {
  spatial <- make_cpp_backend_safety_object(n_cells = 9L)
  reference <- make_cpp_backend_safety_object(n_cells = 14L, seed = 20260733L)
  reference$cell_type <- factor(rep(c("a", "b"), length.out = ncol(reference)))

  out <- RunCytoSPACE(
    spatial,
    reference = reference,
    reference_label = "cell_type",
    features = rownames(spatial)[1:20],
    cell_fractions = c(a = 0.5, b = 0.5),
    n_cells_per_spot = rep(1L, ncol(spatial)),
    backend = "cpp",
    seed = 19L,
    verbose = FALSE
  )

  parameters <- out@tools$CytoSPACE$parameters
  expect_identical(parameters$fraction_backend, "cpp")
  expect_identical(parameters$assignment_backend, "cpp")
  expect_gt(parameters$dense_working_set_gib_lower_bound, 0)
  expect_equal(nrow(out@tools$CytoSPACE$assigned_locations), ncol(spatial))
  expect_equal(
    unname(rowSums(out@tools$CytoSPACE$fractional_abundances_by_spot)),
    rep(1, ncol(spatial))
  )
})

test_that("PAGA public C++ path records its limited implementation scope", {
  skip_if_not_installed("BiocNeighbors")

  srt <- make_cpp_backend_safety_object()
  out <- RunPAGA(
    srt = srt,
    group.by = "group",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    n_neighbors = 5L,
    backend = "cpp",
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_identical(out@tools$PAGA$backend, "cpp")
  expect_false(out@tools$PAGA$implementation$exact_reference)
  expect_match(out@tools$PAGA$implementation$scope, "connectivity graph")
})

test_that("scVelo public C++ path records memory and result semantics", {
  srt <- make_cpp_backend_safety_object()
  out <- RunSCVELO(
    srt = srt,
    assay_y = c("spliced", "unspliced"),
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    n_neighbors = 5L,
    mode = "stochastic",
    backend = "cpp",
    filter_genes = FALSE,
    compute_velocity_graph = FALSE,
    compute_terminal_states = FALSE,
    compute_pseudotime = FALSE,
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_identical(out@tools$SCVELO$backend, "cpp")
  expect_false(out@tools$SCVELO$implementation$exact_reference)
  expect_gt(
    out@tools$SCVELO$implementation$dense_working_set_gib_lower_bound,
    0
  )
})

test_that("CellRank public C++ path records approximation and memory scope", {
  skip_if_not_installed("BiocNeighbors")

  srt <- make_cpp_backend_safety_object()
  out <- RunCellRank(
    srt = srt,
    group.by = "group",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    kernel_type = "cytotrace",
    estimator_type = "CFLARE",
    n_neighbors = 5L,
    n_macrostates = 3L,
    backend = "cpp",
    allow_approximate = TRUE,
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_identical(out@tools$CellRank$backend, "cpp")
  expect_identical(out@tools$CellRank$estimator, "cflare_approximation")
  expect_false(out@tools$CellRank$implementation$exact_reference)
  expect_gt(
    out@tools$CellRank$implementation$dense_working_set_gib_lower_bound,
    0
  )
})

test_that("CellRank C++ pseudotime path preserves direction before connectivity mixing", {
  skip_if_not_installed("BiocNeighbors")

  srt <- make_cpp_backend_safety_object()
  srt$dpt_pseudotime <- seq(0, 1, length.out = ncol(srt))
  out <- RunCellRank(
    srt = srt,
    group.by = "group",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    kernel_type = "pseudotime",
    n_neighbors = 5L,
    backend = "cpp",
    allow_approximate = TRUE,
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_identical(
    out@tools$CellRank$kernel,
    "pseudotime_connectivity_combined"
  )
  expect_identical(out@tools$CellRank$parameters$n_macrostates, 3L)
})

test_that("CellRank hard pseudotime threshold retains a stochastic directed graph", {
  hard_threshold <- getFromNamespace("cellrank_hard_threshold_kernel", "scop")
  graph <- Matrix::sparseMatrix(
    i = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
    j = c(2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3),
    x = c(3, 2, 1, 3, 2, 1, 2, 2, 1, 1, 1, 1),
    dims = c(4, 4)
  )
  transition <- hard_threshold(
    graph,
    pseudotime = c(0, 0.25, 0.5, 1),
    frac_to_keep = 0.3
  )

  expect_equal(as.numeric(Matrix::rowSums(transition)), rep(1, 4), tolerance = 1e-12)
  expect_gt(transition[1, 2], 0)
  expect_equal(sum(transition[4, 1:3]), 0)
  expect_equal(transition[4, 4], 1)
})

test_that("CellRank connectivity kernel applies CellRank density normalization", {
  density_normalize <- getFromNamespace(
    "cellrank_density_normalize_connectivities",
    "scop"
  )
  graph <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 2, 3, 3, 4),
    j = c(2, 3, 1, 3, 4, 1, 2, 2),
    x = c(1, 2, 1, 3, 1, 2, 3, 1),
    dims = c(4, 4)
  )
  density <- Matrix::colSums(graph)
  inverse <- Matrix::Diagonal(x = 1 / density)
  expected <- inverse %*% graph %*% inverse
  expected <- Matrix::Diagonal(x = 1 / Matrix::rowSums(expected)) %*% expected
  observed <- density_normalize(graph)

  expect_equal(as.matrix(observed), as.matrix(expected), tolerance = 1e-12)
  expect_equal(as.numeric(Matrix::rowSums(observed)), rep(1, 4), tolerance = 1e-12)
})

test_that("CellRank GPCCA solves cell-level absorption for every macrostate", {
  transition <- matrix(0, 12, 12)
  for (i in seq_len(12)) {
    transition[i, i] <- 0.2
    if (i > 1L) transition[i, i - 1L] <- transition[i, i - 1L] + 0.4
    if (i < 12L) transition[i, i + 1L] <- transition[i, i + 1L] + 0.4
    transition[i, ] <- transition[i, ] / sum(transition[i, ])
  }
  result <- cellrank_gpcca_cpp(
    T_ = transition,
    n_states = 3L,
    n_cells_terminal = 1L
  )

  expect_equal(dim(result$absorption_probabilities), c(12L, 3L))
  expect_equal(as.numeric(rowSums(result$absorption_probabilities)), rep(1, 12), tolerance = 1e-8)
  expect_equal(sum(result$terminal_states != 0L), 3L)
  expect_true(result$absorption_converged)
  expect_true(result$absorption_method %in% c("sparse_lu", "fixed_point"))
  expect_true(result$membership_method %in% c(
    "optimized_inner_simplex",
    "inner_simplex_optimization_limit",
    "absolute_schur_fallback"
  ))
  expect_equal(as.numeric(rowSums(result$chi)), rep(1, 12), tolerance = 1e-8)
  expect_true(all(result$chi >= 0))
  if (!identical(result$membership_method, "absolute_schur_fallback")) {
    expect_length(result$simplex_indices, 3L)
    expect_true(all(result$simplex_indices >= 1L & result$simplex_indices <= 12L))
    expect_lte(
      result$membership_objective_final,
      result$membership_objective_initial + 1e-12
    )
  }
  expect_gt(length(unique(result$fate_confidence)), 3L)
})

test_that("CIBERSORT diagnostics are buffered instead of discarded", {
  source_path <- testthat::test_path("..", "..", "src", "cibersort_libsvm.cpp")
  testthat::skip_if_not(
    file.exists(source_path),
    "C++ source tree is unavailable in installed-package tests"
  )
  source_text <- paste(
    readLines(
      source_path,
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl(
    "cibersort_libsvm_reprintf\\(const char \\*, \\.\\.\\.\\) \\{\\}",
    source_text
  ))
  expect_match(source_text, "cibersort_svm_flush_diagnostics")
  expect_match(source_text, "std::mutex")
})
