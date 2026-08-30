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
    RunCellRank(
      srt = srt,
      backend = "cpp",
      show_plot = FALSE,
      verbose = FALSE
    ),
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

test_that("Palantir C++ supports dm_n_eigs", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_cpp_backend_safety_object()
  out <- RunPalantir(
    srt = srt,
    backend = "cpp",
    allow_approximate = TRUE,
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    n_neighbors = 5L,
    early_cell = colnames(srt)[[1L]],
    dm_n_eigs = 4L,
    verbose = FALSE
  )

  expect_identical(out@tools$Palantir$parameters$dm_n_eigs, 4L)
})

test_that("PAGA C++ keeps embedded_with_PAGA native", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_cpp_backend_safety_object()
  out <- RunPAGA(
    srt, group.by = "group", linear_reduction = "pca",
    nonlinear_reduction = "umap", n_neighbors = 5L,
    embedded_with_PAGA = TRUE, backend = "cpp",
    show_plot = FALSE, verbose = FALSE
  )
  expect_identical(out@tools$PAGA$backend, "cpp")
  expect_true("paga" %in% names(out@reductions))
  expect_true(all(is.finite(SeuratObject::Embeddings(out[["paga"]]))))
})

test_that("PHATE C++ handles gamma, native distance, MDS, and clustering", {
  set.seed(1)
  data <- matrix(stats::rnorm(180), nrow = 30L)
  rownames(data) <- paste0("cell", seq_len(nrow(data)))
  out <- RunPHATE(
    data, assay = "RNA", backend = "cpp", n_pca = NULL, n_landmark = 30L,
    knn = 5L, t = 3L, gamma = 0, knn_dist = "cosine",
    mds = "metric", do_cluster = TRUE, n_clusters = 3L
  )
  expect_identical(SeuratObject::Misc(out, "backend"), "cpp")
  expect_true(all(is.finite(SeuratObject::Embeddings(out))))
  expect_identical(nlevels(SeuratObject::Misc(out, "clusters")), 3L)
})

test_that("scVelo C++ handles native postprocessing parameters", {
  srt <- make_cpp_backend_safety_object()
  out <- RunSCVELO(
    srt, assay_y = c("spliced", "unspliced"), group.by = "group",
    linear_reduction = "pca", nonlinear_reduction = "umap",
    n_neighbors = 5L, backend = "cpp", filter_genes = FALSE,
    magic_impute = TRUE, knn = 3L, t = 2L,
    diff_kinetics = TRUE, denoise = TRUE, denoise_topn = 3L,
    calculate_velocity_genes = TRUE, show_plot = FALSE, verbose = FALSE
  )
  expect_identical(out@tools$SCVELO$backend, "cpp")
  expect_true(!is.null(out@tools$SCVELO$magic_imputed_data))
  expect_true(!is.null(out@tools$SCVELO$stochastic$differential_kinetics))
  expect_lte(out@tools$SCVELO$stochastic$n_velocity_graph_genes, 3L)
})

test_that("CellRank C++ accepts Schur, terminal, lineage, and graph controls", {
  skip_if_not_installed("BiocNeighbors")
  srt <- make_cpp_backend_safety_object()
  srt$custom_time <- seq(0, 1, length.out = ncol(srt))
  out <- RunCellRank(
    srt, group.by = "group", linear_reduction = "pca",
    nonlinear_reduction = "umap", kernel_type = "pseudotime",
    time_key = "custom_time", n_neighbors = 5L, n_macrostates = 3L,
    schur_method = "krylov", schur_n_components = 6L,
    terminal_states = "g3", terminal_state_agg = "union",
    driver_lineages = c(1L, 2L), recompute_neighbors = FALSE,
    backend = "cpp", allow_approximate = TRUE,
    show_plot = FALSE, verbose = FALSE
  )
  expect_identical(out@tools$CellRank$backend, "cpp")
  expect_identical(out@tools$CellRank$parameters$schur_n_components, 6L)
  expect_false(out@tools$CellRank$parameters$recompute_neighbors)
  expect_identical(ncol(out@tools$CellRank$lineage_drivers), 4L)
  expect_identical(
    sum(startsWith(colnames(out@meta.data), "cellrank_fate_Lineage")),
    ncol(out@tools$CellRank$fate_probabilities)
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
  expect_true(out@tools$SCVELO$parameters$normalize_per_cell)
})

test_that("scVelo C++ honors disabled per-cell normalization", {
  srt <- make_cpp_backend_safety_object()
  out <- RunSCVELO(
    srt = srt,
    assay_y = c("spliced", "unspliced"),
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    n_neighbors = 5L,
    mode = "deterministic",
    backend = "cpp",
    filter_genes = FALSE,
    normalize_per_cell = FALSE,
    compute_velocity_graph = FALSE,
    compute_terminal_states = FALSE,
    compute_pseudotime = FALSE,
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_false(out@tools$SCVELO$parameters$normalize_per_cell)
  expect_true(all(is.finite(Seurat::Embeddings(out[["deterministic_umap"]]))))
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

test_that("CellRank C++ honors a custom pseudotime key and computes drivers", {
  skip_if_not_installed("BiocNeighbors")

  srt <- make_cpp_backend_safety_object()
  srt$custom_time <- seq(0, 1, length.out = ncol(srt))
  out <- RunCellRank(
    srt = srt,
    group.by = "group",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    kernel_type = "pseudotime",
    time_key = "custom_time",
    n_neighbors = 5L,
    backend = "cpp",
    allow_approximate = TRUE,
    show_plot = FALSE,
    verbose = FALSE
  )

  expect_identical(out@tools$CellRank$parameters$pseudotime_source, "custom_time")
  expect_identical(out@tools$CellRank$parameters$time_key, "custom_time")
  expect_true(any(grepl("_corr$", colnames(out@tools$CellRank$lineage_drivers))))
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

test_that("CellRank rebuilds Scanpy-compatible UMAP fuzzy connectivities", {
  build_connectivities <- getFromNamespace("cellrank_umap_connectivities", "scop")
  knn_idx <- matrix(c(2L, 3L, 1L, 4L, 4L, 1L, 3L, 2L), nrow = 4, byrow = TRUE)
  knn_dist <- matrix(c(0.2, 0.8, 0.2, 0.8, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE)
  graph <- as.matrix(build_connectivities(knn_idx, knn_dist))
  far_membership <- log2(3) - 1

  expect_equal(graph, t(graph), tolerance = 1e-12)
  expect_equal(diag(graph), rep(0, 4), tolerance = 1e-12)
  expect_equal(graph[cbind(1:4, c(2, 1, 4, 3))], rep(1, 4), tolerance = 1e-12)
  expect_equal(graph[1, 3], 2 * far_membership - far_membership^2, tolerance = 1e-5)
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
