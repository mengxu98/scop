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
