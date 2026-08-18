make_python_backend_test_srt <- function(
  cells = c("cell1", "cell2"),
  genes = c("gene1", "gene2")
) {
  counts <- Matrix::sparseMatrix(
    i = c(1L, 2L),
    j = c(1L, 2L),
    x = c(1, 1),
    dims = c(length(genes), length(cells)),
    dimnames = list(genes, cells)
  )
  Seurat::CreateSeuratObject(counts = counts)
}

test_that("RunPAGA python backend promotes misc$paga into tools$PAGA", {
  adata_in <- list(obs = data.frame(cluster = factor(c("a", "b"))))
  paga_srt <- make_python_backend_test_srt()
  paga_srt$cluster <- factor(c("a", "b"))
  paga_srt@misc[["paga"]] <- list(
    connectivities = diag(2),
    connectivities_tree = diag(2),
    pos = matrix(c(0, 0, 1, 1), ncol = 2)
  )

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) NULL,
    py_to_r2 = function(x) x,
    palette_colors = function(...) c("a" = "#111111", "b" = "#222222"),
    scop_python_import = function(...) {
      list(PAGA = function(...) "fake_adata_out")
    },
    adata_to_srt = function(...) paga_srt
  )

  out <- RunPAGA(
    adata = adata_in,
    group.by = "cluster",
    backend = "python",
    return_seurat = TRUE,
    verbose = FALSE
  )

  expect_null(out@misc[["paga"]])
  expect_identical(out@tools$PAGA$backend, "python")
  expect_identical(out@tools$PAGA$groups, "cluster")
  expect_equal(out@tools$PAGA$connectivities, diag(2))
})

test_that("RunSCVELO python backend records result locations in tools", {
  adata_in <- list(obs = data.frame(cluster = factor(c("a", "b"))))
  scvelo_srt <- make_python_backend_test_srt()
  scvelo_srt$cluster <- factor(c("a", "b"))
  scvelo_srt$velocity_confidence <- c(0.1, 0.2)
  scvelo_srt$velocity_length <- c(0.3, 0.4)
  scvelo_srt$velocity_pseudotime <- c(0.5, 0.6)
  scvelo_srt[["velocity_umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      ncol = 2,
      dimnames = list(colnames(scvelo_srt), c("UMAP_1", "UMAP_2"))
    ),
    assay = "RNA",
    key = "velocityumap_"
  )
  scvelo_srt@misc[["velocity_graph"]] <- diag(2)

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) NULL,
    py_to_r2 = function(x) x,
    scop_python_import = function(...) {
      list(SCVELO = function(...) "fake_adata_out")
    },
    adata_to_srt = function(...) scvelo_srt
  )

  out <- RunSCVELO(
    adata = adata_in,
    group.by = "cluster",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    mode = "stochastic",
    backend = "python",
    return_seurat = TRUE,
    verbose = FALSE
  )

  expect_identical(out@tools$SCVELO$backend, "python")
  expect_identical(out@tools$SCVELO$group.by, "cluster")
  expect_equal(unlist(out@tools$SCVELO$mode, use.names = FALSE), "stochastic")
  expect_identical(
    out@tools$SCVELO$stochastic$velocity_reduction,
    "velocity_umap"
  )
  expect_identical(
    out@tools$SCVELO$stochastic$confidence_key,
    "velocity_confidence"
  )
  expect_identical(
    out@tools$SCVELO$stochastic$length_key,
    "velocity_length"
  )
  expect_identical(
    out@tools$SCVELO$stochastic$pseudotime_key,
    "velocity_pseudotime"
  )
  expect_equal(out@tools$SCVELO$stochastic$velocity_graph, diag(2))
})

test_that("RunPalantir python backend records converted outputs in tools", {
  adata_in <- list(obs = data.frame(cluster = factor(c("a", "b"))))
  palantir_srt <- make_python_backend_test_srt()
  palantir_srt$cluster <- factor(c("a", "b"))
  palantir_srt$palantir_pseudotime <- c(0.1, 0.9)
  palantir_srt$palantir_diff_potential <- c(0.3, 0.4)
  palantir_srt$branchA_diff_potential <- c(0.8, 0.2)
  palantir_srt[["palantir_dm"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      c(0, 0, 1, 1),
      ncol = 2,
      dimnames = list(colnames(palantir_srt), c("DM_1", "DM_2"))
    ),
    assay = "RNA",
    key = "palantirdm_"
  )
  palantir_srt@misc[["dm_kernel"]] <- diag(2)

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) NULL,
    check_python = function(...) TRUE,
    py_to_r2 = function(x) x,
    scop_python_import = function(...) {
      list(Palantir = function(...) "fake_adata_out")
    },
    adata_to_srt = function(...) palantir_srt
  )

  out <- RunPalantir(
    adata = adata_in,
    group.by = "cluster",
    linear_reduction = "pca",
    nonlinear_reduction = "umap",
    early_group = "a",
    backend = "python",
    return_seurat = TRUE,
    verbose = FALSE
  )

  expect_identical(out@tools$Palantir$backend, "python")
  expect_identical(out@tools$Palantir$group.by, "cluster")
  expect_equal(out@tools$Palantir$pseudotime, c(0.1, 0.9))
  expect_equal(out@tools$Palantir$diff_potential, c(0.3, 0.4))
  expect_equal(
    out@tools$Palantir$branch_probs,
    matrix(c(0.8, 0.2), ncol = 1, dimnames = list(colnames(palantir_srt), "branchA_diff_potential"))
  )
  expect_identical(
    out@tools$Palantir$diffusion_map_reduction,
    "palantir_dm"
  )
  expect_equal(out@tools$Palantir$dm_kernel, diag(2))
})
