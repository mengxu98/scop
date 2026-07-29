test_that("RunStandardWorkflow defaults to 2D UMAP", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("uwot")

  set.seed(1)
  counts <- Matrix::rsparsematrix(120, 50, density = 0.08)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)

  out <- suppressWarnings(RunStandardWorkflow(
    srt,
    nHVF = 40,
    linear_reduction_dims = 10,
    linear_reduction_dims_use = 1:5,
    neighbor_k = 10,
    verbose = FALSE,
    cores = 1
  ))

  expect_true("StandardpcaUMAP2D" %in% SeuratObject::Reductions(out))
  expect_false("StandardpcaUMAP3D" %in% SeuratObject::Reductions(out))
})

test_that("RunStandardWorkflow runs explicit 3D UMAP through RunUMAP2 path", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("uwot")

  set.seed(1)
  counts <- Matrix::rsparsematrix(120, 50, density = 0.08)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)

  out <- suppressWarnings(RunStandardWorkflow(
    srt,
    nHVF = 40,
    linear_reduction_dims = 10,
    linear_reduction_dims_use = 1:5,
    nonlinear_reduction_dims = c(2, 3),
    neighbor_k = 10,
    verbose = FALSE,
    cores = 1
  ))

  expect_true("Standardpca_SNN" %in% names(out@graphs))
  expect_true("StandardpcaUMAP2D" %in% SeuratObject::Reductions(out))
  expect_true("StandardpcaUMAP3D" %in% SeuratObject::Reductions(out))
})

test_that("RunStandardWorkflow records a compact FindClusters command", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("uwot")

  set.seed(11)
  counts <- Matrix::rsparsematrix(120, 50, density = 0.08)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  out <- suppressWarnings(RunStandardWorkflow(
    Seurat::CreateSeuratObject(counts),
    nHVF = 40,
    linear_reduction_dims = 10,
    linear_reduction_dims_use = 1:5,
    neighbor_k = 10,
    verbose = FALSE,
    cores = 1
  ))

  command <- out@commands[["FindClusters"]]
  expect_identical(command@params$resolution, 0.6)
  expect_lt(sum(nchar(command@call.string)), 1000)
})

test_that("legacy workflow aliases warn and forward until version 1.0.0", {
  testthat::local_mocked_bindings(
    RunStandardWorkflow = function(...) list(workflow = "standard", ...),
    RunIntegration = function(...) list(workflow = "integration", ...),
    .package = "scop"
  )

  expect_warning(
    standard <- standard_scop(value = 1),
    "deprecated.*RunStandardWorkflow.*removed in scop 1\\.0\\.0"
  )
  expect_identical(standard, list(workflow = "standard", value = 1))

  expect_warning(
    integration <- integration_scop(value = 2),
    "deprecated.*RunIntegration.*removed in scop 1\\.0\\.0"
  )
  expect_identical(integration, list(workflow = "integration", value = 2))
})
