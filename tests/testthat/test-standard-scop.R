test_that("standard_scop defaults to 2D UMAP", {
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

  out <- suppressWarnings(standard_scop(
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

test_that("standard_scop runs explicit 3D UMAP through RunUMAP2 path", {
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

  out <- suppressWarnings(standard_scop(
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

test_that("standard_scop records a compact FindClusters command", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("uwot")

  set.seed(11)
  counts <- Matrix::rsparsematrix(120, 50, density = 0.08)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  out <- suppressWarnings(standard_scop(
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
