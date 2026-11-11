test_that("RunUMAP native uwot path matches Seurat exactly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("uwot")

  set.seed(41)
  counts <- Matrix::rsparsematrix(80, 60, density = 0.12)
  counts@x <- abs(counts@x * 10)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  object <- Seurat::CreateSeuratObject(counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = 50, verbose = FALSE)
  object <- Seurat::ScaleData(
    object,
    features = SeuratObject::VariableFeatures(object),
    verbose = FALSE
  )
  object <- Seurat::RunPCA(
    object,
    features = SeuratObject::VariableFeatures(object),
    npcs = 10,
    verbose = FALSE
  )

  actual <- RunUMAP(
    object,
    dims = 1:10,
    n.epochs = 30L,
    reduction.name = "umap_scop",
    verbose = FALSE
  )
  expected <- suppressWarnings(Seurat::RunUMAP(
    object,
    dims = 1:10,
    n.epochs = 30L,
    reduction.name = "umap_seurat",
    verbose = FALSE
  ))

  expect_identical(
    unname(SeuratObject::Embeddings(actual[["umap_scop"]])),
    unname(SeuratObject::Embeddings(expected[["umap_seurat"]]))
  )
})
