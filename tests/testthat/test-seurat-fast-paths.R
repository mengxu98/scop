make_fast_path_object <- function(seed = 7) {
  set.seed(seed)
  mat <- Matrix::rsparsematrix(120, 80, density = 0.08)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 60, verbose = FALSE)
  obj <- ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
  RunPCA(obj, features = SeuratObject::VariableFeatures(obj), npcs = 10, verbose = FALSE)
}

test_that("FindNeighbors fast path stores Seurat-compatible graphs", {
  obj <- make_fast_path_object()
  out <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:10,
    k.param = 10,
    graph.name = c("RNA_nn", "RNA_snn"),
    verbose = FALSE
  )

  expect_true(all(c("RNA_nn", "RNA_snn") %in% names(out@graphs)))
  expect_s4_class(out@graphs$RNA_nn, "Graph")
  expect_s4_class(out@graphs$RNA_snn, "Graph")
  expect_equal(dim(out@graphs$RNA_snn), rep(ncol(out), 2))
  expect_equal(out@graphs$RNA_snn@assay.used, "RNA")
})

test_that("RunPCA fast path stores total variance and cell embeddings", {
  obj <- make_fast_path_object()
  expect_true("total.variance" %in% names(obj@reductions$pca@misc))
  expect_true(is.finite(obj@reductions$pca@misc$total.variance))
  expect_equal(nrow(SeuratObject::Embeddings(obj[["pca"]])), ncol(obj))
})
