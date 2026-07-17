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

test_that("NormalizeData uses the fast path for Assay5 objects", {
  set.seed(11)
  mat <- Matrix::rsparsematrix(20, 10, density = 0.2)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)

  expect_true(inherits(obj[["RNA"]], "StdAssay"))
  expect_false("data" %in% SeuratObject::Layers(obj[["RNA"]], search = NA))
  out <- NormalizeData(obj, verbose = FALSE)

  expect_true(inherits(out[["RNA"]], "StdAssay"))
  expect_true("data" %in% SeuratObject::Layers(out[["RNA"]], search = NA))
  expect_equal(
    as.matrix(GetAssayData5(out, assay = "RNA", layer = "data")),
    as.matrix(Seurat::NormalizeData(mat, verbose = FALSE)),
    tolerance = 1e-8
  )
})

test_that("GetAssayData5 bypasses JoinLayers for an exact single layer", {
  obj <- make_fast_path_object()
  expected <- SeuratObject::GetAssayData(obj[["RNA"]], layer = "data")
  testthat::local_mocked_bindings(
    JoinLayers = function(...) stop("JoinLayers should not run for a single exact layer"),
    .package = "SeuratObject"
  )

  actual <- GetAssayData5(obj, assay = "RNA", layer = "data")
  expect_identical(actual, expected)
})

test_that("GetAssayData5 still joins split Assay5 layers", {
  obj <- make_fast_path_object()
  obj[["RNA"]] <- split(obj[["RNA"]], f = rep(c("A", "B"), length.out = ncol(obj)))
  expected <- SeuratObject::GetAssayData(
    SeuratObject::JoinLayers(obj[["RNA"]]),
    layer = "data"
  )

  actual <- GetAssayData5(obj, assay = "RNA", layer = "data")
  expect_identical(actual, expected)
})

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

test_that("NormalizeData supports legacy Assay objects", {
  mat <- Matrix::rsparsematrix(20, 10, density = 0.2)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  suppressWarnings(
    obj[["RNA"]] <- SeuratObject::CreateAssayObject(counts = mat)
  )

  expect_s4_class(obj[["RNA"]], "Assay")
  expect_no_error(out <- NormalizeData(obj, verbose = FALSE))
  expect_s4_class(out[["RNA"]], "Assay")
  expect_gt(length(out[["RNA"]]@data@x), 0)
  expect_equal(
    as.matrix(out[["RNA"]]@data),
    as.matrix(Seurat::NormalizeData(mat, verbose = FALSE)),
    tolerance = 1e-8
  )
})

test_that("Seurat fast paths support legacy Assay preprocessing", {
  set.seed(12)
  mat <- Matrix::rsparsematrix(80, 40, density = 0.15)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  suppressWarnings(
    obj[["RNA"]] <- SeuratObject::CreateAssayObject(counts = mat)
  )
  SeuratObject::Idents(obj) <- rep(c("A", "B"), length.out = ncol(obj))

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
  obj <- ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(
    obj,
    features = SeuratObject::VariableFeatures(obj),
    npcs = 10,
    verbose = FALSE
  )
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:10, k.param = 10, verbose = FALSE)
  markers <- FindMarkers(
    obj,
    ident.1 = "A",
    ident.2 = "B",
    verbose = FALSE
  )
  all_markers <- suppressWarnings(FindAllMarkers(obj, verbose = FALSE))

  expect_s4_class(obj[["RNA"]], "Assay")
  expect_equal(nrow(obj[["RNA"]]@scale.data), length(SeuratObject::VariableFeatures(obj)))
  expect_true("pca" %in% SeuratObject::Reductions(obj))
  expect_true(all(c("RNA_nn", "RNA_snn") %in% names(obj@graphs)))
  expect_s3_class(markers, "data.frame")
  expect_s3_class(all_markers, "data.frame")
})
