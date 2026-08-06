make_consistency_object <- function(seed = 21, n_genes = 60, n_cells = 50) {
  set.seed(seed)
  mat <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.15)
  mat@x <- abs(mat@x * 10)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  obj <- Seurat::CreateSeuratObject(counts = mat)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 40, verbose = FALSE)
  obj <- ScaleData(obj, features = SeuratObject::VariableFeatures(obj), verbose = FALSE)
  obj
}

test_that("RunPCA default output is identical to Seurat", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)

  scop_pca <- RunPCA(
    obj,
    features = hvf,
    npcs = 8,
    reduction.name = "pca_scop",
    seed.use = 42,
    verbose = FALSE
  )
  seurat_pca <- Seurat::RunPCA(
    obj,
    features = hvf,
    npcs = 8,
    reduction.name = "pca_seurat",
    seed.use = 42,
    verbose = FALSE
  )

  expect_identical(
    SeuratObject::Embeddings(scop_pca[["pca_scop"]]),
    SeuratObject::Embeddings(seurat_pca[["pca_seurat"]])
  )
})

test_that("RunPCA cpp backend is opt-in and does not alter the default", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)
  expect_no_error(
    RunPCA(obj, features = hvf, npcs = 8, backend = "cpp", verbose = FALSE)
  )
})

test_that("FindNeighbors nn graph is identical to Seurat (self edges preserved)", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)
  obj <- RunPCA(obj, features = hvf, npcs = 8, verbose = FALSE)

  scop_nn <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_scop", "snn_scop"),
    verbose = FALSE
  )
  seurat_nn <- Seurat::FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_seurat", "snn_seurat"),
    verbose = FALSE
  )

  expect_identical(
    as.matrix(scop_nn@graphs[["nn_scop"]]),
    as.matrix(seurat_nn@graphs[["nn_seurat"]])
  )
  expect_equal(
    sum(Matrix::diag(scop_nn@graphs[["nn_scop"]])),
    ncol(obj),
    label = "nn graph keeps self edges like Seurat"
  )
})

test_that("FindNeighbors SNN graph is identical to Seurat", {
  obj <- make_consistency_object()
  hvf <- SeuratObject::VariableFeatures(obj)
  obj <- RunPCA(obj, features = hvf, npcs = 8, verbose = FALSE)

  scop_nn <- FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_scop", "snn_scop"),
    verbose = FALSE
  )
  seurat_nn <- Seurat::FindNeighbors(
    obj,
    reduction = "pca",
    dims = 1:8,
    k.param = 10,
    graph.name = c("nn_seurat", "snn_seurat"),
    verbose = FALSE
  )

  expect_identical(
    as.matrix(scop_nn@graphs[["snn_scop"]]),
    as.matrix(seurat_nn@graphs[["snn_seurat"]])
  )
})

test_that("VST on split Assay5 layers matches Seurat's layered reference", {
  obj <- make_consistency_object()
  obj[["RNA"]] <- split(
    obj[["RNA"]],
    f = rep(c("A", "B"), length.out = ncol(obj))
  )
  joined <- SeuratObject::JoinLayers(obj, assay = "RNA")

  layered_hvf <- SeuratObject::VariableFeatures(
    FindVariableFeatures(obj, selection.method = "vst", nfeatures = 40, verbose = FALSE)
  )
  seurat_layered_hvf <- SeuratObject::VariableFeatures(
    Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 40, verbose = FALSE)
  )
  joined_hvf <- SeuratObject::VariableFeatures(
    FindVariableFeatures(joined, selection.method = "vst", nfeatures = 40, verbose = FALSE)
  )

  # scop follows Seurat's per-layer VST semantics: the layered HVF must be
  # byte-identical to the reference implementation on the same layered input.
  expect_identical(layered_hvf, seurat_layered_hvf)
  # Seurat's layered behaviour itself differs from the joined reference
  # (per-layer vst + consensus merge vs one-shot joined vst); scop mirrors
  # that behaviour, so the layered result may differ from the joined one.
  expect_identical(
    identical(layered_hvf, joined_hvf),
    identical(seurat_layered_hvf, joined_hvf)
  )
})

test_that("VST aligns Assay5 layers with partially overlapping features", {
  set.seed(42)
  counts_a <- Matrix::Matrix(
    matrix(stats::rpois(80 * 30, lambda = 2), nrow = 80),
    sparse = TRUE,
    dimnames = list(paste0("g", 1:80), paste0("a", 1:30))
  )
  counts_b <- Matrix::Matrix(
    matrix(stats::rpois(80 * 30, lambda = 2), nrow = 80),
    sparse = TRUE,
    dimnames = list(paste0("g", 41:120), paste0("b", 1:30))
  )
  obj <- suppressWarnings(suppressMessages(
    Seurat::CreateSeuratObject(counts = list(A = counts_a, B = counts_b))
  ))

  scop_result <- suppressWarnings(FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  seurat_result <- suppressWarnings(Seurat::FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))

  expect_identical(
    SeuratObject::VariableFeatures(scop_result),
    SeuratObject::VariableFeatures(seurat_result)
  )
  expect_identical(
    rownames(scop_result[[SeuratObject::DefaultAssay(scop_result)]]),
    rownames(seurat_result[[SeuratObject::DefaultAssay(seurat_result)]])
  )
  scop_assay <- scop_result[[SeuratObject::DefaultAssay(scop_result)]]
  seurat_assay <- seurat_result[[SeuratObject::DefaultAssay(seurat_result)]]
  scop_meta <- methods::slot(scop_assay, "meta.data")
  seurat_meta <- methods::slot(seurat_assay, "meta.data")
  hvf_columns <- grep("^vf_vst_", colnames(scop_meta), value = TRUE)
  expect_identical(rownames(scop_meta), rownames(seurat_meta))
  expect_identical(
    unname(scop_meta[, hvf_columns, drop = FALSE]),
    unname(seurat_meta[, hvf_columns, drop = FALSE])
  )
  expect_true(isTRUE(methods::validObject(scop_assay, test = TRUE)))

  reordered_obj <- suppressWarnings(suppressMessages(
    Seurat::CreateSeuratObject(
      counts = list(A = counts_a, B = counts_b[rev(seq_len(nrow(counts_b))), ])
    )
  ))
  reordered_scop <- suppressWarnings(FindVariableFeatures(
    reordered_obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  reordered_seurat <- suppressWarnings(Seurat::FindVariableFeatures(
    reordered_obj,
    selection.method = "vst",
    nfeatures = 20,
    verbose = FALSE
  ))
  expect_identical(
    SeuratObject::VariableFeatures(reordered_scop),
    SeuratObject::VariableFeatures(reordered_seurat)
  )
})
