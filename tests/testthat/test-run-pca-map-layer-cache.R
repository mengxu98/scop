make_pcamap_cache_objects <- function() {
  counts <- Matrix::Matrix(
    matrix(
      c(4, 1, 2, 3, 1, 5, 2, 4, 3, 6, 5, 2),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  query <- Seurat::CreateSeuratObject(counts)
  reference <- Seurat::CreateSeuratObject(counts)
  loadings <- matrix(
    c(0.4, -0.2, 0.1, 0.3, 0.2, 0.5, -0.4, 0.1),
    nrow = 4,
    dimnames = list(rownames(reference), c("PC_1", "PC_2"))
  )
  reference[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(6), ncol = 2,
      dimnames = list(colnames(reference), c("PC_1", "PC_2"))
    ),
    loadings = loadings,
    assay = "RNA",
    key = "PC_"
  )
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(6), ncol = 2,
      dimnames = list(colnames(reference), c("UMAP_1", "UMAP_2"))
    ),
    assay = "RNA",
    key = "UMAP_"
  )
  list(query = query, reference = reference, loadings = loadings)
}

test_that("RunPCAMap reuses checked data layers for PCA projection", {
  objects <- make_pcamap_cache_objects()
  query_data <- matrix(
    c(1, 2, 3, 2, 4, 5, 3, 5, 7, 4, 6, 8),
    nrow = 4,
    dimnames = list(paste0("Gene", 1:4), colnames(objects$query))
  )
  ref_data <- matrix(
    c(2, 3, 5, 4, 6, 7, 8, 5, 4, 9, 3, 6),
    nrow = 4,
    dimnames = list(paste0("Gene", 1:4), colnames(objects$reference))
  )
  calls <- character()

  testthat::local_mocked_bindings(
    GetAssayData5 = function(object, layer, assay, ...) {
      calls <<- c(calls, paste(if (identical(object, objects$query)) "query" else "reference", layer, sep = ":"))
      if (identical(object, objects$query)) query_data else ref_data
    },
    CheckDataType = function(object, ...) {
      force(object)
      "lognormalized"
    },
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )

  out <- RunPCAMap(
    objects$query,
    objects$reference,
    ref_pca = "pca",
    ref_umap = "umap",
    ref_dims = 1:2,
    projection_method = "knn",
    verbose = FALSE
  )
  expected <- scale(
    t(query_data),
    center = rowMeans(ref_data),
    scale = apply(ref_data, 1, stats::sd)
  ) %*% objects$loadings

  expect_identical(calls, c("query:data", "reference:data"))
  expect_equal(SeuratObject::Embeddings(out[["ref.pca"]]), expected, tolerance = 1e-12)
})
