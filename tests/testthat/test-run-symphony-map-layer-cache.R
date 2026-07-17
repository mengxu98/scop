make_symphony_cache_objects <- function() {
  counts <- methods::as(Matrix::Matrix(
    matrix(
      c(5, 0, 2, 1, 3, 0, 4, 1, 2, 1, 3, 2),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  ), "dgCMatrix")
  query <- Seurat::CreateSeuratObject(counts)
  reference <- Seurat::CreateSeuratObject(counts)
  embeddings <- matrix(
    seq_len(6),
    ncol = 2,
    dimnames = list(colnames(reference), c("DIM_1", "DIM_2"))
  )
  reference[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings, assay = "RNA", key = "PC_")
  harmony <- SeuratObject::CreateDimReducObject(embeddings, assay = "RNA", key = "HARMONY_")
  harmony@misc$reduction_dims <- 1:2
  reference[["harmony"]] <- harmony
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings, assay = "RNA", key = "UMAP_")
  list(query = query, reference = reference)
}

test_that("RunSymphonyMap reuses checked query data in mapQuery", {
  objects <- make_symphony_cache_objects()
  query_data <- Matrix::Matrix(matrix(seq_len(12), nrow = 4), sparse = TRUE)
  ref_data <- Matrix::Matrix(matrix(seq_len(12) + 10, nrow = 4), sparse = TRUE)
  rownames(query_data) <- rownames(ref_data) <- paste0("Gene", 1:4)
  colnames(query_data) <- colnames(objects$query)
  colnames(ref_data) <- colnames(objects$reference)
  calls <- character()

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    GetAssayData5 = function(object, layer, assay, ...) {
      calls <<- c(calls, paste(if (identical(object, objects$query)) "query" else "reference", layer, sep = ":"))
      if (identical(object, objects$query)) query_data else ref_data
    },
    CheckDataType = function(object, ...) {
      force(object)
      "lognormalized"
    },
    buildReferenceFromSeurat = function(...) list(loadings = matrix(1, 4, 2)),
    mapQuery = function(exp_query, ...) {
      expect_identical(exp_query, query_data)
      pca <- matrix(seq_len(6), nrow = 2)
      harmony <- matrix(seq_len(6) + 10, nrow = 2)
      rownames(pca) <- c("PC_1", "PC_2")
      rownames(harmony) <- c("harmony_1", "harmony_2")
      colnames(pca) <- colnames(harmony) <- colnames(objects$query)
      list(
        Z_pca_query = pca,
        Zq_corr = harmony,
        R_query = matrix(1, nrow = 2, ncol = 3)
      )
    },
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )

  out <- RunSymphonyMap(
    objects$query,
    objects$reference,
    ref_pca = "pca",
    ref_harmony = "harmony",
    ref_umap = "umap",
    projection_method = "knn",
    verbose = FALSE
  )

  expect_identical(calls, c("query:data", "reference:data"))
  expect_true(all(c("ref.pca", "ref.harmony") %in% SeuratObject::Reductions(out)))
  expect_equal(
    out[["ref.pca"]]@stdev,
    apply(matrix(seq_len(6), nrow = 2), 1, stats::sd),
    tolerance = 1e-12
  )
  expect_equal(
    out[["ref.harmony"]]@stdev,
    apply(matrix(seq_len(6) + 10, nrow = 2), 1, stats::sd),
    tolerance = 1e-12
  )
})
