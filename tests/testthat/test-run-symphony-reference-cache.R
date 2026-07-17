test_that("buildReferenceFromSeurat reads the variable-gene layer once", {
  counts <- Matrix::Matrix(
    matrix(
      c(2, 1, 3, 4, 5, 2, 1, 4, 3, 6, 2, 5),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  srt <- Seurat::CreateSeuratObject(counts)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)
  embeddings <- matrix(
    seq_len(6), ncol = 2,
    dimnames = list(colnames(srt), c("DIM_1", "DIM_2"))
  )
  srt[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    loadings = matrix(
      seq_len(8), ncol = 2,
      dimnames = list(rownames(srt), c("PC_1", "PC_2"))
    ),
    assay = "RNA",
    key = "PC_"
  )
  harmony <- SeuratObject::CreateDimReducObject(embeddings, assay = "RNA", key = "HARMONY_")
  harmony@misc$R <- matrix(1, nrow = ncol(srt), ncol = 2)
  srt[["harmony"]] <- harmony
  umap <- SeuratObject::CreateDimReducObject(embeddings, assay = "RNA", key = "UMAP_")
  umap@misc$model <- list(mock = TRUE)
  srt[["umap"]] <- umap
  expr <- matrix(
    c(2, 3, 5, 4, 6, 7, 8, 5, 4, 9, 3, 6),
    nrow = 4,
    dimnames = list(rownames(srt), colnames(srt))
  )
  get_calls <- 0L

  testthat::local_mocked_bindings(
    GetAssayData5 = function(...) {
      get_calls <<- get_calls + 1L
      expr
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "symphony")
      switch(name,
        rowSDs = function(A, row_means) apply(A, 1, stats::sd),
        cosine_normalize = function(V, dim) V,
        compute_ref_cache = function(...) list(matrix(1), matrix(1)),
        stop("unexpected symphony function: ", name)
      )
    },
    .package = "scop"
  )

  ref <- buildReferenceFromSeurat(srt, pca = "pca", harmony = "harmony", umap = "umap", verbose = FALSE)

  expect_identical(get_calls, 1L)
  expect_equal(unname(ref$vargenes$mean), unname(rowMeans(expr)), tolerance = 1e-12)
  expect_equal(unname(ref$vargenes$stddev), unname(apply(expr, 1, stats::sd)), tolerance = 1e-12)
})
