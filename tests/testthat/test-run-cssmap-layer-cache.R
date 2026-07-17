make_cssmap_cache_objects <- function() {
  counts <- methods::as(Matrix::Matrix(
    matrix(
      c(5, 0, 2, 1, 3, 0, 4, 1, 2, 1, 3, 2),
      nrow = 4,
      dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  ), "dgCMatrix")
  query <- Seurat::NormalizeData(Seurat::CreateSeuratObject(counts), verbose = FALSE)
  reference <- Seurat::NormalizeData(Seurat::CreateSeuratObject(counts), verbose = FALSE)
  css <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(6), ncol = 2, dimnames = list(colnames(reference), c("CSS_1", "CSS_2"))),
    assay = "RNA",
    key = "CSS_"
  )
  css@misc$model <- list(sim2profiles = list())
  reference[["css"]] <- css
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(6), ncol = 2, dimnames = list(colnames(reference), c("UMAP_1", "UMAP_2"))),
    assay = "RNA",
    key = "UMAP_"
  )
  list(query = query, reference = reference)
}

test_that("RunCSSMap reuses the checked data layers for projection", {
  objects <- make_cssmap_cache_objects()
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
    CheckDataType = function(...) "lognormalized",
    invoke_fun = function(.fn, .args) {
      expect_identical(.args$object, query_data)
      matrix(seq_len(6), ncol = 2)
    },
    RunKNNMap = function(srt_query, ...) srt_query,
    .package = "scop"
  )

  out <- RunCSSMap(
    objects$query,
    objects$reference,
    ref_css = "css",
    ref_umap = "umap",
    projection_method = "knn",
    verbose = FALSE
  )

  expect_identical(calls, c("query:data", "reference:data"))
  expect_true("css_proj" %in% SeuratObject::Reductions(out))
})
