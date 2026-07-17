make_scmap_cache_objects <- function() {
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
  reference$celltype <- c("A", "B", "A")
  list(query = query, reference = reference)
}

test_that("RunScmap reuses checked logcounts when constructing SCE inputs", {
  objects <- make_scmap_cache_objects()
  query_data <- Matrix::Matrix(matrix(seq_len(12), nrow = 4), sparse = TRUE)
  ref_data <- Matrix::Matrix(matrix(seq_len(12) + 10, nrow = 4), sparse = TRUE)
  query_counts <- Matrix::Matrix(matrix(seq_len(12) + 20, nrow = 4), sparse = TRUE)
  ref_counts <- Matrix::Matrix(matrix(seq_len(12) + 30, nrow = 4), sparse = TRUE)
  for (x in list(query_data, ref_data, query_counts, ref_counts)) {
    rownames(x) <- paste0("Gene", 1:4)
  }
  colnames(query_data) <- colnames(query_counts) <- colnames(objects$query)
  colnames(ref_data) <- colnames(ref_counts) <- colnames(objects$reference)
  calls <- character()

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    GetAssayData5 = function(object, layer, assay, ...) {
      source <- if (identical(object, objects$query)) "query" else "reference"
      calls <<- c(calls, paste(source, layer, sep = ":"))
      switch(paste(source, layer, sep = ":"),
        "query:data" = query_data,
        "reference:data" = ref_data,
        "query:counts" = query_counts,
        "reference:counts" = ref_counts
      )
    },
    CheckDataType = function(...) "lognormalized",
    .package = "scop"
  )
  testthat::local_mocked_bindings(
    selectFeatures = function(x, ...) {
      SummarizedExperiment::rowData(x)$scmap_features <- TRUE
      x
    },
    indexCluster = function(x, ...) {
      S4Vectors::metadata(x)$scmap_cluster_index <- list(mock = TRUE)
      x
    },
    scmapCluster = function(projection, ...) {
      expect_identical(SummarizedExperiment::assay(projection, "logcounts"), query_data)
      list(
        scmap_cluster_labs = matrix(c("A", "B", "A"), ncol = 1),
        scmap_cluster_siml = matrix(c(0.9, 0.8, 0.7), ncol = 1)
      )
    },
    .package = "scmap"
  )

  out <- RunScmap(objects$query, objects$reference, ref_group = "celltype", nfeatures = 4, verbose = FALSE)

  expect_identical(calls, c("query:data", "reference:data", "query:counts", "reference:counts"))
  expect_identical(unname(out$scmap_annotation), c("A", "B", "A"))
})
