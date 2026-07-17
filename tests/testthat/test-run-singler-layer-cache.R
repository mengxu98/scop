test_that("RunSingleR reuses checked logcounts when constructing SCE inputs", {
  counts <- Matrix::Matrix(
    matrix(
      c(2, 0, 1, 3, 1, 4),
      nrow = 2,
      dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
    ),
    sparse = TRUE
  )
  query <- Seurat::CreateSeuratObject(counts)
  reference <- Seurat::CreateSeuratObject(counts)
  reference$celltype <- c("A", "B", "A")
  query_data <- matrix(seq_len(6), nrow = 2, dimnames = dimnames(counts))
  ref_data <- matrix(seq_len(6) + 10, nrow = 2, dimnames = dimnames(counts))
  calls <- character()

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    GetAssayData5 = function(object, layer, assay, ...) {
      source <- if (identical(object, query)) "query" else "reference"
      calls <<- c(calls, paste(source, layer, sep = ":"))
      switch(paste(source, layer, sep = ":"),
        "query:data" = query_data,
        "reference:data" = ref_data,
        "query:counts" = counts,
        "reference:counts" = counts
      )
    },
    CheckDataType = function(object, ...) {
      force(object)
      "lognormalized"
    },
    .package = "scop"
  )
  testthat::local_mocked_bindings(
    SingleR = function(test, ref, ...) {
      expect_identical(SummarizedExperiment::assay(test, "logcounts"), query_data)
      expect_identical(SummarizedExperiment::assay(ref, "logcounts"), ref_data)
      list(
        labels = c("A", "B", "A"),
        pruned.labels = c("A", "B", "A"),
        scores = matrix(c(0.9, 0.1, 0.2, 0.8, 0.7, 0.3), nrow = 3,
          dimnames = list(NULL, c("A", "B")))
      )
    },
    .package = "SingleR"
  )

  out <- RunSingleR(query, reference, ref_group = "celltype", verbose = FALSE)

  expect_identical(calls, c("query:data", "reference:data", "query:counts", "reference:counts"))
  expect_identical(unname(out$singler_annotation), c("A", "B", "A"))
  expect_equal(unname(out$singler_score), c(0.9, 0.7, 0.2), tolerance = 1e-12)
})
