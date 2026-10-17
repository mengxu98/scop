test_that("RecoverCounts preserves expm1 and untransformed recovery with independent lookup", {
  counts <- matrix(c(1, 2, 3, 1, 3, 5, 1, 4, 2), nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3)))
  normalized <- sweep(counts, 2, colSums(counts), "/") * 10000
  testthat::local_mocked_bindings(CheckDataType = function(...) "log_normalized_counts", .package = "scop")
  for (trans in c("expm1", "none")) {
    values <- switch(trans, expm1 = log1p(normalized), none = normalized)
    object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(values, sparse = TRUE))
    actual <- RecoverCounts(object, trans = trans, sf = 10000, min_count = 1, verbose = FALSE)
    expect_equal(as.matrix(GetAssayData5(actual, layer = "counts")), counts)
    expect_equal(unname(actual$nCount_RNA), unname(colSums(counts)))
  }
  expect_error(RecoverCounts(object, trans = "invalid", verbose = FALSE), "arg")
  logged <- SeuratObject::CreateSeuratObject(Matrix::Matrix(log1p(normalized), sparse = TRUE))
  actual <- local({
    expm1 <- function(...) stop("caller transformation was used")
    RecoverCounts(logged, trans = "expm1", sf = 10000, min_count = 1, verbose = FALSE)
  })
  expect_equal(as.matrix(GetAssayData5(actual, layer = "counts")), counts)
})

test_that("Scanorama binding preserves embeddings and dense or sparse corrections", {
  make_object <- function(prefix) {
    counts <- matrix(seq_len(9), nrow = 3,
      dimnames = list(paste0("g", 1:3), paste0(prefix, 1:3)))
    Seurat::NormalizeData(SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE)), verbose = FALSE)
  }
  objects <- list(make_object("a"), make_object("b"))
  cells <- unlist(lapply(objects, colnames), use.names = FALSE)
  genes <- rownames(objects[[1]])
  embeddings <- matrix(seq_len(12), nrow = 6,
    dimnames = list(cells, paste0("Scanorama_", 1:2)))
  expression <- matrix(seq_len(18), nrow = 6, dimnames = list(cells, genes))
  sparse <- FALSE
  blocks <- function(value, use_sparse = FALSE) {
    if (use_sparse) value <- Matrix::Matrix(value, sparse = TRUE)
    stats::setNames(list(value[1:3, , drop = FALSE], value[4:6, , drop = FALSE]), c("first", "second"))
  }
  backend <- list(
    integrate = function(...) list(blocks(embeddings)),
    correct = function(...) list(blocks(embeddings), blocks(expression, sparse), genes)
  )
  testthat::local_mocked_bindings(import = function(...) backend, .package = "reticulate")
  testthat::local_mocked_bindings(
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) invisible(TRUE),
    python_cells_by_features = function(x) Matrix::t(x),
    CheckDataList = function(srt_list, ...) list(srt_list = srt_list, HVF = genes, assay = "RNA"),
    find_neighbors_and_clusters = function(srt, ...) srt,
    run_nonlinear_reduction = function(srt, ...) srt,
    .package = "scop"
  )
  for (sparse in c(FALSE, TRUE)) {
    for (corrected in c(FALSE, TRUE)) {
      actual <- suppressWarnings(Scanorama_integrate(srt_list = objects, HVF = genes, return_corrected = corrected, verbose = FALSE))
      actual_embeddings <- SeuratObject::Embeddings(actual, "Scanorama")
      expect_equal(unname(actual_embeddings), unname(embeddings))
      expect_identical(as.character(rownames(actual_embeddings)), cells)
      expect_identical(colnames(actual_embeddings), colnames(embeddings))
      expect_identical(SeuratObject::VariableFeatures(actual), genes)
      if (corrected) {
        actual_expression <- as.matrix(GetAssayData5(actual, assay = "Scanoramacorrected", layer = "data"))
        expect_equal(unname(actual_expression), unname(t(expression)))
        expect_identical(rownames(actual_expression), genes)
        expect_identical(as.character(colnames(actual_expression)), cells)
        expect_identical(SeuratObject::VariableFeatures(actual[["Scanoramacorrected"]]), genes)
      }
    }
  }
})
