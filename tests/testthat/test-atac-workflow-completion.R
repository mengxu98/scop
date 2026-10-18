make_atac_completion_object <- function() {
  counts <- withr::with_seed(17, Matrix::Matrix(matrix(
    stats::rpois(100 * 80, rep(seq(1, 3, length.out = 80), each = 100)), 100, 80,
    dimnames = list(paste0("p", 1:100), paste0("cell", 1:80))
  ), sparse = TRUE))
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 200, length.out = 100), width = 100))
  ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
}

test_that("ATAC representation records the actual dimensions and native depth correlations", {
  object <- make_atac_completion_object()
  result <- RunStandardWorkflow(object,
    prefix = "ATACFinal", assay = "peaks",
    nHVF = 80, linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nonlinear_reduction = character(), neighbor_k = 10, verbose = FALSE
  )
  reduction <- result[["ATACFinallsi"]]
  receipt <- SeuratObject::Misc(reduction, "scop_atac_representation")
  expected <- stats::cor(SeuratObject::Embeddings(reduction), result$nCount_peaks)
  expect_equal(unname(receipt$depth_correlation), as.vector(expected))
  expect_identical(receipt$dimensions, 2:5)
  expect_identical(receipt$cells, rownames(SeuratObject::Embeddings(reduction)))
  expect_identical(
    SeuratObject::LayerData(object, layer = "counts"),
    SeuratObject::LayerData(result, layer = "counts")
  )
  expect_identical(receipt, SeuratObject::Misc(result[["ATACFinalsvd"]], "scop_atac_representation"))
})

test_that("ATAC graph failures cannot return a successful partial workflow", {
  object <- make_atac_completion_object()
  testthat::local_mocked_bindings(FindNeighbors = function(...) stop("injected graph failure"), .package = "scop")
  expect_error(RunStandardWorkflow(object,
    assay = "peaks", nHVF = 80,
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nonlinear_reduction = character(), verbose = FALSE
  ), "[Ii]njected graph failure")
  expect_length(SeuratObject::Reductions(object), 0)
})

test_that("ATAC UMAP failures cannot return a successful partial workflow", {
  object <- make_atac_completion_object()
  testthat::local_mocked_bindings(RunUMAP = function(...) stop("injected UMAP failure"), .package = "scop")
  expect_error(RunStandardWorkflow(object,
    assay = "peaks", nHVF = 80,
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    neighbor_k = 10, verbose = FALSE
  ), "[Ii]njected UMAP failure")
  expect_length(SeuratObject::Reductions(object), 0)
})

test_that("ATAC clustering failures cannot return a successful partial workflow", {
  object <- make_atac_completion_object()
  testthat::local_mocked_bindings(FindClusters = function(...) stop("injected clustering failure"), .package = "scop")
  expect_error(RunStandardWorkflow(object,
    assay = "peaks", nHVF = 80,
    linear_reduction_dims = 6, linear_reduction_dims_use = 2:5,
    nonlinear_reduction = character(), neighbor_k = 10, verbose = FALSE
  ), "[Ii]njected clustering failure")
  expect_length(SeuratObject::Reductions(object), 0)
})
