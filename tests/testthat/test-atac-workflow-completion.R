make_atac_completion_object <- function() {
  skip_if_not_installed("Signac")
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

test_that("ATAC default dimensions fit a small measured object", {
  object <- subset(make_atac_completion_object(), cells = paste0("cell", 1:20))
  result <- suppressWarnings(RunStandardWorkflow(object,
    nHVF = 80, neighbor_k = 5, nonlinear_reduction = character(), verbose = FALSE
  ))
  receipt <- SeuratObject::Misc(result[["ATAClsi"]], "scop_atac_representation")
  expect_identical(receipt$dimensions, 2:19)
  expect_identical(receipt$features, rownames(SeuratObject::Loadings(result[["ATAClsi"]])))
  expect_true("ATAClsiclusters" %in% colnames(result[[]]))
  expect_error(suppressWarnings(RunStandardWorkflow(object,
    linear_reduction_dims_use = 2:30, nHVF = 80, neighbor_k = 5,
    nonlinear_reduction = character(), verbose = FALSE
  )), "dimensions are unavailable")
  for (dimensions in list(numeric(), c(1, NA), c(0, 1), c(1, 1), c(1, 1.5))) {
    expect_error(RunStandardWorkflow(object,
      linear_reduction_dims_use = dimensions,
      verbose = FALSE
    ), "unique positive integer")
  }
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
