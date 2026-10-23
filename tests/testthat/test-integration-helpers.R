test_that("resolve_cluster_algorithm_index preserves Seurat algorithm indices", {
  expect_identical(resolve_cluster_algorithm_index("louvain"), 1)
  expect_identical(resolve_cluster_algorithm_index("slm"), 3)
  expect_identical(resolve_cluster_algorithm_index("leiden"), 4)
  expect_error(
    resolve_cluster_algorithm_index("louvain_refined"),
    "must be one of"
  )
  expect_error(
    resolve_cluster_algorithm_index("LEIDEN"),
    "must be one of"
  )
  expect_error(
    resolve_cluster_algorithm_index("unknown"),
    "must be one of"
  )
})

test_that("validate_nonlinear_reductions supports shared and restricted sets", {
  expect_invisible(
    validate_nonlinear_reductions(c("umap", "phate", "fr"))
  )
  expect_error(
    validate_nonlinear_reductions("invalid"),
    "must be one of"
  )

  restricted <- c("umap", "umap-naive", "fr")
  expect_invisible(
    validate_nonlinear_reductions("fr", allowed = restricted)
  )
  expect_error(
    validate_nonlinear_reductions("phate", allowed = restricted),
    "must be one of"
  )
})

test_that("validate_integration_input_cells preserves input checks", {
  srt_list <- list(
    batch_a = matrix(
      0,
      nrow = 1,
      ncol = 2,
      dimnames = list("gene", c("cell1", "cell2"))
    )
  )
  matching_merge <- matrix(
    0,
    nrow = 1,
    ncol = 2,
    dimnames = list("gene", c("cell2", "cell1"))
  )
  mismatched_merge <- matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list("gene", "cell3")
  )

  expect_true(validate_integration_input_cells(srt_list, NULL))
  expect_true(validate_integration_input_cells(NULL, matching_merge))
  expect_true(validate_integration_input_cells(srt_list, matching_merge))
  expect_error(
    validate_integration_input_cells(NULL, NULL),
    "were all empty"
  )
  expect_error(
    validate_integration_input_cells(srt_list, mismatched_merge),
    "different cells"
  )
})

test_that("Coralysis row-sparse input preserves native models and probabilities", {
  skip_if_not_installed("Coralysis")
  withr::local_options(lifecycle_verbosity = "quiet")
  set.seed(41)
  counts <- matrix(stats::rpois(100 * 240, lambda = 2), nrow = 100)
  dimnames(counts) <- list(paste0("g", 1:100), paste0("c", 1:240))
  object <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::ScaleData(object, features = rownames(object), verbose = FALSE)
  object$batch <- rep(c("a", "b"), each = 120)
  native_invoke <- invoke_fun
  checked <- FALSE

  testthat::local_mocked_bindings(
    invoke_fun = function(fun, params, ...) {
      if (!identical(fun, Coralysis::RunParallelDivisiveICP)) {
        return(native_invoke(fun, params, ...))
      }
      expect_s4_class(
        SummarizedExperiment::assay(params$object, "logcounts"),
        if (is.null(params$scale.by)) "dgRMatrix" else "dgCMatrix"
      )
      reference_params <- params
      SummarizedExperiment::assay(reference_params$object, "logcounts") <- methods::as(
        SummarizedExperiment::assay(reference_params$object, "logcounts"), "CsparseMatrix"
      )
      set.seed(11)
      expected <- native_invoke(fun, reference_params, ...)
      set.seed(11)
      actual <- native_invoke(fun, params, ...)
      expected_fit <- S4Vectors::metadata(expected)$coralysis
      actual_fit <- S4Vectors::metadata(actual)$coralysis
      expect_equal(actual_fit$joint.probability, expected_fit$joint.probability, tolerance = 1e-12)
      expect_equal(actual_fit$metrics, expected_fit$metrics, tolerance = 1e-12)
      expect_equal(
        lapply(actual_fit$models, function(model) model$W),
        lapply(expected_fit$models, function(model) model$W),
        tolerance = 1e-12
      )
      checked <<- TRUE
      actual
    },
    .package = "scop"
  )

  for (scale_by in list(NULL, "cell")) {
    checked <- FALSE
    result <- Coralysis_integrate(
      srt_merge = object, batch = "batch", nHVF = 80,
      nonlinear_reduction = character(), cores = 1,
      RunParallelDivisiveICP_params = list(
        L = 2, k = 4, scale.by = scale_by,
        build.train.params = list(nhvg = 80, p = 10)
      ),
      RunPCA_params = list(p = 2), verbose = FALSE
    )
    expect_true(checked)
    expect_identical(colnames(result), colnames(object))
    expect_true(all(is.finite(SeuratObject::Embeddings(result, "Coralysis"))))
  }
})
