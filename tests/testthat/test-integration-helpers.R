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
