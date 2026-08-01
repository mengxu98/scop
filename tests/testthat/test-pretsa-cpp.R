test_that("PreTSA C++ block agrees with the R reference", {
  fit <- getFromNamespace("pretsa_one_block", "scop")
  set.seed(51)
  pseudotime <- sort(runif(180))
  names(pseudotime) <- paste0("c", seq_along(pseudotime))
  expression <- matrix(
    rnorm(75 * length(pseudotime)),
    nrow = 75,
    dimnames = list(paste0("g", seq_len(75)), names(pseudotime))
  )
  expression[1, ] <- 0
  expression[2, ] <- 3

  reference <- fit(
    expression,
    rownames(expression),
    pseudotime,
    pseudotime,
    knot = "auto",
    max_knot_allowed = 5,
    backend = "r",
    padjust_method = "fdr"
  )
  native <- fit(
    expression,
    rownames(expression),
    pseudotime,
    pseudotime,
    knot = "auto",
    max_knot_allowed = 5,
    backend = "cpp",
    padjust_method = "fdr"
  )

  expect_identical(native$fit_mat, reference$fit_mat)
  expect_identical(
    native$DynamicFeatures[rownames(reference$DynamicFeatures), ],
    reference$DynamicFeatures
  )
})
