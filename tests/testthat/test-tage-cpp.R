# Tests for tAge C++ backend — prediction
#
# tAge C++ code (src/tAge.cpp, 43 lines) contains:
#   tage_elastic_net_predict_cpp — Elastic Net prediction with imputation

# ---------------------------------------------------------------------------
# tage_elastic_net_predict_cpp
# ---------------------------------------------------------------------------

test_that("tage_elastic_net_predict_cpp returns correct structure", {
  n_samples <- 5
  n_features <- 3
  expr <- matrix(runif(n_features * n_samples), nrow = n_features, ncol = n_samples)

  out <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = as.integer(c(1, 2, 3)),
    imputer = c(0, 0, 0),
    center = c(0, 0, 0),
    scale = c(1, 1, 1),
    coef = c(0.2, 0.3, 0.5),
    intercept = 0.1
  )

  expect_type(out, "double")
  expect_length(out, n_samples)
  expect_true(all(is.finite(out)))
})

test_that("tage_elastic_net_predict_cpp with zero coefficients returns intercept", {
  n_samples <- 5
  n_features <- 3
  expr <- matrix(runif(n_features * n_samples), nrow = n_features, ncol = n_samples)

  out <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = as.integer(c(1, 2, 3)),
    imputer = c(0, 0, 0),
    center = c(0, 0, 0),
    scale = c(1, 1, 1),
    coef = c(0, 0, 0),
    intercept = 0.5
  )

  expect_equal(out, rep(0.5, n_samples), tolerance = 1e-10)
})

test_that("tage_elastic_net_predict_cpp applies centering and scaling", {
  n_samples <- 3
  # Single feature: values 2, 4, 6
  expr <- matrix(c(2, 4, 6), nrow = 1, ncol = n_samples)

  # center=2, scale=2: (x-2)/2, coef=1: pred = intercept + (x-2)/2
  # x=2: 0 + (2-2)/2 = 0
  # x=4: 0 + (4-2)/2 = 1
  # x=6: 0 + (6-2)/2 = 2
  out <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = 1L,
    imputer = 0,
    center = 2.0,
    scale = 2.0,
    coef = 1.0,
    intercept = 0.0
  )

  expect_equal(out, c(0, 1, 2), tolerance = 1e-10)
})

test_that("tage_elastic_net_predict_cpp handles missing features via imputation", {
  n_samples <- 3
  # Only 2 features in expression matrix
  expr <- matrix(c(
    1, 2, 3,
    4, 5, 6
  ), nrow = 2, ncol = n_samples, byrow = TRUE)

  # feature_match[3] = NA_INTEGER → use imputer[3] = 10
  out <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = as.integer(c(1, 2, NA)),
    imputer = c(0, 0, 10),
    center = c(0, 0, 0),
    scale = c(1, 1, 1),
    coef = c(0, 0, 1),  # only the missing feature has non-zero coef
    intercept = 0.0
  )

  # pred = 1 * 10 = 10 for all samples
  expect_equal(out, rep(10, n_samples), tolerance = 1e-10)
})

test_that("tage_elastic_net_predict_cpp multi-feature prediction", {
  n_samples <- 2
  # Two features: gene_1 (row 1), gene_2 (row 2)
  expr <- matrix(c(
    3, 7,   # gene_1
    5, 1    # gene_2
  ), nrow = 2, ncol = n_samples, byrow = TRUE)

  # pred = intercept + coef1*(x1-c1)/s1 + coef2*(x2-c2)/s2
  # center=0, scale=1, coef1=2, coef2=3, intercept=1
  # Sample 1: 1 + 2*3 + 3*5 = 1 + 6 + 15 = 22
  # Sample 2: 1 + 2*7 + 3*1 = 1 + 14 + 3 = 18
  out <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = as.integer(c(1, 2)),
    imputer = c(0, 0),
    center = c(0, 0),
    scale = c(1, 1),
    coef = c(2, 3),
    intercept = 1.0
  )

  expect_equal(out, c(22, 18), tolerance = 1e-10)
})

test_that("tage_elastic_net_predict_cpp empty scale vector works", {
  n_samples <- 3
  expr <- matrix(c(1, 2, 3), nrow = 1, ncol = n_samples)

  # Empty scale → treated as all 1.0
  out <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = 1L,
    imputer = 0,
    center = 0,
    scale = numeric(0),
    coef = 1.0,
    intercept = 0.0
  )

  expect_equal(out, c(1, 2, 3), tolerance = 1e-10)
})

test_that("tage_elastic_net_predict_cpp is deterministic", {
  n_samples <- 8
  expr <- matrix(runif(5 * n_samples), nrow = 5, ncol = n_samples)

  args <- list(
    expr = expr,
    feature_match = as.integer(c(1, 2, 3, 4, 5)),
    imputer = runif(5, 0, 1),
    center = runif(5, 0, 5),
    scale = runif(5, 0.5, 2),
    coef = runif(5, -1, 1),
    intercept = 0.3
  )

  out1 <- do.call(tage_elastic_net_predict_cpp, args)
  out2 <- do.call(tage_elastic_net_predict_cpp, args)
  expect_equal(out1, out2)
})

test_that("tage_elastic_net_predict_cpp rejects mismatched sizes", {
  expr <- matrix(1, 3, 5)
  expect_error(
    tage_elastic_net_predict_cpp(
      expr = expr,
      feature_match = as.integer(c(1, 2)),  # 2 features
      imputer = c(0, 0, 0),                  # 3 imputers (mismatch)
      center = c(0, 0),
      scale = c(1, 1),
      coef = c(0.5, 0.5),
      intercept = 0
    ),
    "same length"
  )
})
