test_that("SciBet lightweight prediction matches full probability path", {
  set.seed(1)
  ref <- matrix(stats::rpois(60 * 40, lambda = 1), nrow = 60)
  query <- matrix(stats::rpois(60 * 25, lambda = 1), nrow = 60)
  labels <- rep(seq_len(4), each = 10)

  full <- scibet_fit_predict(
    ref = ref,
    query = query,
    labels = labels,
    n_labels = 4L,
    n_top = 20L,
    additional_per_label = 2L,
    return_probabilities = TRUE
  )
  light <- scibet_fit_predict(
    ref = ref,
    query = query,
    labels = labels,
    n_labels = 4L,
    n_top = 20L,
    additional_per_label = 2L,
    return_probabilities = FALSE
  )

  expect_null(light$probabilities)
  expect_identical(light$predicted_index, full$predicted_index)
  expect_equal(light$max_probability, apply(full$probabilities, 1, max), tolerance = 1e-12)
})

test_that("SciBet sparse lightweight prediction matches full probability path", {
  set.seed(2)
  ref <- matrix(stats::rpois(60 * 40, lambda = 0.8), nrow = 60)
  query <- matrix(stats::rpois(60 * 25, lambda = 0.8), nrow = 60)
  labels <- rep(seq_len(4), each = 10)
  ref <- methods::as(Matrix::Matrix(ref, sparse = TRUE), "dgCMatrix")
  query <- methods::as(Matrix::Matrix(query, sparse = TRUE), "dgCMatrix")

  full <- scibet_fit_predict_sparse(
    ref = ref,
    query = query,
    labels = labels,
    n_labels = 4L,
    n_top = 20L,
    additional_per_label = 2L,
    return_probabilities = TRUE
  )
  light <- scibet_fit_predict_sparse(
    ref = ref,
    query = query,
    labels = labels,
    n_labels = 4L,
    n_top = 20L,
    additional_per_label = 2L,
    return_probabilities = FALSE
  )

  expect_null(light$probabilities)
  expect_identical(light$predicted_index, full$predicted_index)
  expect_equal(light$max_probability, apply(full$probabilities, 1, max), tolerance = 1e-12)
})
