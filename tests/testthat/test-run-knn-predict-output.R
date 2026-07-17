test_that("RunKNNPredict probability maxima retain row labels and missing values", {
  probabilities <- matrix(
    c(0.1, 0.7, NA_real_, 0.5, 0.4, NA_real_),
    nrow = 3,
    dimnames = list(c("Q1", "Q2", "Q3"), c("A", "B"))
  )
  testthat::local_mocked_bindings(check_r = function(...) TRUE, .package = "scop")

  expect_identical(
    knn_match_prob_max(probabilities),
    apply(probabilities, 1, max)
  )
})

test_that("KNN prediction labels retain legacy ties and non-finite values", {
  finite_probabilities <- rbind(
    Query1 = c(B = 0.8, A = 0.8, C = 0.1),
    Query2 = c(B = 0.2, A = 0.9, C = 0.3)
  )
  legacy_finite <- apply(
    finite_probabilities,
    1L,
    function(x) names(x)[order(x, decreasing = TRUE)][1L]
  )
  expect_identical(knn_match_best_labels(finite_probabilities), legacy_finite)

  probabilities <- rbind(
    finite_probabilities,
    Query3 = c(B = NA_real_, A = 0.5, C = 0.5)
  )
  legacy <- apply(
    probabilities,
    1L,
    function(x) names(x)[order(x, decreasing = TRUE)][1L]
  )

  expect_identical(knn_match_best_labels(probabilities), legacy)
})
