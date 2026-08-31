# Reference implementation kept only as a test oracle for the C++ path.
knn_match_best_labels <- function(match_prob) {
  match_prob <- as.matrix(match_prob)
  if (!all(is.finite(match_prob))) {
    return(apply(
      match_prob,
      1L,
      function(x) names(x)[order(x, decreasing = TRUE)][1L]
    ))
  }
  out <- colnames(match_prob)[max.col(match_prob, ties.method = "first")]
  names(out) <- rownames(match_prob)
  out
}
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

test_that("matrix neighbor dispatch remains available to KNN reference mapping", {
  set.seed(17)
  reference <- matrix(
    stats::rnorm(80),
    nrow = 20,
    dimnames = list(paste0("ref", seq_len(20)), paste0("dim", seq_len(4)))
  )
  query <- matrix(
    stats::rnorm(40),
    nrow = 10,
    dimnames = list(paste0("query", seq_len(10)), paste0("dim", seq_len(4)))
  )

  observed <- Seurat::FindNeighbors(
    object = reference,
    query = query,
    k.param = 3,
    nn.method = "annoy",
    annoy.metric = "cosine",
    return.neighbor = TRUE,
    verbose = FALSE
  )

  expect_s4_class(observed, "Neighbor")
  expect_identical(dim(observed@nn.idx), c(10L, 3L))
  expect_identical(observed@cell.names, rownames(query))
  expect_true(all(observed@nn.idx >= 1L & observed@nn.idx <= nrow(reference)))
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
