test_that("KNN index-to-name conversion preserves the apply reference result", {
  indices <- matrix(
    c(4L, 1L, 2L, 3L, 2L, 4L),
    nrow = 3,
    dimnames = list(paste0("query_", 1:3), paste0("neighbor_", 1:2))
  )
  reference_names <- paste0("ref_", 1:4)

  expected <- apply(indices, c(1, 2), function(i) reference_names[i])
  actual <- scop:::knn_indices_to_names(indices, reference_names)

  expect_identical(actual, expected)
})

test_that("KNN index-to-name conversion is vectorized for large index matrices", {
  set.seed(1)
  indices <- matrix(sample.int(5000L, 5000L * 30L, replace = TRUE), ncol = 30L)
  reference_names <- paste0("ref_", seq_len(5000L))

  expect_identical(
    scop:::knn_indices_to_names(indices, reference_names),
    apply(indices, c(1, 2), function(i) reference_names[i])
  )
})
