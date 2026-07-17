test_that("metacell native helpers preserve KNN graph and count aggregation", {
  knn_idx <- rbind(
    c(2L, 3L),
    c(1L, 3L),
    c(1L, 2L)
  )
  actual_graph <- metacell_knn_adjacency(knn_idx, n_cells = 3L)
  expected_graph <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L, 2L, 3L, 3L),
    j = c(2L, 3L, 1L, 3L, 1L, 2L),
    x = 1,
    dims = c(3L, 3L)
  )
  expect_equal(as.matrix(actual_graph), as.matrix(expected_graph))

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 4, 0, 0, 5, 6), nrow = 3L),
    sparse = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("c", 1:3))
  )
  membership <- c("MC2", "MC1", "MC2")
  actual_counts <- metacell_aggregate_counts(counts, membership)
  expected_counts <- cbind(
    MC2 = Matrix::rowSums(counts[, c(1L, 3L), drop = FALSE]),
    MC1 = Matrix::rowSums(counts[, 2L, drop = FALSE])
  )

  expect_equal(as.matrix(actual_counts), expected_counts)
  expect_identical(colnames(actual_counts), c("MC2", "MC1"))
})
