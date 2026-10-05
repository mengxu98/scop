test_that("RCTD rejects invalid selected counts before rounding", {
  for (bad in c(-1, NA_real_, NaN, Inf, -Inf)) {
    for (label in c("Spatial", "Reference")) {
      for (rounding in c(TRUE, FALSE)) {
        mat <- Matrix::Matrix(matrix(c(bad, 2, 3, 4), 2,
          dimnames = list(c("g1", "g2"), c("s1", "s2"))), sparse = TRUE)
        testthat::local_mocked_bindings(GetAssayData5 = function(...) mat)
        expect_error(rctd_get_count_matrix(NULL, "RNA", "counts", c("g1", "g2"),
          data_label = label, round_counts = rounding, verbose = FALSE),
          paste0(label, ".*counts must contain only finite, non-negative values"))
      }
    }
  }
})

test_that("RCTD keeps the documented fractional rounding policy and identifiers", {
  mat <- Matrix::Matrix(matrix(c(1.2, 0, 2.8, 4), 2,
    dimnames = list(c("g1", "g2"), c("s2", "s1"))), sparse = TRUE)
  testthat::local_mocked_bindings(GetAssayData5 = function(...) mat)
  out <- rctd_get_count_matrix(NULL, "RNA", "counts", c("g2", "g1"), verbose = FALSE)
  expect_equal(as.matrix(out), round(as.matrix(mat[c("g2", "g1"), , drop = FALSE])))
  expect_error(rctd_get_count_matrix(NULL, "RNA", "counts", c("g1", "g2"),
    round_counts = FALSE), "non-integer")
  expect_equal(as.matrix(mat)[1, 1], 1.2)
})
