test_that("numeric matrix validation retains plain double matrices and result schema", {
  score <- matrix(
    c(0, -0.0, -2, 3), nrow = 2L,
    dimnames = list(samples = c("s2", "s1"), methods = c("a", "b"))
  )
  # Includes dimnames labels, row and column order, and preservation of signed zero.
  expect_identical(numeric_matrix_result(score), score, num.eq = FALSE)
  corrupted <- score
  corrupted[] <- c(Inf, NaN, -Inf, 5)
  expected <- score
  expected[] <- c(NA_real_, NA_real_, NA_real_, 5)
  expect_identical(numeric_matrix_result(corrupted), expected)
  expect_identical(corrupted[[1L]], Inf)
  ints <- matrix(1:4, nrow = 2L, dimnames = list(c("s2", "s1"), NULL))
  expect_identical(
    numeric_matrix_result(ints),
    matrix(as.double(1:4), nrow = 2L,
           dimnames = list(c("s2", "s1"), c("score_1", "score_2")))
  )
  marked <- score
  attr(marked, "extra") <- "must not leak"
  expect_identical(numeric_matrix_result(marked), score, num.eq = FALSE)
  expect_false("numeric_matrix_result" %in% getNamespaceExports("scop"))
})
