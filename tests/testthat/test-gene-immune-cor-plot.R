test_that("GeneImmuneCorPlot variable-column filtering matches legacy variance checks", {
  values <- cbind(
    variable = c(1, 2, NA_real_, 4),
    constant = c(3, 3, 3, 3),
    all_missing = c(NA_real_, NA_real_, NA_real_, NA_real_)
  )

  legacy <- apply(values, 2, function(x) stats::var(x, na.rm = TRUE) > 0)
  fast <- scop:::gene_immune_variable_columns(values)

  expect_identical(fast, legacy)
})
