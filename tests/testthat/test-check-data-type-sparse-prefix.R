sparse_prefix_type_capture <- function(code) {
  warnings <- character()
  value <- withCallingHandlers(code, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

test_that("sparse data type preserves bounds, missing data and late fractions", {
  mat <- function(v) Matrix::sparseMatrix(
    i = rep.int(1L, length(v)), j = seq_along(v), x = v,
    dims = c(1L, length(v))
  )
  status <- function(v) suppressMessages(CheckDataType(mat(v), verbose = FALSE))
  for (n in c(4096L, 4097L)) {
    x <- rep(1, n)
    expect_identical(status(x), "raw_counts")
    x[1L] <- 1.25
    expect_identical(status(x), "log_normalized_counts")
    x[1L] <- 1
    x[n] <- 1.25
    expect_identical(status(x), "log_normalized_counts")
  }
  x <- rep(1, 8301L)
  x[c(2L, 4L, 10L, 6300L, 7000L)] <- c(NA_real_, -0.0, .Machine$double.xmin, NaN, 0)
  expect_identical(status(x), "log_normalized_counts")
  x[6301L] <- -1
  expect_identical(suppressWarnings(status(x)), "unknown")
  x[6301L] <- Inf
  expect_identical(suppressWarnings(status(x)), "unknown")

  x <- rep(1, 8301L)
  x[1L] <- 1.25
  for (last in c(2^40 - 2^-13, 2^40, 2^40 + 2^-12, 2^52)) {
    x[length(x)] <- last
    expect_identical(status(x), "raw_normalized_counts")
  }

  x[length(x)] <- 1e30
  reference <- sparse_prefix_type_capture(any(x %% 1 != 0, na.rm = TRUE))
  result <- sparse_prefix_type_capture(status(x))
  expect_identical(result$value, "raw_normalized_counts")
  expect_identical(result$warnings, reference$warnings)
})

test_that("sparse data type preserves modulus warnings for spare-capacity matrices", {
  x <- rep(1, 8301L)
  x[1L] <- 1.25
  x[length(x)] <- 1e30
  overallocated <- Matrix::sparseMatrix(
    i = rep.int(1L, length(x)), j = seq_along(x), x = x,
    dims = c(1L, length(x))
  )
  j <- length(overallocated@p)
  overallocated@p[j] <- overallocated@p[j] - 1L
  if (!isTRUE(methods::validObject(overallocated, test = TRUE))) {
    skip("This Matrix version does not support spare capacity in the fixture")
  }
  reference <- sparse_prefix_type_capture(any(x %% 1 != 0, na.rm = TRUE))
  result <- sparse_prefix_type_capture(CheckDataType(overallocated, verbose = FALSE))
  expect_identical(result$value, "log_normalized_counts")
  expect_identical(result$warnings, reference$warnings)
})
