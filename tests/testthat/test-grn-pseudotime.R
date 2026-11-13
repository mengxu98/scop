test_that("time methods align by cell ID and preserve measured time", {
  expression <- matrix(seq_len(24), 4L, dimnames = list(letters[1:4], paste0("c", 1:6)))
  time <- matrix(c(2, 1, 3, 2, 4, 5), ncol = 1L, dimnames = list(colnames(expression), "trajectory"))
  input <- list(expression = expression, pseudotime = time, regulators = NULL, targets = NULL)
  config <- list(method = "scode", method_args = list())
  expected <- grn_time_data(input, config)
  input$expression <- expression[, c(6, 2, 5, 1, 4, 3)]
  input$pseudotime <- time[c(4, 3, 5, 2, 1, 6), , drop = FALSE]
  actual <- grn_time_data(input, config)
  expect_identical(actual, expected)
  expect_identical(as.numeric(actual$pseudotime), c(1, 2, 2, 3, 4, 5))
})

test_that("time methods reject missing or duplicate cell identities", {
  expression <- matrix(seq_len(24), 4L, dimnames = list(letters[1:4], paste0("c", 1:6)))
  time <- matrix(seq_len(6), ncol = 1L, dimnames = list(colnames(expression), "trajectory"))
  input <- list(expression = expression, pseudotime = time)
  config <- list(method = "scode", method_args = list())
  input$pseudotime <- time[-1L, , drop = FALSE]
  expect_error(grn_time_data(input, config), "matching expression")
  input$pseudotime <- time
  rownames(input$pseudotime)[2L] <- rownames(input$pseudotime)[1L]
  expect_error(grn_time_data(input, config), "unique cell IDs")
  input$pseudotime <- time
  input$pseudotime[2L, 1L] <- NA_real_
  expect_error(grn_time_data(input, config), "Every cell")
  input$pseudotime <- time
  input$pseudotime[] <- 1
  expect_error(grn_time_data(input, config), "two distinct")
})

test_that("LEAP preserves branches while single-trajectory methods require selection", {
  expression <- matrix(seq_len(24), 4L, dimnames = list(letters[1:4], paste0("c", 1:6)))
  time <- cbind(a = c(1:3, rep(NA_real_, 3)), b = c(rep(NA_real_, 3), 4:6))
  rownames(time) <- colnames(expression)
  input <- list(expression = expression, pseudotime = time)
  expect_identical(grn_time_data(input, list(method = "leap", method_args = list()))$pseudotime, time)
  expect_error(grn_time_data(input, list(method = "scode", method_args = list())), "Select one trajectory")
  expect_error(grn_time_data(input, list(method = "sincerities", method_args = list(trajectory = "a"))), "every selected cell")
  expect_error(grn_time_data(input, list(method = "scode", method_args = list(trajectory = "missing"))), "must name")
})
