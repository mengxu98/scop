test_that("RunIntegration requires append mode for multiple methods", {
  expect_error(
    RunIntegration(
      srt_merge = list(),
      batch = "batch",
      append = FALSE,
      integration_methods = c("Uncorrected", "Harmony"),
      verbose = FALSE
    ),
    "append = TRUE",
    fixed = TRUE
  )
})

test_that("RunIntegration dispatches the new plural argument", {
  testthat::local_mocked_bindings(
    run_integration_methods = function(args, integration_methods) {
      expect_false("integration_method" %in% names(args))
      expect_identical(args$integration_methods, integration_methods)
      integration_methods
    },
    .package = "scop"
  )

  result <- RunIntegration(
    srt_merge = list(),
    batch = "batch",
    integration_methods = c("Uncorrected", "Harmony"),
    verbose = FALSE
  )
  expect_identical(result, c("Uncorrected", "Harmony"))
})

test_that("multiple integration methods are run sequentially", {
  calls <- list()
  runner <- function(
    srt_merge = NULL,
    srt_list = NULL,
    integration_methods,
    append = TRUE,
    ...
  ) {
    calls[[length(calls) + 1L]] <<- list(
      srt_merge = srt_merge,
      srt_list = srt_list,
      integration_methods = integration_methods,
      append = append
    )
    paste(c(srt_merge, integration_methods), collapse = "+")
  }

  result <- scop:::run_integration_methods(
    args = list(
      srt_merge = "raw",
      srt_list = list("split"),
      append = TRUE
    ),
    integration_methods = c("Uncorrected", "Harmony", "CSS"),
    runner = runner
  )

  expect_identical(result, "raw+Uncorrected+Harmony+CSS")
  expect_identical(
    vapply(calls, `[[`, character(1), "integration_methods"),
    c("Uncorrected", "Harmony", "CSS")
  )
  expect_identical(calls[[1]]$srt_list, list("split"))
  expect_null(calls[[2]]$srt_list)
  expect_null(calls[[3]]$srt_list)
})

test_that("deprecated integration_method forwards until 1.0.0", {
  testthat::local_mocked_bindings(
    run_integration_methods = function(args, integration_methods) {
      integration_methods
    },
    .package = "scop"
  )

  expect_warning(
    result <- RunIntegration(
      srt_merge = list(),
      batch = "batch",
      integration_method = c("Uncorrected", "Harmony"),
      verbose = FALSE
    ),
    "removed in scop 1.0.0",
    fixed = TRUE
  )
  expect_identical(result, c("Uncorrected", "Harmony"))
})

test_that("new and deprecated integration arguments cannot be combined", {
  expect_error(
    RunIntegration(
      srt_merge = list(),
      batch = "batch",
      integration_methods = "Uncorrected",
      integration_method = "Harmony",
      verbose = FALSE
    ),
    "Supply only one",
    fixed = TRUE
  )
})
