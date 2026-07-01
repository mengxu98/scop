pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)

test_that("external_wrappers expands to external wrapper modules", {
  modules <- getFromNamespace("normalize_env_modules", "scop")("external_wrappers")

  expect_setequal(
    modules,
    c("scmalignantfinder", "secact", "scpagwas", "scanpy")
  )
})

test_that("scMalignantFinder module declares Python runtime requirements", {
  req <- env_requirements(modules = "scmalignantfinder")

  expect_identical(
    unname(req$packages[["scMalignantFinder"]]),
    "git+https://github.com/Jonyyqn/scMalignantFinder.git"
  )
  expect_identical(unname(req$packages[["xgboost"]]), "xgboost")
  expect_identical(unname(req$install_methods[["scMalignantFinder"]]), "pip")
  expect_identical(unname(req$install_methods[["xgboost"]]), "pip")
  expect_true("scanpy" %in% names(req$packages))
})

test_that("external wrapper R packages remain optional explicit installs", {
  calls <- list()
  testthat::local_mocked_bindings(
    check_r = function(packages, dependencies, verbose, ...) {
      calls[[length(calls) + 1]] <<- list(
        packages = packages,
        dependencies = dependencies,
        verbose = verbose
      )
      invisible(TRUE)
    }
  )

  ok <- getFromNamespace("ensure_external_wrapper_r_packages", "scop")(
    modules = c("scmalignantfinder", "secact", "scpagwas"),
    verbose = FALSE
  )

  expect_equal(ok, TRUE)
  expect_equal(
    vapply(calls, `[[`, character(1), "packages"),
    c("data2intelligence/SecAct", "sulab-wmu/scPagwas")
  )
  expect_equal(
    vapply(calls, `[[`, logical(1), "dependencies"),
    c(NA, NA)
  )
})
