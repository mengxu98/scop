test_that("external_wrappers expands to external wrapper modules", {
  modules <- getFromNamespace("normalize_env_modules", "scop")("external_wrappers")

  expect_setequal(
    modules,
    c("scmalignantfinder", "secact", "scpagwas", "choir", "scanpy")
  )
})

test_that("scMalignantFinder module declares Python runtime requirements", {
  req <- scop::env_requirements(modules = "scmalignantfinder")

  expect_identical(
    unname(req$packages[["scMalignantFinder"]]),
    paste0(
      "git+https://github.com/Jonyyqn/scMalignantFinder.git@",
      "6ed094279c7adbbca30f2c73245fc074894d3715"
    )
  )
  expect_identical(unname(req$install_methods[["scMalignantFinder"]]), "pip")
  expect_identical(unname(req$packages[["numpy"]]), "numpy==1.26.4")
  expect_identical(unname(req$packages[["pyarrow"]]), "pyarrow==14.0.2")
  expect_false("xgboost" %in% names(req$packages))
  expect_false("scanpy" %in% names(req$packages))
})

test_that("external wrapper R packages remain optional explicit installs", {
  calls <- list()
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, dependencies, verbose, ...) {
      calls[[length(calls) + 1]] <<- list(
        packages = packages,
        dependencies = dependencies,
        verbose = verbose
      )
      invisible(TRUE)
    },
  )

  ok <- getFromNamespace("ensure_external_wrapper_r_packages", "scop")(
    modules = c("scmalignantfinder", "secact", "scpagwas", "choir"),
    verbose = FALSE
  )

  expect_equal(ok, TRUE)
  expect_equal(
    vapply(calls, `[[`, character(1), "packages"),
    c(
      "data2intelligence/SecAct",
      "sulab-wmu/scPagwas",
      "corceslab/CHOIR"
    )
  )
  expect_equal(
    vapply(calls, `[[`, logical(1), "dependencies"),
    c(NA, NA, NA)
  )
})

test_that("PrepareEnv exposes CHOIR as an optional R module", {
  packages <- getFromNamespace("env_r_packages", "scop")("choir")

  expect_true("corceslab/CHOIR" %in% packages)
})

test_that("PrepareEnv supports explicit R and Python components", {
  expect_identical(formals(scop::PrepareEnv)$components, "python")
  normalize <- getFromNamespace("norm_env_components", "scop")
  expect_setequal(normalize(c("python", "r")), c("python", "r"))
  expect_setequal(normalize("all"), c("python", "r"))

  packages <- getFromNamespace("env_r_packages", "scop")("secact")
  expect_true(all(c(
    "reticulate", "RcppAnnoy", "UCell",
    "jinworks/SpatialCellChat", "data2intelligence/SecAct"
  ) %in% packages))

  calls <- character()
  requested_cores <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, dependencies, cores, verbose, ...) {
      calls <<- c(calls, packages)
      requested_cores <<- cores
      stats::setNames(rep(TRUE, length(packages)), packages)
    }
  )
  status <- getFromNamespace("ensure_env_r_packages", "scop")(
    modules = "secact", cores = 2, verbose = FALSE
  )
  expect_named(status, packages)
  expect_true(is.logical(status))
  expect_setequal(calls, packages)
  expect_identical(requested_cores, 2)
})
