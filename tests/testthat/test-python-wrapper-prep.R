test_that("RunGNIPLR bootstraps its Python module before checking packages", {
  modules <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) {
      modules <<- modules
      stop("prepared")
    }
  )
  expr <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))
  )

  expect_error(
    RunGNIPLR(expr, backend = "python", correlation_threshold = 0),
    "prepared"
  )
  expect_identical(modules, "scanpy")
})

test_that("GRN Python wrappers bootstrap their isolated SCENIC environments", {
  modules <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) {
      modules <<- modules
      stop("prepared")
    }
  )
  expr <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    dimnames = list(paste0("cell", 1:3), c("g1", "g2"))
  )

  expect_error(
    scop:::grnboost_python(expr, regulators = "g1", work_dir = tempdir()),
    "prepared"
  )
  expect_identical(modules, "scenic")

  expect_error(
    scop:::regdiffusion_python(expr, regulators = "g1", work_dir = tempdir()),
    "prepared"
  )
  expect_identical(modules, c("scenic", "regdiffusion"))
})

test_that("MDIC3 Python backend selects the prepared Python runtime", {
  modules <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) {
      modules <<- modules
      stop("prepared")
    }
  )
  expr <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3"))
  )
  grn <- diag(2)
  dimnames(grn) <- list(rownames(expr), rownames(expr))

  expect_error(
    RunMDIC3(
      expr,
      labels = c("a", "a", "b"),
      grn = grn,
      backend = "python"
    ),
    "prepared"
  )
  expect_identical(modules, "scanpy")
})

test_that("SCENICPlus only prepares Python when an official object is supplied", {
  modules <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) {
      modules <<- modules
      stop("prepared")
    }
  )

  expect_error(
    scop:::run_scenicplus_python(
      srt = NULL,
      envname = "scenicplus_env",
      conda = "auto",
      scplus_object = "result.pkl"
    ),
    "prepared"
  )
  expect_identical(modules, "scenicplus")
})
