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

test_that("GNIPLR C++ preserves Python stable-order Granger semantics", {
  regulator <- rep(0:5, each = 5)
  within_tie <- rep(c(0, 1, 0, 2, 1), 6)
  expression <- rbind(
    regulator,
    2 * regulator + within_tie,
    5 - regulator + 0.5 * within_tie,
    2 * (regulator %% 3) + within_tie,
    c(tail(regulator, 3), head(regulator, -3)) + 0.25 * within_tie,
    ifelse((0:29) %% 2 == 0, regulator, 5 - regulator) + within_tie
  )

  result <- gniplr_cpp(
    expression = expression,
    target_idx = 1:3,
    correlation_threshold = 0.3,
    lasso_degree = 2L,
    lasso_alpha = 0.1,
    max_lag = 3L
  )
  adjacency <- as.data.frame(result$adjacency)
  keys <- paste(adjacency$regulator, adjacency$target, sep = "->")
  python_reference <- c(
    "3->1" = 1.243170864792076e-10,
    "3->2" = 3.065279319705049e-10,
    "4->1" = 5.467714782994445e-05,
    "4->2" = 2.632483710594411e-05,
    "4->3" = 3.107402195966683e-04,
    "5->1" = 1.149484643121255e-02,
    "5->2" = 7.376714193596356e-03,
    "5->3" = 2.246415159551971e-02,
    "6->2" = 9.744633084677119e-02,
    "6->3" = 3.857503077608526e-01
  )

  expect_true(all(names(python_reference) %in% keys))
  expect_equal(
    adjacency$pvalue[match(names(python_reference), keys)],
    unname(python_reference),
    tolerance = 1e-11
  )
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
    grnboost_python(expr, regulators = "g1", work_dir = tempdir()),
    "prepared"
  )
  expect_identical(modules, "scenic")

  expect_error(
    regdiffusion_python(expr, regulators = "g1", work_dir = tempdir()),
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
    run_scenicplus_python(
      srt = NULL,
      envname = "scenicplus_env",
      conda = "auto",
      scplus_object = "result.pkl"
    ),
    "prepared"
  )
  expect_identical(modules, "scenicplus")
})

test_that("CellRank fate confidence is derived from absorption probabilities", {
  fate_confidence <- getFromNamespace(
    "cellrank_fate_confidence_from_absorption",
    "scop"
  )
  absorption <- matrix(
    c(
      0.1, 0.7, 0.2,
      0.5, 0.25, 0.25,
      0.05, 0.15, 0.8
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(fate_confidence(absorption), c(0.7, 0.5, 0.8))
  expect_error(fate_confidence(matrix(NA_real_, 1, 1)), "non-finite")
})

test_that("Palantir restores existing and unset environment variables on failure", {
  withr::local_envvar(c(OMP_NUM_THREADS = "7", OPENBLAS_NUM_THREADS = NA, nm = NA))
  before <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "nm"), unset = NA_character_)
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) stop("controlled preparation failure")
  )
  expect_error(RunPalantir(backend = "python"), "controlled preparation failure")
  expect_identical(Sys.getenv(names(before), unset = NA_character_), before)
})

test_that("rejected Python runtime does not change thread configuration", {
  withr::local_envvar(c(OMP_NUM_THREADS = "7", NUMBA_NUM_THREADS = NA))
  before <- Sys.getenv(c("OMP_NUM_THREADS", "NUMBA_NUM_THREADS"), unset = NA_character_)
  python <- tempfile()
  file.create(python)
  on.exit(unlink(python), add = TRUE)
  testthat::local_mocked_bindings(
    .package = "scop",
    assert_python_runtime_switchable = function(...) stop("interpreter mismatch")
  )
  expect_error(configure_python_runtime(python), "interpreter mismatch")
  expect_identical(Sys.getenv(names(before), unset = NA_character_), before)
})

test_that("trajectory Python wrappers prepare through the public entry point", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) {
      seen <<- modules
      stop("public preparation")
    }
  )
  expect_error(RunPHATE(matrix(1:6, nrow = 3), backend = "python"), "public preparation")
  expect_identical(seen, "phate")
  expect_error(RunPAGA(adata = list(), group.by = "cluster", backend = "python"), "public preparation")
  expect_identical(seen, "scanpy")
  expect_error(RunSCVELO(backend = "python", magic_impute = FALSE), "public preparation")
  expect_identical(seen, "scvelo")
  expect_error(RunCellRank(backend = "python", magic_impute = FALSE), "public preparation")
  expect_true("cellrank" %in% seen)
  expect_error(RunPalantir(backend = "python"), "public preparation")
  expect_identical(seen, c("scanpy", "palantir"))
})

test_that("PrepareEnv reuses its managed environment cache", {
  python <- tempfile()
  file.create(python)
  on.exit(unlink(python), add = TRUE)
  withr::local_options(list(scop_env_cache = list(python = python)))
  configured <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    resolve_conda = function(conda) conda,
    is_cached_env_valid = function(spec) TRUE,
    assert_python_runtime_switchable = function(...) NULL,
    configure_python_runtime = function(python_path) configured <<- python_path,
    remember_python_environment = function(...) NULL,
    ensure_external_wrapper_r_packages = function(...) NULL,
    log_message = function(...) NULL
  )
  expect_null(PrepareEnv(envname = "test", conda = "test-conda",
    version = "3.10-1", modules = "phate", verbose = FALSE))
  expect_identical(configured, python)
})
