test_that("PrepareEnv module specs keep scanpy and scvelo distinct", {
  specs <- getFromNamespace("env_module_requirements", "scop")()
  scanpy <- specs$scanpy
  scvelo <- specs$scvelo

  expect_true("scanpy" %in% names(scanpy$packages))
  expect_false("scvelo" %in% names(scanpy$packages))
  expect_true("scvelo" %in% names(scvelo$packages))
  expect_false("scanpy" %in% names(scvelo$packages))
})

test_that("PrepareEnv resolves scvelo through modules rather than a dedicated helper", {
  req <- env_requirements(modules = "scvelo")
  expect_true("scanpy" %in% names(req$packages))
  expect_true("scvelo" %in% names(req$packages))
  expect_false(exists("scvelo_python_requirements", envir = asNamespace("scop"), inherits = FALSE))
  expect_false(exists("scanpy_python_requirements", envir = asNamespace("scop"), inherits = FALSE))
})

test_that("PrepareEnv retains Imports without per-method checks", {
  packages <- env_r_packages()
  expect_true(all(c(
    "reticulate", "ComplexHeatmap", "igraph", "Signac",
    "SummarizedExperiment", "S4Vectors", "thisplot"
  ) %in% packages))
})

test_that("scprep modules keep a separate compatible dependency profile", {
  expect_false(any(c("phate", "magic") %in% default_env_modules()))
  expect_identical(env_module_environment("sccoda"), "tensorflow_env")
  expect_identical(env_requirements(modules = "sccoda")$packages[["rpy2"]], "rpy2==3.5.17")
  expect_identical(env_module_environment("scomm"), "tensorflow_env")
  expect_identical(env_module_environment(c("scvelo", "magic")), "scprep_env")
  expect_identical(env_requirements(modules = c("scvelo", "magic"))$packages[["pandas"]], "pandas==2.0.3")
  expect_identical(env_requirements(modules = "phate")$packages[["phate"]], "phate==1.0.11")
  expect_identical(env_requirements(modules = "scrublet")$packages[["pandas"]], "pandas==2.2.0")
})

test_that("macOS repairs the affected SciPy wheel and verifies the replacement", {
  calls <- list()
  testthat::local_mocked_bindings(
    .package = "scop",
    is_osx = function() TRUE,
    runner_system2 = function(command, args, env, stdout, stderr) {
      args <- gsub("^'|'$", "", args)
      calls[[length(calls) + 1L]] <<- list(args = args, env = env)
      if (length(calls) == 1L) {
        writeLines(c("arm64 1.15.3", "section '__DATA/__thread_bss' has a zero-fill section type, but offset field is not zero"), stdout)
        return(1L)
      }
      if ("download" %in% args) {
        path <- args[match("--dest", args) + 1L]
        file.create(file.path(gsub("^'|'$", "", path), "scipy-1.15.3-cp311-cp311-macosx_12_0_arm64.whl"))
      }
      writeLines("arm64 1.15.3", stdout)
      0L
    }
  )
  expect_invisible(ensure_scipy_runtime("python", pip_options = c("--index-url", "https://example.org/simple", "--force-reinstall", "--upgrade", "--upgrade-strategy", "eager", "-f/local/wheels", "-ihttps://example.org/simple")))
  expect_length(calls, 4L)
  expect_true(all(c("--platform", "macosx_12_0_arm64", "scipy==1.15.3", "--index-url") %in% calls[[2L]]$args))
  expect_true(all(c("--no-deps", "--force-reinstall") %in% calls[[3L]]$args))
  expect_identical(unname(calls[[1L]]$env["KMP_DUPLICATE_LIB_OK"]), "FALSE")
  expect_true("-I" %in% calls[[1L]]$args)
  expect_identical(unname(calls[[1L]]$env["PYTHONPATH"]), "")
  expect_false(any(c("--force-reinstall", "--upgrade", "--upgrade-strategy", "eager") %in% calls[[2L]]$args))
  expect_true("https://example.org/simple" %in% calls[[2L]]$args)
  expect_true(all(c("-f/local/wheels", "-ihttps://example.org/simple") %in% calls[[2L]]$args))
})

test_that("reading an existing full SCENIC+ environment preserves its upstream pins", {
  python <- withr::local_tempfile()
  file.create(python)
  seen <- NULL
  testthat::local_mocked_bindings(
    .package = "scop",
    resolve_conda = function(...) "conda",
    conda_python = function(...) python,
    exist_python_pkgs = function(...) TRUE,
    env_requirements = function(modules, ...) {
      seen <<- modules
      stop("requirements captured")
    }
  )
  expect_error(PrepareEnv(envname = "existing", conda = "conda", modules = "scenicplus_io", version = "3.11-1"), "requirements captured")
  expect_identical(seen, "scenicplus")
  expect_error(PrepareEnv(envname = "existing", conda = "conda", modules = "scenicplus_io", version = "3.11-1", force = TRUE), "requirements captured")
  expect_identical(seen, "scenicplus_io")
})

test_that("subprocess preparation leaves an existing embedded runtime alone", {
  python <- withr::local_tempfile()
  file.create(python)
  withr::local_options(scop_env_cache = list(python = python))
  previous_path <- Sys.getenv("PATH", unset = NA_character_)
  testthat::local_mocked_bindings(
    .package = "scop",
    resolve_conda = function(conda) conda,
    is_cached_env_valid = function(...) TRUE,
    ensure_scipy_runtime = function(...) { Sys.setenv(PATH = "installer-path"); NULL },
    assert_python_runtime_switchable = function(...) stop("embedded runtime touched"),
    configure_python_runtime = function(...) stop("embedded runtime touched"),
    remember_python_environment = function(...) stop("embedded runtime touched"),
    ensure_external_wrapper_r_packages = function(...) NULL,
    log_message = function(...) NULL
  )
  expect_null(PrepareEnv(envname = "other", conda = "conda", modules = "tage", version = "3.10-1", configure_runtime = FALSE))
  expect_identical(Sys.getenv("PATH", unset = NA_character_), previous_path)
})

test_that("SCENIC+ keeps its upstream dependency stack separate", {
  requirements <- env_requirements(modules = "scenicplus")
  expect_identical(env_module_environment("scenicplus"), "scenicplus_env")
  expect_false(any(c("scipy", "pandas", "numba") %in% names(requirements$packages)))
  expect_true("scenicplus" %in% names(requirements$packages))
  expect_error(env_requirements(modules = c("scenicplus", "scanpy")), "standalone")
  io <- env_requirements(modules = "scenicplus_io")
  expect_identical(env_module_environment("scenicplus_io"), "scenicplus_io_env")
  expect_identical(io$packages[["scipy"]], "scipy==1.15.3")
  expect_true("mudata" %in% names(io$packages))
  expect_false(any(c("scenicplus", "tbb", "numba", "leidenalg") %in% names(io$packages)))
  expect_false("scenicplus_io" %in% default_env_modules())
})
