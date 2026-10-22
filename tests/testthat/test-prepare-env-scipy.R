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

test_that("healthy SciPy installations and other platforms are left alone", {
  calls <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    is_osx = function() TRUE,
    runner_system2 = function(command, args, env, stdout, stderr) {
      calls <<- calls + 1L
      writeLines("arm64 1.15.3", stdout)
      0L
    }
  )
  ensure_scipy_runtime("python")
  expect_identical(calls, 1L)
  testthat::local_mocked_bindings(.package = "scop", is_osx = function() FALSE)
  ensure_scipy_runtime("python")
  expect_identical(calls, 1L)
})

test_that("unrelated SciPy failures and failed replacements remain errors", {
  calls <- 0L
  output <- "ImportError: unrelated binary failure"
  testthat::local_mocked_bindings(
    .package = "scop",
    is_osx = function() TRUE,
    runner_system2 = function(command, args, env, stdout, stderr) {
      calls <<- calls + 1L
      writeLines(output, stdout)
      1L
    }
  )
  expect_error(ensure_scipy_runtime("python"), "unrelated binary failure")
  expect_identical(calls, 1L)
  output <- "arm64 1.15.3\n__thread_bss zero-fill section type"
  expect_error(ensure_scipy_runtime("python"), "SciPy wheel download failed")
  expect_identical(calls, 3L)
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

test_that("SciPy wheel installation and import failures cannot mark the runtime ready", {
  for (failed_step in c(3L, 4L)) {
    calls <- 0L
    testthat::local_mocked_bindings(
      .package = "scop",
      is_osx = function() TRUE,
      runner_system2 = function(command, args, env, stdout, stderr) {
        calls <<- calls + 1L
        args <- gsub("^'|'$", "", args)
        if (calls == 1L) {
          writeLines("arm64 1.15.3\n__thread_bss zero-fill section type", stdout)
          return(1L)
        }
        if ("download" %in% args) {
          file.create(file.path(args[match("--dest", args) + 1L], "scipy-1.15.3-cp311-cp311-macosx_12_0_arm64.whl"))
        }
        writeLines("controlled replacement failure", stdout)
        if (calls == failed_step) 1L else 0L
      }
    )
    expect_error(ensure_scipy_runtime("python"), "controlled replacement failure")
    expect_identical(calls, failed_step)
  }
})
