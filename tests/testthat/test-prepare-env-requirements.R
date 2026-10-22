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
      args <- gsub("^[\"']|[\"']$", "", args)
      calls[[length(calls) + 1L]] <<- list(args = args, env = env)
      if (length(calls) == 1L) {
        writeLines(c("arm64 1.15.3", "section '__DATA/__thread_bss' has a zero-fill section type, but offset field is not zero"), stdout)
        return(1L)
      }
      if ("download" %in% args) {
        path <- args[match("--dest", args) + 1L]
        file.create(file.path(path, "scipy-1.15.3-cp311-cp311-macosx_12_0_arm64.whl"))
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
  testthat::local_mocked_bindings(.package = "scop",
    env_requirements = function(...) list(python = "3.11.8", packages = c(scenicplus = "pinned", rust = "rust")),
    build_env_cache_spec = function(requirements, ...) { seen <<- requirements; stop("runtime captured") }
  )
  expect_error(PrepareEnv(envname = "existing", conda = "conda", modules = "scenicplus_io", version = "3.11-1"), "runtime captured")
  expect_identical(seen$python, "3.11-1")
  expect_identical(seen$packages, c(scenicplus = "scenicplus"))
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
  expect_identical(requirements$python, "3.11.8")
  expect_identical(env_module_environment("scenicplus"), "scenicplus_env")
  expect_false(any(c("scipy", "pandas", "numba") %in% names(requirements$packages)))
  expect_true("scenicplus" %in% names(requirements$packages))
  git <- requirements$packages[grepl(" @ git\\+", requirements$packages)]
  expect_length(git, 5L)
  expect_true(all(grepl("@[0-9a-f]{40}$", git)))
  expect_error(env_requirements(modules = c("scenicplus", "scanpy")), "standalone")
  io <- env_requirements(modules = "scenicplus_io")
  expect_identical(env_module_environment("scenicplus_io"), "scenicplus_io_env")
  expect_identical(io$packages[["scipy"]], "scipy==1.15.3")
  expect_true("mudata" %in% names(io$packages))
  expect_false(any(c("scenicplus", "tbb", "numba", "leidenalg") %in% names(io$packages)))
  expect_false("scenicplus_io" %in% default_env_modules())
})

test_that("native Windows rejects full SCENIC+ before creating an environment", {
  testthat::local_mocked_bindings(.package = "scop", is_windows = function() TRUE,
    resolve_conda = function(...) stop("environment manager touched"))
  expect_error(PrepareEnv(modules = "scenicplus"), "full official SCENIC\\+ environment requires Linux")
})

test_that("Python preparation disables duplicate OpenMP loading by default", {
  withr::local_envvar(c(KMP_DUPLICATE_LIB_OK = NA_character_))
  configure_python_thread_env()
  expect_identical(Sys.getenv("KMP_DUPLICATE_LIB_OK"), "FALSE")
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "FALSE")
  configure_python_thread_env()
  expect_identical(Sys.getenv("KMP_DUPLICATE_LIB_OK"), "FALSE")
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  configure_python_thread_env()
  expect_identical(Sys.getenv("KMP_DUPLICATE_LIB_OK"), "TRUE")
})

test_that("full SCENIC+ installation activates conda and preserves dependency pins", {
  calls <- list()
  withr::local_envvar(c(UV_CACHE_DIR = NA_character_))
  previous <- Sys.getenv(c("UV_BUILD_CONSTRAINT", "PIP_BUILD_CONSTRAINT", "PIP_CONSTRAINT", "UV_CACHE_DIR", "UV_HTTP_TIMEOUT", "UV_CONCURRENT_DOWNLOADS", "CFLAGS", "CPATH", "LIBRARY_PATH", "PKG_CONFIG_PATH", "CC", "CXX"), unset = NA_character_)
  testthat::local_mocked_bindings(
    .package = "scop",
    find_uv = function(pip_options, ...) { expect_false("--overrides" %in% pip_options); "uv" },
    log_message = function(...) NULL,
    get_namespace_fun = function(...) function(command, args) {
      calls[[length(calls) + 1L]] <<- list(command = command, args = gsub("^[\"']|[\"']$", "", args),
        constraints = readLines(Sys.getenv("UV_BUILD_CONSTRAINT")), include = Sys.getenv("CPATH"), cache = Sys.getenv("UV_CACHE_DIR"))
      if (length(calls) == 1L) 0L else 1L
    }
  )
  expect_true(install_python_packages("scenicplus @ git+https://github.com/aertslab/scenicplus.git", "/env/bin/python", "/env", "conda"))
  expect_error(install_python_packages("scenicplus @ git+https://github.com/aertslab/scenicplus.git", "/env/bin/python", "/env", "micromamba"), "Pinned SCENIC\\+ installation via uv failed")
  expect_length(calls, 2L)
  expect_identical(vapply(calls, `[[`, "", "command"), c("conda", "micromamba"))
  expect_true(all(vapply(calls, function(x) all(c("run", "--prefix", "/env") %in% x$args), logical(1))))
  expect_true("--no-capture-output" %in% calls[[1L]]$args)
  expect_false("--no-capture-output" %in% calls[[2L]]$args)
  expect_true("--overrides" %in% calls[[1L]]$args)
  expect_true("setuptools<81" %in% calls[[1L]]$constraints)
  expect_true(startsWith(calls[[1L]]$include, "/env/include"))
  expect_identical(calls[[1L]]$cache, "/env/.scop-uv-cache")
  expect_identical(Sys.getenv(names(previous), unset = NA_character_), previous)
})

test_that("full SCENIC+ preparation verifies Git pins in an existing environment", {
  root <- withr::local_tempdir()
  testthat::local_mocked_bindings(.package = "scop",
    is_windows = function() FALSE,
    resolve_conda = function(...) "conda",
    conda_python = function(...) file.path(root, "bin", "python"),
    python_minor_version = function(...) "3.11",
    python_version_matches = function(...) TRUE,
    build_env_cache_spec = function(...) list(),
    is_cached_env_valid = function(...) FALSE,
    get_conda_envs_dir = function(...) root,
    conda_env_path = function(...) root,
    env_exist = function(...) TRUE,
    find_uv = function(...) "uv",
    log_message = function(...) NULL,
    installed_python_pkgs_pip = function(...) data.frame(package = character()),
    check_python = function(pip, ...) { if (pip) stop("Git pins skipped") else TRUE },
    install_python_packages = function(packages, ...) {
      expect_length(packages[grepl(" @ git\\+", packages)], 5L)
      stop("Git pins verified")
    }
  )
  expect_error(PrepareEnv(envname = root, conda = "conda", modules = "scenicplus", configure_runtime = FALSE), "Git pins verified")
})
