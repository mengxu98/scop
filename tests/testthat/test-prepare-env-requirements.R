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

mock_dependency_check <- function(output, status = 1L, closure = c("cellrank", "pygpcca", "jinja2", "pandas"), .env = parent.frame()) {
  testthat::local_mocked_bindings(
    get_namespace_fun = function(package, name) {
      expect_identical(c(package, name), c("reticulate", "system2t"))
      function(command, args, stdout = "", stderr = "") {
        if (identical(args, c("-m", "pip", "check"))) {
          writeLines(output, stdout)
          return(status)
        }
        writeLines(closure, stdout)
        writeLines(character(), stderr)
        0L
      }
    },
    .package = "scop",
    .env = .env
  )
}

test_that("Python conflicts in the transitive closure remain fatal", {
  mock_dependency_check(c(
    "pygpcca 1.0.4 has requirement jinja2==3.0.3, but you have jinja2 3.1.6.",
    "numpyro 0.20.1 has requirement jax>=0.7.0, but you have jax 0.4.38."
  ))
  expect_error(
    check_python_dependencies("python", "cellrank==2.0.7", "test-env"),
    "pygpcca.*jinja2"
  )
})

test_that("Python conflicts in directly requested packages remain fatal", {
  mock_dependency_check("cellrank 2.0.7 requires pygpcca, which is not installed.")
  expect_error(
    check_python_dependencies("python", "cellrank==2.0.7", "test-env"),
    "cellrank.*pygpcca"
  )
})

test_that("unrelated Python conflicts warn even when a dependency is shared", {
  withr::local_options(log_message.verbose = FALSE)
  mock_dependency_check(c(
    "numpyro 0.20.1 has requirement jax>=0.7.0, but you have jax 0.4.38.",
    "scprep 1.2.3 has requirement pandas<2.1,>=0.25, but you have pandas 2.2.0."
  ))
  expect_warning(
    result <- check_python_dependencies("python", "cellrank==2.0.7", "test-env", verbose = TRUE),
    "(?s)numpyro.*scprep.*PrepareEnv\\(force[[:space:]]*=[[:space:]]*TRUE\\)",
    perl = TRUE
  )
  expect_true(result)
})

test_that("Python dependency checks pass without conflicts", {
  mock_dependency_check("No broken requirements found.", status = 0L)
  expect_no_warning(
    result <- check_python_dependencies("python", "cellrank==2.0.7", "test-env")
  )
  expect_true(result)
})

test_that("unrecognized pip failures cannot be downgraded", {
  mock_dependency_check("ERROR: invalid distribution metadata")
  expect_error(
    check_python_dependencies("python", "cellrank==2.0.7", "test-env"),
    "invalid distribution metadata"
  )
})

test_that("unexpected pip exit codes remain fatal", {
  mock_dependency_check(
    "scprep 1.2.3 has requirement pandas<2.1, but you have pandas 2.2.0.",
    status = 2L
  )
  expect_error(
    check_python_dependencies("python", "cellrank==2.0.7", "test-env"),
    "Python dependency check failed"
  )
})

test_that("dependency metadata failures remain fatal", {
  testthat::local_mocked_bindings(
    get_namespace_fun = function(package, name) {
      function(command, args, stdout, stderr) {
        writeLines("", stdout)
        writeLines("Cannot read Requires-Dist", stderr)
        1L
      }
    },
    .package = "scop"
  )
  expect_error(
    python_dependency_closure("python", "cellrank"),
    "Cannot[[:space:]]+read[[:space:]]+Requires-Dist"
  )
})

test_that("Python installation preserves original exceptions with verbose disabled", {
  testthat::local_mocked_bindings(
    resolve_conda = identity,
    unique_python_requirements = function(packages, envname, conda) packages,
    env_exist = function(conda, envname) TRUE,
    exist_python_pkgs = function(packages, envname, conda, verbose) stats::setNames(FALSE, packages),
    conda_install = function(conda, packages, envname, pip, pip_options, ...) {
      stop("Original resolver failure: incompatible jaxlib wheel")
    },
    .package = "scop"
  )
  expect_error(
    check_python("jaxlib==0.4.38", envname = "test-env", conda = "conda", verbose = FALSE),
    "Original[[:space:]]+resolver[[:space:]]+failure:[[:space:]]+incompatible[[:space:]]+jaxlib[[:space:]]+wheel"
  )
})

test_that("scVI and Cell2location share a compatible JAX profile", {
  scvi <- env_requirements(version = "3.11-1", modules = "scvi")$packages
  location <- env_requirements(version = "3.11-1", modules = "cell2location")$packages
  expected <- c(
    jax = "jax[cpu]==0.4.38", jaxlib = "jaxlib==0.4.38",
    numpyro = "numpyro==0.17.0", optax = "optax==0.2.5",
    flax = "flax==0.10.4", chex = "chex==0.1.89",
    "orbax-checkpoint" = "orbax-checkpoint==0.8.0", tensorstore = "tensorstore==0.1.65"
  )
  expect_identical(scvi[names(expected)], expected)
  expect_identical(location[names(expected)], expected)
  fate <- env_requirements(modules = "cell2fate")$packages
  expect_identical(fate[["jax"]], "jax==0.4.10")
  expect_identical(fate[["jaxlib"]], "jaxlib==0.4.10")
})

test_that("moscot cannot silently override the scVI JAX profile", {
  for (version in c("3.10-1", "3.11-1", "3.12-1")) {
    expect_error(
      PrepareEnv(version = version, modules = c("scvi", "moscot")),
      "separate environments"
    )
    expect_error(
      PrepareEnv(version = version, modules = c("cell2location", "moscot")),
      "separate environments"
    )
  }
  for (version in c("3.11-1", "3.12-1")) {
    requirements <- env_requirements(version = version, modules = "moscot")$packages
    expect_identical(sub("jax==", "", requirements[["jax"]]), sub("jaxlib==", "", requirements[["jaxlib"]]))
  }
})

test_that("ordinary package lookup does not combine incompatible JAX profiles", {
  testthat::local_mocked_bindings(
    is_windows = function() FALSE,
    resolve_conda = identity,
    ensure_conda = function(conda) TRUE,
    env_exist = function(conda, envname) TRUE,
    conda_python = function(envname, conda) "python",
    python_minor_version = function(python_path) "3.11",
    installed_python_pkgs = function(envname, conda, verbose) {
      data.frame(package = "cellrank", version = "2.0.7")
    },
    build_env_cache_spec = function(...) stop("requirements resolved"),
    .package = "scop"
  )
  expect_identical(
    unname(resolve_requested_python_packages("jax", envname = "test-env", conda = "conda")),
    "jax[cpu]==0.4.38"
  )
  expect_identical(
    unname(resolve_requested_python_packages("moscot", envname = "test-env", conda = "conda")),
    "moscot==0.5.2"
  )
  expect_identical(
    exist_python_pkgs("cellrank==2.0.7", envname = "test-env", conda = "conda", verbose = FALSE),
    c("cellrank==2.0.7" = TRUE)
  )
  for (module in c("cellrank", "scvi", "moscot")) {
    expect_error(
      PrepareEnv(envname = "test-env", conda = "conda", version = "3.11-1", modules = module),
      "requirements resolved"
    )
  }
})

test_that("installed metadata follows transitive dependencies, extras, markers and cycles", {
  python <- Sys.which("python3")
  if (!nzchar(python)) skip("python3 is unavailable")
  probe <- suppressWarnings(system2(python, c("-c", shQuote("import packaging")), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(probe, "status"))) skip("Python packaging is unavailable")
  directory <- withr::local_tempdir()
  requires <- list(
    "closure-root" = c("closure-middle[feature]", "closure-middle", "closure-inactive; python_version < '2'"),
    "closure-middle" = c("closure-leaf", "closure-extra; extra == 'feature'", "closure-root"),
    "closure-leaf" = character(),
    "closure-extra" = character()
  )
  for (name in names(requires)) {
    path <- file.path(directory, paste0(gsub("-", "_", name), "-1.0.dist-info"))
    dir.create(path)
    writeLines(c(
      "Metadata-Version: 2.1", paste("Name:", name), "Version: 1.0",
      sprintf("Requires-Dist: %s", requires[[name]])
    ), file.path(path, "METADATA"))
  }
  withr::local_envvar(c(PYTHONPATH = directory))
  expect_setequal(
    python_dependency_closure(python, "Closure_Root==1.0"),
    c("closure-root", "closure-middle", "closure-leaf", "closure-extra")
  )
})
