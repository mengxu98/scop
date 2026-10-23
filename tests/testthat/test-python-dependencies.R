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
  }
  for (version in c("3.11-1", "3.12-1")) {
    requirements <- env_requirements(version = version, modules = "moscot")$packages
    expect_identical(sub("jax==", "", requirements[["jax"]]), sub("jaxlib==", "", requirements[["jaxlib"]]))
  }
})

test_that("ordinary package lookup does not combine incompatible JAX profiles", {
  testthat::local_mocked_bindings(
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
