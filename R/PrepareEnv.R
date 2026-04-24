#' @title Prepare the python environment
#'
#' @description
#' Prepare the python environment by installing the required dependencies and setting up the environment.
#'
#' @md
#' @param envname The name of the conda environment.
#' If `NULL`, the environment name will be set to `"scop_env"`.
#' Default is `NULL`.
#' @param conda The path to a conda executable.
#' Use `"auto"` to allow automatically finding an appropriate conda binary.
#' @param miniconda_repo Repository URL for miniconda.
#' Default is \url{https://repo.anaconda.com/miniconda}.
#' @param force Whether to force recreation of the environment.
#' If `TRUE`, the existing environment will be removed and recreated.
#' Default is `FALSE`.
#' @param version The Python version.
#' Default is `"3.10-1"` on macOS and Unix and `"3.11-1"` on Windows.
#' @param modules Optional Python dependency modules to install in addition to
#' the default scientific stack.
#' Supported values are `"scanpy"`, `"scvi"`, `"glue"`, `"scanorama"`, `"bbknn"`,
#' `"celltypist"`, `"cellphonedb"`, `"magic"`, `"scrublet"`,
#' `"doubletdetection"`, `"doublet"`, `"palantir"`, `"scvelo"`,
#' `"cellrank"`, `"wot"`, `"phate"`, `"pacmap"`, `"trimap"`, `"multimap"`,
#' and `"scomm"`.
#' If `NULL` or omitted in [PrepareEnv()], the complete environment is installed.
#' @param pip_options Additional command line arguments to be passed to `uv`/`pip` when installing pip packages.
#' @param ... Additional arguments passed to package installation functions.
#'
#' @export
PrepareEnv <- function(
  envname = NULL,
  conda = "auto",
  miniconda_repo = "https://repo.anaconda.com/miniconda",
  version = if (is_windows()) "3.11-1" else "3.10-1",
  force = FALSE,
  modules = NULL,
  pip_options = character(),
  ...
) {
  log_message(
    "Preparing python environment...",
    text_color = "blue",
    message_type = "running"
  )

  if (!is.null(envname)) {
    options(scop_envname = envname)
  }

  envname <- get_envname(envname)
  modules <- normalize_env_modules(modules = modules)
  pip_options <- normalize_cli_args(pip_options)
  requirements <- env_requirements(
    version = version,
    modules = modules
  )
  python_version <- requirements[["python"]]

  log_message(
    "Environment name: {.file {envname}} and python version: {.pkg {python_version}}"
  )

  conda_info <- get_namespace_fun("reticulate", "conda_info")

  if (!is.null(conda)) {
    if (identical(conda, "auto")) {
      log_message("Auto-detecting conda...")
    }
    conda <- resolve_conda(conda)
  }

  cache_spec <- build_env_cache_spec(
    envname = envname,
    python_version = python_version,
    conda = conda,
    modules = modules,
    pip_options = pip_options,
    requirements = requirements
  )
  if (is_cached_env_valid(cache_spec) && isFALSE(force)) {
    python_cached <- getOption("scop_env_cache", default = NULL)[["python"]] %||%
      tryCatch(
        conda_python(envname = envname, conda = conda),
        error = function(...) NULL
      )
    if (!is.null(python_cached) && nzchar(python_cached) && file.exists(python_cached)) {
      configure_python_runtime(python_cached)
    }
    if ("scomm" %in% modules) {
      ensure_scomm_runtime_support(
        envname = envname,
        conda = conda,
        pip_options = pip_options,
        keep_jax = "scvi" %in% modules
      )
    }
    log_message(
      "{cli::col_green('Python environment completed')}\n",
      "{cli::col_grey('Until next loading, the environment will be cached')}",
      message_type = "success"
    )
    return(invisible(NULL))
  }

  if (is.null(conda)) {
    env <- FALSE
    log_message(
      "Conda not found, will install miniconda",
      message_type = "warning"
    )
  } else {
    envs_dir <- get_conda_envs_dir(conda = conda)
    env <- env_exist(
      conda = conda,
      envname = envname,
      envs_dir = envs_dir
    )

    if (isTRUE(force) && isTRUE(env)) {
      log_message("Force recreating environment...")
      unlink(paste0(envs_dir, "/", envname), recursive = TRUE)
      env <- FALSE
    }

    if (isTRUE(env)) {
      log_message(
        "Using existing environment: {.file {paste0(envs_dir, '/', envname)}}"
      )
    }
  }

  if (isFALSE(env)) {
    force <- TRUE

    if (is.null(conda)) {
      log_message("Installing miniconda...")
      conda <- install_miniconda2(miniconda_repo)
      cache_spec <- build_env_cache_spec(
        envname = envname,
        python_version = python_version,
        conda = conda,
        modules = modules,
        pip_options = pip_options,
        requirements = requirements
      )
    }

    log_message(
      "Creating conda environment with python {.val {python_version}}..."
    )

    accept_conda_tos(conda = conda)

    base_packages <- c("pip", "setuptools", "wheel")

    python_path <- reticulate::conda_create(
      conda = conda,
      envname = envname,
      python_version = python_version,
      packages = base_packages
    )

    python <- conda_python(envname = envname, conda = conda)

    envs_dir <- get_conda_envs_dir(conda = conda)
    env_path <- paste0(envs_dir, "/", envname)
    env <- file.exists(env_path)

    if (isFALSE(env)) {
      log_message(
        "Environment creation failed",
        message_type = "warning"
      )
      print(conda_info(conda = conda))
      print(reticulate::conda_list(conda = conda))
      log_message(
        "Unable to find environment under the expected path: {.file {env_path}}\n",
        "conda: {.pkg {conda}}\n",
        "python: {.file {python_path}}",
        message_type = "error"
      )
    } else {
      log_message(
        "Environment created successfully: {.file {env_path}}",
        message_type = "success"
      )
    }
  }

  install_methods <- requirements[["install_methods"]]
  packages <- requirements[["packages"]]

  python <- conda_python(envname = envname, conda = conda)
  uv <- find_uv(
    python = python,
    envname = envname,
    conda = conda,
    auto_install = TRUE,
    pip_options = pip_options
  )

  conda_packages <- packages[install_methods == "conda"]
  pip_packages <- packages[install_methods == "pip"]

  log_message(
    "{.val {length(packages)}} package{?s} to install, {.val {length(conda_packages)}} conda packages and {.val {length(pip_packages)}} pip packages"
  )

  if (length(conda_packages) > 0) {
    log_message(
      "Installing {.val {length(conda_packages)}} {.pkg conda} packages"
    )
    check_python(
      packages = conda_packages,
      envname = envname,
      conda = conda,
      force = force,
      pip = FALSE,
      verbose = FALSE,
      ...
    )
  }

  if (length(pip_packages) > 0) {
    if (!is.null(uv)) {
      log_message(
        "Installing {.val {length(pip_packages)}} packages using {.pkg uv}"
      )
    } else {
      log_message(
        "Installing {.val {length(pip_packages)}} {.pkg pip} packages"
      )
    }
    check_python(
      packages = pip_packages,
      envname = envname,
      conda = conda,
      force = force,
      pip = TRUE,
      pip_options = pip_options,
      verbose = FALSE,
      ...
    )
  }

  ensure_windows_scvi_support(
    envname = envname,
    conda = conda,
    force = force,
    pip_options = pip_options
  )

  if ("scomm" %in% modules) {
    ensure_scomm_runtime_support(
      envname = envname,
      conda = conda,
      pip_options = pip_options,
      keep_jax = "scvi" %in% modules
    )
  }

  configure_python_runtime(python)

  log_message(
    "{cli::col_green('Python environment ready')}\n",
    "{cli::col_grey('Until next loading, the environment will be cached')}",
    message_type = "success"
  )

  env_info(
    conda = conda,
    envname = envname,
    verbose = isTRUE(getOption("log_message.verbose", TRUE))
  )

  set_env_cache(cache_spec, python)
}

normalize_cli_args <- function(args) {
  if (length(args) == 0) {
    return(character())
  }

  args <- unlist(args, use.names = FALSE)
  if (length(args) == 1 && grepl("\\s", args)) {
    args <- strsplit(args, "\\s+")[[1]]
  }

  args <- trimws(as.character(args))
  args[nzchar(args)]
}

supported_env_modules <- function() {
  c(
    "scanpy",
    "scvi",
    "glue",
    "scanorama",
    "bbknn",
    "celltypist",
    "cellphonedb",
    "magic",
    "scrublet",
    "doubletdetection",
    "palantir",
    "scvelo",
    "cellrank",
    "wot",
    "phate",
    "pacmap",
    "trimap",
    "multimap",
    "scomm"
  )
}

optional_env_modules <- function() {
  c(
    "celltypist",
    "cellphonedb",
    "magic",
    "scrublet",
    "doubletdetection"
  )
}

env_module_dependencies <- function() {
  list(
    celltypist = "scanpy",
    cellphonedb = "scanpy",
    palantir = "scanpy",
    scvelo = "scanpy",
    cellrank = c("scanpy", "scvelo"),
    wot = "scanpy",
    scvi = "scanpy",
    glue = "scanpy",
    scanorama = "scanpy",
    bbknn = "scanpy",
    multimap = "scanpy"
  )
}

env_module_requirements <- function() {
  list(
    scanpy = scanpy_python_requirements(),
    scvi = scvi_python_requirements(),
    glue = glue_python_requirements(),
    scanorama = scanorama_python_requirements(),
    bbknn = bbknn_python_requirements(),
    celltypist = celltypist_python_requirements(),
    cellphonedb = cellphonedb_python_requirements(),
    magic = magic_python_requirements(),
    scrublet = scrublet_python_requirements(),
    doubletdetection = doubletdetection_python_requirements(),
    palantir = palantir_python_requirements(),
    scvelo = scvelo_python_requirements(),
    cellrank = cellrank_python_requirements(),
    wot = wot_python_requirements(),
    phate = phate_python_requirements(),
    pacmap = pacmap_python_requirements(),
    trimap = trimap_python_requirements(),
    multimap = multimap_python_requirements(),
    scomm = scomm_python_requirements()
  )
}

normalize_env_modules <- function(modules = NULL, include_optional = FALSE) {
  modules <- modules %||% character(0)
  modules <- unlist(modules, use.names = FALSE)
  modules <- trimws(as.character(modules))
  modules <- modules[nzchar(modules)]

  full_modules <- supported_env_modules()

  if (length(modules) == 0) {
    modules <- full_modules
  }

  if (isTRUE(include_optional)) {
    modules <- c(modules, optional_env_modules())
  }

  if ("optional" %in% modules) {
    modules <- c(setdiff(modules, "optional"), optional_env_modules())
  }

  if ("doublet" %in% modules) {
    modules <- c(
      setdiff(modules, "doublet"),
      "scrublet",
      "doubletdetection"
    )
  }

  module_dependencies <- env_module_dependencies()
  for (module in modules) {
    deps <- module_dependencies[[module]]
    if (!is.null(deps)) {
      modules <- c(modules, deps)
    }
  }

  modules <- unique(modules)
  invalid_modules <- setdiff(modules, c(supported_env_modules(), "optional", "doublet"))
  if (length(invalid_modules) > 0) {
    log_message(
      "{.arg modules} contains unsupported values: {.val {invalid_modules}}",
      message_type = "error"
    )
  }
  modules
}

build_env_cache_spec <- function(
  envname,
  python_version,
  conda,
  modules,
  pip_options,
  requirements = NULL
) {
  conda_path <- if (is.null(conda) || identical(conda, "")) {
    NULL
  } else {
    normalizePath(conda, mustWork = FALSE)
  }

  list(
    envname = get_envname(envname),
    python_version = python_version,
    conda = conda_path,
    modules = normalize_env_modules(modules = modules),
    pip_options = normalize_cli_args(pip_options),
    requirements = requirements[c("packages", "install_methods")]
  )
}

clear_env_cache <- function() {
  options(scop_env_cache = NULL)
  invisible(NULL)
}

set_env_cache <- function(spec, python = NULL) {
  cache <- spec
  if (!is.null(python) && nzchar(python) && file.exists(python)) {
    cache[["python"]] <- normalizePath(python, mustWork = FALSE)
  }
  options(scop_env_cache = cache)
  invisible(cache)
}

is_cached_env_valid <- function(spec) {
  cache <- getOption("scop_env_cache", default = NULL)
  required_fields <- c(
    "envname",
    "python_version",
    "conda",
    "modules",
    "pip_options",
    "requirements"
  )

  if (!is.list(cache) || !all(required_fields %in% names(cache))) {
    return(FALSE)
  }

  if (
    !identical(unname(cache[required_fields]), unname(spec[required_fields]))
  ) {
    return(FALSE)
  }

  if (is.null(spec[["conda"]]) || !file.exists(spec[["conda"]])) {
    return(FALSE)
  }

  if (!env_exist(conda = spec[["conda"]], envname = spec[["envname"]])) {
    return(FALSE)
  }

  python <- cache[["python"]] %||%
    tryCatch(
      conda_python(envname = spec[["envname"]], conda = spec[["conda"]]),
      error = function(e) NULL
    )

  !is.null(python) && file.exists(python)
}

prepend_path_var <- function(var, values) {
  values <- values[nzchar(values)]
  values <- values[file.exists(values)]
  if (length(values) == 0) {
    return(invisible(NULL))
  }

  values <- unique(normalizePath(values, mustWork = FALSE))
  current <- Sys.getenv(var, unset = "")
  current_values <- strsplit(current, .Platform$path.sep, fixed = TRUE)[[1]]
  current_values <- current_values[nzchar(current_values)]
  do.call(
    Sys.setenv,
    stats::setNames(
      list(paste(unique(c(values, current_values)), collapse = .Platform$path.sep)),
      var
    )
  )
  invisible(NULL)
}

configure_python_runtime <- function(python_path) {
  python_path <- normalizePath(python_path, mustWork = FALSE)
  if (!nzchar(python_path) || !file.exists(python_path)) {
    return(invisible(FALSE))
  }

  python_dir <- dirname(python_path)
  env_path <- dirname(python_dir)

  Sys.setenv(
    RETICULATE_PYTHON = python_path,
    PYTHONNOUSERSITE = "1",
    PIP_USER = "0"
  )
  Sys.unsetenv("PYTHONPATH")

  prepend_path_var(
    "PATH",
    c(
      python_dir,
      file.path(env_path, "bin"),
      file.path(env_path, "Library", "bin"),
      file.path(env_path, "Scripts")
    )
  )

  if (is_linux()) {
    prepend_path_var(
      "LD_LIBRARY_PATH",
      c(
        file.path(env_path, "lib"),
        file.path(env_path, "lib64"),
        file.path(env_path, "Library", "lib")
      )
    )
  }

  if (is_osx()) {
    prepend_path_var(
      "DYLD_FALLBACK_LIBRARY_PATH",
      c(
        file.path(env_path, "lib"),
        file.path(env_path, "lib64")
      )
    )
  }

  invisible(TRUE)
}

set_python_env <- function(conda, envname, verbose = TRUE) {
  Sys.unsetenv("RETICULATE_PYTHON")
  options(reticulate.miniconda.enabled = FALSE)

  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
  Sys.setenv(NUMEXPR_NUM_THREADS = "1")
  Sys.setenv(KMP_WARNINGS = "0")
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  Sys.setenv(NUMBA_NUM_THREADS = "1")

  python_path <- conda_python(
    conda = conda,
    envname = envname
  )

  configure_python_runtime(python_path)

  reticulate::use_python(
    python_path,
    required = TRUE
  )

  tryCatch(
    {
      reticulate::py_run_string(
        "
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['KMP_WARNINGS'] = '0'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
os.environ['NUMBA_NUM_THREADS'] = '1'
try:
    import numba
    if hasattr(numba, 'set_num_threads'):
        try:
            numba.set_num_threads(1)
        except RuntimeError:
            pass
    if hasattr(numba, 'config'):
        numba.config.NUMBA_NUM_THREADS = 1
except ImportError:
    pass
except Exception:
    pass
"
      )
    },
    error = function(e) {}
  )
}

ensure_scomm_runtime_support <- function(
  envname,
  conda,
  pip_options = character(),
  keep_jax = FALSE,
  verbose = TRUE
) {
  installed <- tryCatch(
    installed_python_pkgs(envname = envname, conda = conda),
    error = function(...) NULL
  )
  if (is.null(installed) || is.null(installed$package)) {
    return(invisible(FALSE))
  }

  removable <- if (isTRUE(keep_jax)) {
    character()
  } else {
    intersect(c("jax", "jaxlib"), installed$package)
  }
  if (length(removable) > 0) {
    log_message(
      "Removing conflicting packages for {.pkg scOMM}: {.pkg {removable}}",
      verbose = verbose
    )
    remove_python(
      packages = removable,
      envname = envname,
      conda = conda,
      pip = TRUE,
      force = TRUE,
      verbose = verbose
    )
  }

  check_python(
    packages = if (isTRUE(keep_jax)) {
      "ml-dtypes>=0.3.2"
    } else {
      "ml-dtypes==0.3.2"
    },
    envname = envname,
    conda = conda,
    force = FALSE,
    pip = TRUE,
    pip_options = pip_options,
    verbose = verbose
  )

  invisible(TRUE)
}

install_miniconda2 <- function(
  miniconda_repo,
  timeout = 600
) {
  log_message("Installing miniconda...")
  options(timeout = timeout)

  info <- as.list(Sys.info())

  if (is_apple_silicon()) {
    base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
    name <- "Miniforge3-MacOSX-arm64.sh"
    url <- file.path(base, name)
  } else {
    version <- "3"
    arch <- get_namespace_fun(
      "reticulate",
      "miniconda_installer_arch"
    )(info)

    name <- if (is_windows()) {
      sprintf("Miniconda%s-latest-Windows-%s.exe", version, arch)
    } else if (is_osx()) {
      sprintf("Miniconda%s-latest-MacOSX-%s.sh", version, arch)
    } else if (is_linux()) {
      sprintf("Miniconda%s-latest-Linux-%s.sh", version, arch)
    } else {
      log_message(
        "Unsupported platform: {.val {info$sysname}}",
        message_type = "error"
      )
    }

    url <- file.path(miniconda_repo, name)
  }

  options(reticulate.miniconda.url = url)

  if (!is.na(Sys.getenv("USER", unset = NA))) {
    miniconda_path <- gsub(
      pattern = "\\$USER",
      replacement = Sys.getenv("USER"),
      reticulate::miniconda_path()
    )
  } else {
    miniconda_path <- reticulate::miniconda_path()
  }

  if (dir.exists(miniconda_path)) {
    log_message("Removing existing miniconda installation...")
    unlink(miniconda_path, recursive = TRUE)
  }

  log_message("Downloading and installing miniconda from: {.url {url}}")
  reticulate::install_miniconda(
    path = miniconda_path,
    force = TRUE,
    update = FALSE
  )

  conda <- get_namespace_fun(
    "reticulate",
    "miniconda_conda"
  )(miniconda_path)
  log_message("Miniconda installed at: {.file {miniconda_path}}")

  conda
}

#' @title Print environment information
#' @inheritParams PrepareEnv
#' @param verbose Whether to print environment information.
#' @md
env_info <- function(conda, envname, verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  envs_dir <- get_conda_envs_dir(conda = conda)

  python_path <- tryCatch(
    conda_python(envname = envname, conda = conda),
    error = function(...) NULL
  )
  python_path <- python_path %||% paste0(envs_dir, "/", get_envname(envname), "/bin/python")
  py_version <- if (!is.null(python_path) && nzchar(python_path) && file.exists(python_path)) {
    out <- tryCatch(
      suppressWarnings(system2(python_path, "--version", stdout = TRUE, stderr = TRUE)),
      error = function(...) character()
    )
    out[[1]] %||% "unknown"
  } else {
    "unknown"
  }

  py_info_mesg <- c(
    cli::col_blue(
      "Conda config:"
    ),
    cli::col_grey(
      " Conda:         ",
      conda
    ),
    cli::col_grey(
      " Environment:   ",
      envs_dir,
      "/",
      get_envname(envname)
    ),
    cli::col_blue(
      "Python config:"
    ),
    cli::col_grey(
      " python:        ",
      python_path
    ),
    cli::col_grey(
      " version:       ",
      py_version
    )
  )
  invisible(lapply(py_info_mesg, packageStartupMessage))
}

#' @title Python environment requirements
#'
#' @description
#' The function returns a list of requirements including the required Python version,
#' package versions, and package name aliases for platform-specific packages.
#' All packages will be installed using uv as the primary tool.
#'
#' @md
#' @param version The Python version of the environment.
#' Default is `"3.10-1"`.
#' @param include_optional Whether to include optional Python dependencies.
#' @param modules Optional requirement modules to include. Supported values are
#' `"scanpy"`, `"scvi"`, `"scanorama"`, `"bbknn"`, `"celltypist"`,
#' `"cellphonedb"`, `"magic"`, `"scrublet"`, `"doubletdetection"`,
#' `"doublet"`, `"palantir"`, `"scvelo"`, `"cellrank"`, `"wot"`, `"phate"`,
#' `"pacmap"`, `"trimap"`, `"multimap"`, and `"scomm"`. If `NULL`, the
#' complete environment is returned.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{python}{Python version string}
#'   \item{packages}{Named vector of package version specifications}
#'   \item{package_aliases}{Named list mapping logical package names to actual installed names}
#' }
#'
#' @export
#' @examples
#' env_requirements("3.10-1")
env_requirements <- function(
  version = "3.10-1",
  include_optional = FALSE,
  modules = NULL
) {
  version <- match.arg(
    version,
    choices = c("3.10-1", "3.11-1")
  )
  modules <- normalize_env_modules(
    modules = modules,
    include_optional = include_optional
  )

  base_requirements <- core_python_requirements()
  package_install_methods <- base_requirements$install_methods
  package_versions <- base_requirements$packages
  package_aliases <- base_requirements$package_aliases

  module_requirements <- env_module_requirements()

  for (module in modules) {
    req_i <- module_requirements[[module]]
    if (is.null(req_i)) {
      next
    }
    package_install_methods <- c(package_install_methods, req_i$install_methods)
    package_versions <- c(package_versions, req_i$packages)
    package_aliases <- c(package_aliases, req_i$package_aliases)
  }

  if (all(c("scvi", "scomm") %in% modules) && "ml_dtypes" %in% names(package_versions)) {
    package_versions[["ml_dtypes"]] <- "ml-dtypes>=0.3.2"
  }

  requirements <- list(
    python = version,
    packages = package_versions,
    install_methods = package_install_methods,
    package_aliases = package_aliases
  )

  return(requirements)
}

core_python_requirements <- function() {
  tbb_install_method <- if (is_apple_silicon()) {
    "pip"
  } else {
    "conda"
  }
  list(
    packages = c(
      "leidenalg" = "leidenalg==0.10.2",
      "tbb" = if (is_apple_silicon()) {
        "pxr-tbb==2022.2.0"
      } else {
        "tbb==2022.2.0"
      },
      "python-igraph" = "python-igraph==0.11.9",
      "matplotlib" = "matplotlib==3.10.8",
      "numba" = "numba==0.59.1",
      "llvmlite" = "llvmlite==0.42.0",
      "numpy" = "numpy==1.26.4",
      "packaging" = "packaging>=24.0",
      "pandas" = "pandas==2.0.3",
      "scikit-learn" = "scikit-learn==1.7.0",
      "scipy" = "scipy==1.15.3",
      if (is_apple_silicon()) {
        c("llvm-openmp" = "llvm-openmp>=17")
      } else {
        NULL
      }
    ),
    install_methods = c(
      "leidenalg" = "conda",
      "tbb" = tbb_install_method,
      "python-igraph" = if (is_apple_silicon()) {
        "pip"
      } else {
        "conda"
      },
      "matplotlib" = "pip",
      "numba" = "pip",
      "llvmlite" = "pip",
      "numpy" = "pip",
      "packaging" = "pip",
      "pandas" = "pip",
      "scikit-learn" = "pip",
      "scipy" = "pip",
      if (is_apple_silicon()) {
        c("llvm-openmp" = "conda")
      } else {
        NULL
      }
    ),
    package_aliases = c(
      list("python-igraph" = "igraph"),
      if (is_apple_silicon()) {
        list("tbb" = "pxr-tbb")
      } else {
        list()
      }
    )
  )
}

scanpy_python_requirements <- function() {
  list(
    packages = c(
      "scanpy" = "scanpy==1.11.3"
    ),
    install_methods = c(
      "scanpy" = "pip"
    ),
    package_aliases = list()
  )
}

scvi_python_requirements <- function() {
  list(
    packages = c(
      "scvi-tools" = "scvi-tools==1.2.1",
      "jax" = "jax[cpu]==0.4.38"
    ),
    install_methods = c(
      "scvi-tools" = "conda",
      "jax" = "pip"
    ),
    package_aliases = list()
  )
}

glue_python_requirements <- function() {
  list(
    packages = c(
      "scglue" = "scglue==0.4.0",
      "bedtools" = "bedtools"
    ),
    install_methods = c(
      "scglue" = "pip",
      "bedtools" = "conda"
    ),
    package_aliases = list()
  )
}

scanorama_python_requirements <- function() {
  list(
    packages = c(
      "scanorama" = "scanorama==1.7.4"
    ),
    install_methods = c(
      "scanorama" = "pip"
    ),
    package_aliases = list()
  )
}

bbknn_python_requirements <- function() {
  list(
    packages = c(
      "bbknn" = "bbknn==1.6.0"
    ),
    install_methods = c(
      "bbknn" = "pip"
    ),
    package_aliases = list()
  )
}

celltypist_python_requirements <- function() {
  list(
    packages = c(
      "celltypist" = "celltypist==1.7.1"
    ),
    install_methods = c(
      "celltypist" = "pip"
    ),
    package_aliases = list()
  )
}

cellphonedb_python_requirements <- function() {
  list(
    packages = c(
      "cellphonedb" = "cellphonedb==5.0.1"
    ),
    install_methods = c(
      "cellphonedb" = "pip"
    ),
    package_aliases = list()
  )
}

magic_python_requirements <- function() {
  list(
    packages = c(
      "magic-impute" = "magic-impute==3.0.0"
    ),
    install_methods = c(
      "magic-impute" = "pip"
    ),
    package_aliases = list()
  )
}

scrublet_python_requirements <- function() {
  list(
    packages = c(
      "scrublet" = "scrublet==0.2.3"
    ),
    install_methods = c(
      "scrublet" = "pip"
    ),
    package_aliases = list()
  )
}

doubletdetection_python_requirements <- function() {
  list(
    packages = c(
      "celltypist" = "celltypist==1.7.1",
      "cellphonedb" = "cellphonedb==5.0.1",
      "sccoda" = "sccoda>=0.1.9",
      "magic-impute" = "magic-impute==3.0.0",
      "scrublet" = "scrublet==0.2.3",
      "doubletdetection" = "doubletdetection==4.3.0.post1",
      "louvain" = "louvain==0.8.2"
    ),
    install_methods = c(
      "celltypist" = "pip",
      "cellphonedb" = "pip",
      "sccoda" = "pip",
      "magic-impute" = "pip",
      "scrublet" = "pip",
      "doubletdetection" = "pip",
      "louvain" = "pip"
    ),
    package_aliases = list()
  )
}

palantir_python_requirements <- function() {
  list(
    packages = c(
      "palantir" = "palantir==1.4.1"
    ),
    install_methods = c(
      "palantir" = "pip"
    ),
    package_aliases = list()
  )
}

scvelo_python_requirements <- function() {
  list(
    packages = c(
      "scvelo" = "scvelo==0.3.3"
    ),
    install_methods = c(
      "scvelo" = "pip"
    ),
    package_aliases = list()
  )
}

cellrank_python_requirements <- function() {
  list(
    packages = c(
      "cellrank" = "cellrank==2.0.7"
    ),
    install_methods = c(
      "cellrank" = "pip"
    ),
    package_aliases = list()
  )
}

wot_python_requirements <- function() {
  list(
    packages = c(
      "wot" = "wot==1.0.8.post2"
    ),
    install_methods = c(
      "wot" = "pip"
    ),
    package_aliases = list()
  )
}

phate_python_requirements <- function() {
  list(
    packages = c(
      "phate" = "phate==1.0.11"
    ),
    install_methods = c(
      "phate" = "pip"
    ),
    package_aliases = list()
  )
}

pacmap_python_requirements <- function() {
  list(
    packages = c(
      "pacmap" = "pacmap==0.8.0"
    ),
    install_methods = c(
      "pacmap" = "pip"
    ),
    package_aliases = list()
  )
}

trimap_python_requirements <- function() {
  list(
    packages = c(
      "trimap" = "trimap==1.1.4"
    ),
    install_methods = c(
      "trimap" = "pip"
    ),
    package_aliases = list()
  )
}

multimap_python_requirements <- function() {
  list(
    packages = c(
      "multimap" = "git+https://github.com/Teichlab/MultiMAP.git"
    ),
    install_methods = c(
      "multimap" = "pip"
    ),
    package_aliases = list(
      "multimap" = "MultiMAP",
      "MultiMAP" = "multimap"
    )
  )
}

scomm_python_requirements <- function() {
  list(
    packages = c(
      "tensorflow" = "tensorflow==2.16.2",
      "keras" = "keras==3.3.3",
      "tf_keras" = "tf_keras==2.16.0",
      "ml_dtypes" = "ml-dtypes==0.3.2"
    ),
    install_methods = c(
      "tensorflow" = "pip",
      "keras" = "pip",
      "tf_keras" = "pip",
      "ml_dtypes" = "pip"
    ),
    package_aliases = list(
      "tf_keras" = "tf-keras",
      "tf-keras" = "tf_keras",
      "ml_dtypes" = "ml-dtypes",
      "ml-dtypes" = "ml_dtypes"
    )
  )
}

env_exist <- function(
  conda = "auto",
  envname = NULL,
  envs_dir = NULL
) {
  envname <- get_envname(envname)
  conda <- resolve_conda(conda)

  if (!is.null(conda)) {
    if (is.null(envs_dir)) {
      conda_info <- tryCatch(
        {
          get_namespace_fun("reticulate", "conda_info")(conda = conda)
        },
        error = function(e) {
          log_message(
            "Failed to get conda info: {.val {e$message}}",
            message_type = "warning"
          )
          NULL
        }
      )

      if (is.null(conda_info)) {
        envs_dir <- infer_conda_envs_dir(conda = conda)
      } else {
        envs_dir <- conda_info[["envs directories"]][1]
        if (is.null(envs_dir) || length(envs_dir) == 0) {
          envs_dir <- conda_info$envs_dirs[1]
        }
      }
      if (is.null(envs_dir) || length(envs_dir) == 0 || !dir.exists(envs_dir)) {
        envs_dir <- infer_conda_envs_dir(conda = conda)
      }
    }

    if (is.null(envs_dir) || length(envs_dir) == 0) {
      return(FALSE)
    }

    env_path <- paste0(envs_dir, "/", envname)
    exist <- file.exists(env_path)

    if (exist) {
      conda_meta_path <- file.path(env_path, "conda-meta")
      if (!dir.exists(conda_meta_path)) {
        log_message(
          "Environment directory exists but appears invalid: {.file {env_path}}",
          message_type = "warning"
        )
        return(FALSE)
      }
    }

    return(exist)
  } else {
    return(FALSE)
  }
}

get_envname <- function(envname = NULL) {
  if (is.character(envname)) {
    envname <- envname
  } else {
    envname <- getOption("scop_envname", default = "scop_env")
  }
  return(envname)
}

get_conda_envs_dir <- function(conda = "auto") {
  conda_info <- get_namespace_fun(
    "reticulate",
    "conda_info"
  )(conda = conda)
  envs_dir <- conda_info[["envs directories"]][1]
  if (is.null(envs_dir) || length(envs_dir) == 0) {
    envs_dir <- conda_info$envs_dirs[1]
  }
  if (is.null(envs_dir) || length(envs_dir) == 0 || !dir.exists(envs_dir)) {
    envs_dir <- infer_conda_envs_dir(conda = conda)
  }
  return(envs_dir)
}

infer_conda_envs_dir <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (is.null(conda) || !file.exists(conda)) {
    return(NULL)
  }
  envs_dir <- file.path(dirname(dirname(conda)), "envs")
  if (!dir.exists(envs_dir)) {
    return(NULL)
  }
  normalizePath(envs_dir, mustWork = FALSE)
}

find_conda <- function() {
  conda <- tryCatch(
    reticulate::conda_binary(conda = "auto"),
    error = identity
  )
  conda_exist <- !inherits(conda, "error")
  if (isFALSE(conda_exist)) {
    if (!is.na(Sys.getenv("USER", unset = NA))) {
      miniconda_path <- gsub(
        pattern = "\\$USER",
        replacement = Sys.getenv("USER"),
        reticulate::miniconda_path()
      )
    } else {
      miniconda_path <- reticulate::miniconda_path()
    }
    conda_exist <- get_namespace_fun(
      "reticulate",
      "miniconda_exists"
    )(
      miniconda_path
    ) &&
      get_namespace_fun(
        "reticulate",
        "miniconda_test"
      )(miniconda_path)
    if (isTRUE(conda_exist)) {
      conda <- get_namespace_fun(
        "reticulate",
        "miniconda_conda"
      )(miniconda_path)
    } else {
      conda <- NULL
    }
  }
  return(conda)
}

resolve_conda <- function(conda = "auto") {
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  return(conda)
}

ensure_conda <- function(conda, error_if_missing = TRUE) {
  if (is.null(conda)) {
    if (error_if_missing) {
      log_message("Conda not found", message_type = "error")
    }
    return(FALSE)
  }
  return(TRUE)
}

get_conda_version <- function(conda) {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = FALSE)) {
    return(NULL)
  }

  system2t <- get_namespace_fun("reticulate", "system2t")
  version_output <- tryCatch(
    {
      system2t(conda, "--version", stdout = TRUE, stderr = TRUE)
    },
    error = function(e) NULL
  )

  if (is.null(version_output) || length(version_output) == 0) {
    return(NULL)
  }

  version_line <- paste(version_output, collapse = " ")
  version_match <- regexec(
    "\\bconda\\s+([0-9]+(?:\\.[0-9]+){1,})\\b",
    version_line,
    perl = TRUE
  )
  match_parts <- regmatches(version_line, version_match)[[1]]
  if (length(match_parts) < 2) {
    return(NULL)
  }

  tryCatch(
    {
      utils::packageVersion(match_parts[2])
    },
    error = function(e) NULL
  )
}

conda_supports_tos <- function(conda) {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = FALSE)) {
    return(FALSE)
  }

  system2t <- get_namespace_fun("reticulate", "system2t")
  status <- tryCatch(
    {
      system2t(
        conda,
        c("tos", "--help"),
        stdout = FALSE,
        stderr = FALSE
      )
    },
    error = function(e) NULL
  )

  identical(status, 0L)
}

install_uv <- function(
  python = NULL,
  envname = NULL,
  conda = "auto",
  pip_options = character()
) {
  log_message("Attempting to install uv...")
  pip_options <- normalize_cli_args(pip_options)

  python <- resolve_python_executable(
    python = python,
    envname = envname,
    conda = conda,
    error_message = "Python not found, cannot install uv"
  )
  if (is.null(python)) {
    return(NULL)
  }

  system2t <- get_namespace_fun("reticulate", "system2t")

  log_message("Installing {.pkg uv} using {.pkg pip}...")

  install_cmd <- c("-m", "pip", "install", "--upgrade")
  if (is.null(envname) && !is.null(python)) {
    if (!grepl("conda|envs", python)) {
      install_cmd <- c(install_cmd, "--user")
    }
  }
  if (length(pip_options) > 0) {
    install_cmd <- c(install_cmd, pip_options)
  }
  install_cmd <- c(install_cmd, "uv")
  install_status <- system2t(python, shQuote(install_cmd))

  if (install_status == 0L) {
    log_message(
      "{.pkg uv} installed successfully via {.pkg pip}",
      message_type = "success"
    )

    python_dir <- dirname(python)
    uv_path <- file.path(python_dir, "uv")
    if (is_windows()) {
      uv_path <- paste0(uv_path, ".exe")
    }

    if (file.exists(uv_path)) {
      return(uv_path)
    }

    uv_system <- Sys.which("uv")
    if (uv_system != "" && file.exists(uv_system)) {
      return(uv_system)
    }

    test_cmd <- c("-m", "uv", "--version")
    test_status <- system2t(python, shQuote(test_cmd))
    if (test_status == 0L) {
      return("python -m uv")
    }

    log_message(
      "{.pkg uv} installed but not found in expected locations",
      message_type = "warning"
    )
    return(NULL)
  } else {
    log_message(
      "Failed to install {.pkg uv} [error code {.val {install_status}}]",
      message_type = "warning"
    )
    return(NULL)
  }
}

find_uv <- function(
  python = NULL,
  envname = NULL,
  conda = "auto",
  auto_install = TRUE,
  pip_options = character()
) {
  pip_options <- normalize_cli_args(pip_options)
  uv <- Sys.which("uv")
  if (uv != "" && file.exists(uv)) {
    return(uv)
  }

  if (is.null(python)) {
    python <- resolve_python_executable(
      envname = envname,
      conda = conda,
      error_if_missing = FALSE
    )
  }

  if (!is.null(python) && file.exists(python)) {
    python_dir <- dirname(python)
    uv_path <- file.path(python_dir, "uv")
    if (is_windows()) {
      uv_path <- paste0(uv_path, ".exe")
    }
    if (file.exists(uv_path)) {
      return(uv_path)
    }

    system2t <- get_namespace_fun("reticulate", "system2t")
    test_cmd <- c("-m", "uv", "--version")
    test_status <- system2t(python, shQuote(test_cmd))
    if (test_status == 0L) {
      return("python -m uv")
    }
  }

  if (isTRUE(auto_install)) {
    uv <- install_uv(
      python = python,
      envname = envname,
      conda = conda,
      pip_options = pip_options
    )
    if (!is.null(uv)) {
      return(uv)
    }
  }

  return(NULL)
}

conda_install <- function(
  envname = NULL,
  packages,
  forge = TRUE,
  channel = character(),
  pip = FALSE,
  pip_options = character(),
  pip_ignore_installed = FALSE,
  conda = "auto",
  python_version = NULL,
  ...
) {
  envname <- get_envname(envname)
  get_namespace_fun(
    "reticulate",
    "check_forbidden_install"
  )("Python packages")
  conda_args <- get_namespace_fun(
    "reticulate",
    "conda_args"
  )
  system2t <- get_namespace_fun(
    "reticulate",
    "system2t"
  )

  if (missing(packages)) {
    if (!is.null(envname)) {
      fmt <- paste(
        "argument \"packages\" is missing, with no default",
        "- did you mean 'conda_install(<envname>, %1$s)'?",
        "- use 'py_install(%1$s)' to install into the active Python environment",
        sep = "\n"
      )
      log_message(
        sprintf(fmt, deparse1(substitute(envname))),
        message_type = "error"
      )
    } else {
      log_message(
        "{.arg packages} argument is required",
        message_type = "error"
      )
    }
  }

  conda <- reticulate::conda_binary(conda)
  envname <- get_namespace_fun(
    "reticulate",
    "condaenv_resolve"
  )(envname)

  log_message(
    "Installing {.val {length(packages)}} packages into environment: {.file {envname}}"
  )

  python_package <- if (is.null(python_version)) {
    NULL
  } else if (grepl("[><=]", python_version)) {
    paste0("python", python_version)
  } else {
    sprintf("python=%s", python_version)
  }

  python <- tryCatch(
    conda_python(envname = envname, conda = conda),
    error = identity
  )

  if (inherits(python, "error") || !file.exists(python)) {
    log_message(
      "Python environment does not exist, creating: {.file {envname}}"
    )
    accept_conda_tos(conda = conda)
    reticulate::conda_create(
      envname = envname,
      packages = python_package %||% "python",
      forge = forge,
      channel = channel,
      conda = conda
    )
    python <- conda_python(envname = envname, conda = conda)
    log_message(
      "Environment created with python: {.pkg {python}}"
    )
  }

  if (!is.null(python_version)) {
    log_message(
      "Updating python version to: {.pkg {python_version}}"
    )
    args <- conda_args("install", envname, python_package)
    status <- system2t(conda, shQuote(args))
    if (status != 0L) {
      fmt <- "installation of '%s' into environment '%s' failed [error code %i]"
      msg <- sprintf(fmt, python_package, envname, status)
      log_message(msg, message_type = "error")
    } else {
      log_message(
        "Python version updated successfully",
        message_type = "success"
      )
    }
  }

  if (pip) {
    install_python_packages(
      packages = packages,
      python = python,
      envname = envname,
      conda = conda,
      pip_options = pip_options,
      pip_ignore_installed = pip_ignore_installed
    )

    return(invisible(packages))
  }

  log_message("Installing packages via {.pkg conda}...")
  args <- conda_args("install", envname)

  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }

  for (ch in channels) {
    args <- c(args, "-c", ch)
    log_message("Using channel: {.val {ch}}")
  }

  conda_packages <- if (is.character(packages)) {
    gsub("==", "=", packages, fixed = TRUE)
  } else {
    packages
  }

  args <- c(args, python_package, conda_packages)

  log_message("Installing {.val {length(packages)}} packages...")

  result <- system2t(conda, shQuote(args))

  if (result != 0L) {
    log_message(
      "Conda installation failed with error code: {.val {result}}",
      message_type = "warning"
    )
  } else {
    log_message(
      "{.pkg conda} installation completed successfully",
      message_type = "success"
    )
  }

  invisible(packages)
}

resolve_python_executable <- function(
  python = NULL,
  envname = NULL,
  conda = "auto",
  error_if_missing = TRUE,
  error_message = NULL
) {
  if (!is.null(python)) {
    python <- normalizePath(python, mustWork = FALSE)
    if (file.exists(python)) {
      return(python)
    }
  }

  if (!is.null(envname)) {
    python <- tryCatch(
      conda_python(envname = envname, conda = conda),
      error = function(e) {
        if (error_if_missing) {
          log_message(
            "Failed to get Python path for environment {.file {envname}}: {.val {e$message}}",
            message_type = "warning"
          )
        }
        NULL
      }
    )
    if (!is.null(python) && file.exists(python)) {
      return(normalizePath(python, mustWork = FALSE))
    }
  }

  for (candidate in c(Sys.which("python3"), Sys.which("python"))) {
    if (nzchar(candidate) && file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }

  if (error_if_missing) {
    log_message(
      error_message %||% "Python executable not found",
      message_type = "error"
    )
  }

  NULL
}

run_uv_command <- function(uv, python, args) {
  system2t <- get_namespace_fun("reticulate", "system2t")

  if (uv == "python -m uv") {
    return(system2t(python, shQuote(c("-m", "uv", args))))
  }

  system2t(uv, shQuote(args))
}

run_pip_command <- function(python, args) {
  system2t <- get_namespace_fun("reticulate", "system2t")
  system2t(python, shQuote(c("-m", "pip", args)))
}

install_python_packages <- function(
  packages,
  python,
  envname,
  conda,
  pip_options = character(),
  pip_ignore_installed = FALSE
) {
  pip_options <- normalize_cli_args(pip_options)
  uv <- find_uv(
    python = python,
    envname = envname,
    conda = conda,
    auto_install = TRUE,
    pip_options = pip_options
  )

  if (!is.null(uv)) {
    log_message("Installing packages via {.pkg uv}...")

    args <- c("pip", "install")
    if (uv != "python -m uv") {
      args <- c(args, "--python", python)
    }
    if (length(pip_options) > 0) {
      args <- c(args, pip_options)
    }
    if (pip_ignore_installed) {
      args <- c(args, "--reinstall")
    }
    args <- c(args, packages)

    status <- run_uv_command(uv, python, args)
    if (status == 0L) {
      log_message(
        "{.pkg uv} installation completed",
        message_type = "success"
      )
      return(invisible(TRUE))
    }

    log_message(
      "{.pkg uv} installation failed [error code {.val {status}}], using {.pkg pip}",
      message_type = "warning"
    )
  } else {
    log_message(
      "{.pkg uv} not available, using {.pkg pip}...",
      message_type = "warning"
    )
  }

  get_namespace_fun("reticulate", "pip_install")(
    python = python,
    packages = packages,
    pip_options = pip_options,
    ignore_installed = pip_ignore_installed,
    conda = conda,
    envname = envname
  )

  invisible(TRUE)
}

ensure_windows_scvi_support <- function(
  envname,
  conda,
  force = FALSE,
  pip_options = character()
) {
  if (!is_windows()) {
    return(invisible(FALSE))
  }

  if (
    !isTRUE(force) &&
      isTRUE(exist_python_pkgs("jax==0.3.20", envname = envname, conda = conda))
  ) {
    return(invisible(TRUE))
  }

  python <- resolve_python_executable(
    envname = envname,
    conda = conda,
    error_if_missing = FALSE
  )
  if (is.null(python)) {
    return(invisible(FALSE))
  }

  pip_options <- normalize_cli_args(pip_options)
  uv <- find_uv(
    python = python,
    envname = envname,
    conda = conda,
    auto_install = TRUE,
    pip_options = pip_options
  )

  args <- c("pip", "install")
  if (!is.null(uv) && uv != "python -m uv") {
    args <- c(args, "--python", python)
  }
  args <- c(
    args,
    pip_options,
    if (isTRUE(force)) "--reinstall",
    "jax[cpu]==0.3.20",
    "-f",
    "https://whls.blob.core.windows.net/unstable/index.html",
    "--use-deprecated",
    "legacy-resolver"
  )

  if (!is.null(uv)) {
    status <- run_uv_command(uv, python, args)
  } else {
    status <- run_pip_command(
      python,
      c(
        "install",
        pip_options,
        if (isTRUE(force)) "--force-reinstall",
        "jax[cpu]==0.3.20",
        "-f",
        "https://whls.blob.core.windows.net/unstable/index.html",
        "--use-deprecated",
        "legacy-resolver"
      )
    )
  }

  if (identical(status, 0L)) {
    log_message(
      "{.pkg jax[cpu]} installed successfully for Windows {.pkg scvi-tools} support",
      message_type = "success"
    )
  } else {
    log_message(
      "Failed to install {.pkg jax[cpu]} for Windows {.pkg scvi-tools} support [error code {.val {status}}]",
      message_type = "warning"
    )
  }

  invisible(TRUE)
}

conda_python <- function(
  envname = NULL,
  conda = "auto",
  all = FALSE
) {
  envname <- get_envname(envname)
  envname <- get_namespace_fun(
    "reticulate",
    "python_environment_resolve"
  )(envname)
  if (grepl("[/\\\\]", envname)) {
    suffix <- if (is_windows()) "python.exe" else "bin/python"
    path <- file.path(envname, suffix)
    if (file.exists(path)) {
      return(path)
    }
    log_message(
      "no conda environment exists at path {.file {envname}}",
      message_type = "error"
    )
  }
  conda_envs <- reticulate::conda_list(conda = conda)
  envs_dir <- get_conda_envs_dir(conda = conda)
  conda_envs <- conda_envs[
    grep(
      normalizePath(
        envs_dir,
        mustWork = FALSE
      ),
      x = normalizePath(conda_envs$python, mustWork = FALSE),
      fixed = TRUE
    ), ,
    drop = FALSE
  ]
  env <- conda_envs[conda_envs$name == envname, , drop = FALSE]
  if (nrow(env) == 0) {
    env_path <- file.path(envs_dir, envname)
    suffix <- if (is_windows()) "python.exe" else "bin/python"
    python_path <- file.path(env_path, suffix)

    if (file.exists(python_path)) {
      return(normalizePath(python_path, mustWork = FALSE))
    }

    log_message(
      "{.val {envname}} environment not found",
      message_type = "error"
    )
  }
  python <- if (all) env$python else env$python[[1L]]
  return(normalizePath(as.character(python), mustWork = FALSE))
}

#' @title Remove a conda environment
#'
#' @md
#' @inheritParams PrepareEnv
#' @param force Whether to force removal without confirmation.
#' Default is `FALSE`.
#'
#' @return Invisibly returns `TRUE` if successful, `FALSE` otherwise.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove default environment
#' RemoveEnv()
#'
#' # Remove a specific environment
#' RemoveEnv("my_old_env")
#'
#' # Removal without confirmation
#' RemoveEnv("my_old_env", force = TRUE)
#' }
RemoveEnv <- function(
  envname = NULL,
  conda = "auto",
  force = FALSE
) {
  envname <- get_envname(envname)

  log_message(
    "Removing conda environment: {.file {envname}}"
  )

  conda <- resolve_conda(conda)
  if (!ensure_conda(conda)) {
    return(invisible(FALSE))
  }

  env_exists <- env_exist(envname = envname, conda = conda)
  if (isFALSE(env_exists)) {
    log_message(
      "{.file {envname}} environment does not exist",
      message_type = "warning"
    )
    return(invisible(FALSE))
  }

  envs_dir <- get_conda_envs_dir(conda = conda)
  env_path <- file.path(envs_dir, envname)

  if (!force) {
    log_message(
      "Environment path: {.file {env_path}}"
    )
    log_message(
      "This will permanently delete the environment and all its packages"
    )

    if (interactive()) {
      response <- readline(
        paste0(
          "Are you sure you want to remove environment: ",
          envname,
          "? (y/N): "
        )
      )
      if (!tolower(response) %in% c("y", "yes")) {
        log_message(
          "Environment removal cancelled",
          message_type = "warning"
        )
        return(invisible(FALSE))
      }
    } else {
      log_message(
        "Use {.arg force = TRUE} to remove environment in non-interactive mode",
        message_type = "warning"
      )
      return(invisible(FALSE))
    }
  }

  log_message("Removing {.file {envname}} environment...")

  result <- tryCatch(
    {
      reticulate::conda_remove(
        envname = envname,
        packages = NULL,
        conda = conda
      )
      log_message(
        "Environment removed successfully using {.fn reticulate::conda_remove}",
        message_type = "success"
      )
      TRUE
    },
    error = function(e) {
      log_message(
        "{.fn reticulate::conda_remove} failed: {.val {e$message}}",
        message_type = "warning"
      )
      FALSE
    }
  )

  if (!result) {
    log_message(
      "Attempting direct removal of environment directory: {.file {env_path}}"
    )

    result <- tryCatch(
      {
        if (dir.exists(env_path)) {
          unlink(env_path, recursive = TRUE, force = TRUE)

          if (!dir.exists(env_path)) {
            log_message(
              "Environment directory {.file {env_path}} removed successfully",
              message_type = "success"
            )
            TRUE
          } else {
            log_message(
              "Failed to remove environment directory: {.file {env_path}}",
              message_type = "error"
            )
            FALSE
          }
        } else {
          log_message(
            "Environment directory does not exist: {.file {env_path}}",
            message_type = "warning"
          )
          TRUE
        }
      },
      error = function(e) {
        log_message(
          "Direct removal failed: {.val {e$message}}",
          message_type = "error"
        )
        FALSE
      }
    )
  }

  if (result) {
    clear_env_cache()
    Sys.unsetenv("RETICULATE_PYTHON")
    log_message(
      "{.file {envname}} environment removed successfully",
      message_type = "success"
    )
  } else {
    log_message(
      "Failed to remove {.file {envname}} environment",
      message_type = "error"
    )
  }

  return(invisible(result))
}

#' @title List conda environments
#'
#' @md
#' @inheritParams PrepareEnv
#' @return A data frame of conda environments.
#' @export
ListEnv <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = TRUE)) {
    return(invisible(NULL))
  }

  conda_envs <- reticulate::conda_list(conda = conda)

  envs_dir <- get_conda_envs_dir(conda = conda)
  if (!is.null(envs_dir) && dir.exists(envs_dir)) {
    existing_names <- conda_envs$name

    env_dirs <- list.dirs(envs_dir, full.names = FALSE, recursive = FALSE)

    for (env_dir in env_dirs) {
      if (!env_dir %in% existing_names) {
        env_path <- file.path(envs_dir, env_dir)
        python_path <- file.path(
          env_path,
          if (is_windows()) "python.exe" else "bin/python"
        )

        if (
          dir.exists(file.path(env_path, "conda-meta")) &&
            file.exists(python_path)
        ) {
          new_row <- data.frame(
            name = env_dir,
            python = normalizePath(python_path, mustWork = FALSE),
            stringsAsFactors = FALSE
          )
          conda_envs <- rbind(conda_envs, new_row)
        }
      }
    }
  }

  if (nrow(conda_envs) > 0) {
    conda_envs <- conda_envs[order(conda_envs$name), , drop = FALSE]
  }

  return(conda_envs)
}

accept_conda_tos <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = FALSE)) {
    return(invisible(FALSE))
  }

  system2t <- get_namespace_fun("reticulate", "system2t")
  conda_version <- get_conda_version(conda)
  supports_tos <- conda_supports_tos(conda)

  if (!supports_tos) {
    version_text <- if (is.null(conda_version)) {
      "unknown"
    } else {
      as.character(conda_version)
    }
    log_message(
      "Skipping conda ToS acceptance because {.pkg conda} does not support the {.val tos} command (version: {.val {version_text}}).",
      message_type = "warning"
    )
    return(invisible(FALSE))
  }

  channels <- c(
    "https://repo.anaconda.com/pkgs/main",
    "https://repo.anaconda.com/pkgs/r",
    "https://repo.anaconda.com/pkgs/msys2"
  )

  log_message(
    "Accepting conda Terms of Service for required channels..."
  )

  for (channel in channels) {
    tryCatch(
      {
        args <- c("tos", "accept", "--override-channels", "--channel", channel)
        status <- system2t(conda, args, stdout = FALSE, stderr = FALSE)
        if (status == 0L) {
          log_message(
            "Accepted ToS for channel: {.val {channel}}",
            message_type = "success"
          )
        } else {
          log_message(
            "Failed to accept ToS for channel: {.val {channel}} [error code {.val {status}}]",
            message_type = "warning"
          )
        }
      },
      error = function(e) {
        log_message(
          "Failed to accept ToS for channel: {.val {channel}} ({.val {e$message}})",
          message_type = "warning"
        )
      }
    )
  }

  return(invisible(TRUE))
}
