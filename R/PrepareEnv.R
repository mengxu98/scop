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
#' Default is `"3.10-1"`.
#' @param pip_options Additional command line arguments to be passed to `uv`/`pip` when installing pip packages.
#' @param ... Additional arguments passed to package installation functions.
#'
#' @export
PrepareEnv <- function(
    envname = NULL,
    conda = "auto",
    miniconda_repo = "https://repo.anaconda.com/miniconda",
    version = "3.10-1",
    force = FALSE,
    pip_options = character(),
    ...) {
  env_cache <- getOption("scop_env_cache", default = NULL)
  if (isTRUE(env_cache) && isFALSE(force)) {
    log_message(
      "{cli::col_green('Python environment already prepared')}\n",
      "{cli::col_grey('Until next loading, the environment will be cached')}",
      message_type = "success"
    )
    return(invisible(NULL))
  }

  log_message(
    "Preparing python environment...",
    text_color = "orange",
    message_type = "running",
    timestamp_style = FALSE
  )

  if (!is.null(envname)) {
    options(scop_envname = envname)
  }

  envname <- get_envname(envname)
  requirements <- env_requirements(version = version)
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
      ...
    )
  }

  set_python_env(conda = conda, envname = envname)

  log_message(
    "{cli::col_green('Python environment ready')}\n",
    "{cli::col_grey('Until next loading, the environment will be cached')}",
    message_type = "success"
  )

  env_info(conda, envname)

  options(scop_env_cache = TRUE)
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
  reticulate::use_python(
    python_path,
    required = TRUE
  )

  tryCatch(
    {
      reticulate::py_run_string("
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
")
    },
    error = function(e) {
    }
  )
}

install_miniconda2 <- function(
    miniconda_repo,
    timeout = 600) {
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
      "reticulate", "miniconda_installer_arch"
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
    "reticulate", "miniconda_conda"
  )(miniconda_path)
  log_message("Miniconda installed at: {.file {miniconda_path}}")

  conda
}

#' @title Print environment information
#' @inheritParams PrepareEnv
#' @md
env_info <- function(conda, envname) {
  envs_dir <- get_conda_envs_dir(conda = conda)

  py_info <- utils::capture.output(reticulate::py_config())

  py_info_mesg <- c(
    cli::col_blue(
      "Conda config:"
    ),
    cli::col_grey(
      " Conda:         ", conda
    ),
    cli::col_grey(
      " Environment:   ", envs_dir, "/", get_envname(envname)
    ),
    cli::col_blue(
      "Python config:"
    ),
    cli::col_grey(
      " ", py_info
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
env_requirements <- function(version = "3.10-1") {
  version <- match.arg(
    version,
    choices = c("3.8-1", "3.8-2", "3.9-1", "3.10-1", "3.11-1")
  )

  tbb_install_method <- if (is_apple_silicon()) {
    "pip"
  } else {
    "conda"
  }

  package_install_methods <- c(
    "leidenalg" = "conda",
    "tbb" = tbb_install_method,
    "python-igraph" = "conda",
    "scvi-tools" = "conda",
    "matplotlib" = "pip",
    "numba" = "pip",
    "llvmlite" = "pip",
    "numpy" = "pip",
    "packaging" = "pip",
    "palantir" = "pip",
    "pandas" = "pip",
    "scanpy" = "pip",
    "scikit-learn" = "pip",
    "scipy" = "pip",
    "scvelo" = "pip",
    "wot" = "pip",
    "trimap" = "pip",
    "pacmap" = "pip",
    "phate" = "pip",
    "bbknn" = "pip",
    "scanorama" = "pip",
    "cellrank" = "pip",
    "celltypist" = "pip"
  )

  package_versions <- c(
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
    "palantir" = "palantir==1.4.1",
    "pandas" = "pandas==2.0.3",
    "scanpy" = "scanpy==1.11.3",
    "scikit-learn" = "scikit-learn==1.7.0",
    "scipy" = "scipy==1.15.3",
    "scvelo" = "scvelo==0.3.3",
    "wot" = "wot==1.0.8.post2",
    "trimap" = "trimap==1.1.4",
    "pacmap" = "pacmap==0.8.0",
    "phate" = "phate==1.0.11",
    "bbknn" = "bbknn==1.6.0",
    "scanorama" = "scanorama==1.7.4",
    "scvi-tools" = "scvi-tools==1.2.1",
    "cellrank" = "cellrank==2.0.7",
    "celltypist" = "celltypist"
  )

  package_aliases <- list(
    "python-igraph" = "igraph"
  )

  if (is_apple_silicon()) {
    package_aliases[["tbb"]] <- "pxr-tbb"
  }

  requirements <- list(
    python = version,
    packages = package_versions,
    install_methods = package_install_methods,
    package_aliases = package_aliases
  )

  return(requirements)
}

env_exist <- function(
    conda = "auto",
    envname = NULL,
    envs_dir = NULL) {
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
        return(FALSE)
      }

      envs_dir <- conda_info[["envs directories"]][1]
      if (is.null(envs_dir) || length(envs_dir) == 0) {
        envs_dir <- conda_info$envs_dirs[1]
      }
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
    "reticulate", "conda_info"
  )(conda = conda)
  envs_dir <- conda_info[["envs directories"]][1]
  if (is.null(envs_dir) || length(envs_dir) == 0) {
    envs_dir <- conda_info$envs_dirs[1]
  }
  return(envs_dir)
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
      "reticulate", "miniconda_exists"
    )(
      miniconda_path
    ) &&
      get_namespace_fun(
        "reticulate", "miniconda_test"
      )(miniconda_path)
    if (isTRUE(conda_exist)) {
      conda <- get_namespace_fun(
        "reticulate", "miniconda_conda"
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

install_uv <- function(
    python = NULL,
    envname = NULL,
    conda = "auto",
    pip_options = character()) {
  log_message("Attempting to install uv...")

  if (is.null(python)) {
    if (!is.null(envname)) {
      tryCatch(
        {
          python <- conda_python(envname = envname, conda = conda)
        },
        error = function(e) {
          log_message(
            "Failed to get Python path for environment {.file {envname}}: {.val {e$message}}",
            message_type = "warning"
          )
        }
      )
    }

    if (is.null(python)) {
      python <- Sys.which("python3")
      if (python == "" || !file.exists(python)) {
        python <- Sys.which("python")
      }
      if (python == "" || !file.exists(python)) {
        log_message(
          "Python not found, cannot install uv",
          message_type = "error"
        )
        return(NULL)
      }
    }
  }

  if (!file.exists(python)) {
    log_message(
      "Python executable not found: {.file {python}}",
      message_type = "error"
    )
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
    # Split pip_options if it's a single string
    if (length(pip_options) == 1 && grepl(" ", pip_options)) {
      pip_options <- strsplit(pip_options, "\\s+")[[1]]
    }
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
    pip_options = character()) {
  uv <- Sys.which("uv")
  if (uv != "" && file.exists(uv)) {
    return(uv)
  }

  if (is.null(python)) {
    if (!is.null(envname)) {
      tryCatch(
        {
          python <- conda_python(envname = envname, conda = conda)
        },
        error = function(e) NULL
      )
    }
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
    uv <- install_uv(python = python, envname = envname, conda = conda, pip_options = pip_options)
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
    ...) {
  envname <- get_envname(envname)
  get_namespace_fun(
    "reticulate", "check_forbidden_install"
  )("Python packages")
  conda_args <- get_namespace_fun(
    "reticulate", "conda_args"
  )
  system2t <- get_namespace_fun(
    "reticulate", "system2t"
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
    "reticulate", "condaenv_resolve"
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
    uv <- find_uv(
      python = python,
      envname = envname,
      conda = conda,
      auto_install = TRUE,
      pip_options = pip_options
    )

    if (!is.null(uv)) {
      log_message("Installing packages via {.pkg uv}...")

      args <- c("pip", "install", "--python", python)

      if (length(pip_options) > 0) {
        # Split pip_options if it's a single string
        if (length(pip_options) == 1 && grepl(" ", pip_options)) {
          pip_options <- strsplit(pip_options, "\\s+")[[1]]
        }
        args <- c(args, pip_options)
      }

      if (pip_ignore_installed) {
        args <- c(args, "--reinstall")
      }

      args <- c(args, packages)

      if (uv == "python -m uv") {
        args <- c("-m", "uv", args)
        status <- system2t(python, shQuote(args))
      } else {
        status <- system2t(uv, shQuote(args))
      }

      if (status != 0L) {
        log_message(
          "{.pkg uv} installation failed [error code {.val {status}}], using {.pkg pip}",
          message_type = "warning"
        )
        get_namespace_fun("reticulate", "pip_install")(
          python = python,
          packages = packages,
          pip_options = pip_options,
          ignore_installed = pip_ignore_installed,
          conda = conda,
          envname = envname
        )
      } else {
        log_message(
          "{.pkg uv} installation completed",
          message_type = "success"
        )
      }
    } else {
      log_message(
        "{.pkg uv} not available, using {.pkg pip}...",
        message_type = "warning"
      )
      get_namespace_fun("reticulate", "pip_install")(
        python = python,
        packages = packages,
        pip_options = pip_options,
        ignore_installed = pip_ignore_installed,
        conda = conda,
        envname = envname
      )
    }

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

conda_python <- function(
    envname = NULL,
    conda = "auto",
    all = FALSE) {
  envname <- get_envname(envname)
  envname <- get_namespace_fun(
    "reticulate", "python_environment_resolve"
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
    force = FALSE) {
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
          "Are you sure you want to remove environment: ", envname, "? (y/N): "
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
          env_path, if (is_windows()) "python.exe" else "bin/python"
        )

        if (dir.exists(file.path(env_path, "conda-meta")) && file.exists(python_path)) {
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

  version_cmd <- c("--version")
  version_output <- tryCatch(
    {
      system2t(conda, shQuote(version_cmd), stdout = TRUE, stderr = TRUE)
    },
    error = function(e) NULL
  )

  if (!is.null(version_output)) {
    conda_version <- tryCatch(
      {
        version_str <- gsub(".*conda ([0-9]+)\\.[0-9]+.*", "\\1", version_output)
        as.numeric(version_str)
      },
      error = function(e) NULL
    )

    if (!is.null(conda_version) && conda_version >= 23) {
      return(invisible(TRUE))
    }
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
        status <- system2t(conda, shQuote(args), stdout = FALSE, stderr = FALSE)
        if (status == 0L) {
          log_message(
            "Accepted ToS for channel: {.val {channel}}",
            message_type = "success"
          )
        }
      },
      error = function(e) {
      }
    )
  }

  return(invisible(TRUE))
}

is_linux <- function() {
  identical(tolower(Sys.info()[["sysname"]]), "linux")
}

is_osx <- function() {
  identical(tolower(Sys.info()[["sysname"]]), "darwin")
}

is_apple_silicon <- function() {
  is_osx() &&
    identical(tolower(Sys.info()[["machine"]]), "arm64")
}

is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}
