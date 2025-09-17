#' @title Prepare the virtual environment
#'
#' @description
#' Prepare the virtual environment by installing the required dependencies and setting up the environment.
#' This function prepares the virtual environment by checking if conda is installed,
#' creating a new conda environment if needed, installing the required packages,
#' and setting up the Python environment for use with scop.
#' In order to create the environment, this function requires the path to the conda binary.
#' If `conda` is set to `"auto"`, it will attempt to automatically find the conda binary.
#' If a conda environment with the specified name already exists and `force` is set to `FALSE`,
#' the function will use the existing environment.
#' If `force` set to `TRUE`, the existing environment will be recreated.
#' Note that recreating the environment will remove any existing data in the environment.
#' The function also checks if the package versions in the environment meet the requirements specified by the `version` parameter.
#'
#' @md
#' @inheritParams check_python
#' @param miniconda_repo Repositories for miniconda.
#' Default is \url{https://repo.anaconda.com/miniconda}.
#' @param force Whether to force a new environment to be created.
#' If `TRUE`, the existing environment will be recreated.
#' Default is `FALSE`.
#' @param version A character vector specifying the version of the environment.
#' Default is `"3.10-1"`.
#' @export
#'
#' @examples
#' PrepareEnv()
PrepareEnv <- function(
    conda = "auto",
    miniconda_repo = "https://repo.anaconda.com/miniconda",
    envname = NULL,
    version = "3.10-1",
    force = FALSE,
    ...) {
  log_message(
    "{cli::col_blue('Preparing scop Python Environment')}"
  )
  if (!is.null(envname)) {
    options(scop_envname = envname)
  }

  envname <- get_envname(envname)
  log_message(
    "Environment name: {.file {envname}}"
  )

  requirements <- env_requirements(version = version)
  python_version <- requirements[["python"]]
  packages <- requirements[["packages"]]

  log_message(
    "Python version: {.pkg {python_version}}"
  )
  log_message(
    "Number of packages to install: {.val {length(packages)}}"
  )

  if (!is.null(conda)) {
    if (identical(conda, "auto")) {
      log_message("Auto-detecting conda...")
      conda <- find_conda()
    } else {
      options(reticulate.conda_binary = conda)
      conda <- find_conda()
    }
  }

  if (is.null(conda)) {
    env <- FALSE
    log_message(
      "Conda not found, will install miniconda",
      message_type = "warning"
    )
  } else {
    envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
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

    if (python_version < numeric_version("3.8.0") || python_version >= numeric_version("3.12.0")) {
      log_message(
        "scop currently supports Python versions 3.8-3.12. Requested: {.val {python_version}}",
        message_type = "error"
      )
    }

    log_message(
      "Creating conda environment with Python {.val {python_version}}..."
    )

    python_path <- reticulate::conda_create(
      conda = conda,
      envname = envname,
      python_version = python_version,
      packages = c("pip", "setuptools", "wheel")
    )

    envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    env_path <- paste0(envs_dir, "/", envname)
    env <- file.exists(env_path)

    if (isFALSE(env)) {
      log_message("Environment creation failed!")
      print(reticulate:::conda_info(conda = conda))
      print(reticulate::conda_list(conda = conda))
      log_message(
        "Unable to find scop environment under the expected path: {.file {env_path}}\n",
        "conda: {.pkg {conda}}\n",
        "scop python: {.file {python_path}}",
        message_type = "error"
      )
    } else {
      log_message(
        "Environment created successfully: {.file {env_path}}",
        message_type = "success"
      )
    }
  }

  log_message("Checking and installing required packages...")

  install_methods <- requirements[["install_methods"]]

  conda_packages <- packages[install_methods == "conda"]
  pip_packages <- packages[install_methods == "pip"]

  if (length(conda_packages) > 0) {
    log_message(
      "Installing {.pkg conda} packages"
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
    log_message(
      "Installing {.pkg pip} packages"
    )
    check_python(
      packages = pip_packages,
      envname = envname,
      conda = conda,
      force = force,
      pip = TRUE,
      ...
    )
  }
  set_python_env(conda = conda, envname = envname)
  log_message(
    "{cli::col_green('Python Environment Ready')}",
    message_type = "success"
  )
  env_info(conda, envname)
}

set_python_env <- function(conda, envname, verbose = TRUE) {
  Sys.unsetenv("RETICULATE_PYTHON")

  options(reticulate.keras.backend = "tensorflow")
  options(reticulate.miniconda.enabled = FALSE)

  python_path <- conda_python(conda = conda, envname = envname)
  reticulate::use_python(python_path, required = TRUE)

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
")

      if (Sys.info()["sysname"] == "Darwin" && Sys.info()["machine"] == "arm64") {
        Sys.setenv(NUMBA_NUM_THREADS = "1")
        Sys.setenv(NUMBA_DISABLE_JIT = "1")
        Sys.setenv(NUMBA_THREADING_LAYER = "tbb")
        Sys.setenv(NUMBA_DEFAULT_NUM_THREADS = "1")
        Sys.setenv(OMP_NUM_THREADS = "1")
        Sys.setenv(OPENBLAS_NUM_THREADS = "1")
        Sys.setenv(MKL_NUM_THREADS = "1")
        Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
        Sys.setenv(NUMEXPR_NUM_THREADS = "1")
        Sys.setenv(NUMBA_CACHE_DIR = "/tmp/numba_cache")
        Sys.setenv(NUMBA_DEBUG = "0")
      }
    },
    error = function(e) {
      cli::col_red(
        "Could not set python environment variables"
      )
    }
  )
}

#' Enhanced miniconda installation
#' @param miniconda_repo Repository URL for miniconda
#' @export
install_miniconda2 <- function(miniconda_repo) {
  log_message("Installing miniconda...")
  options(timeout = 600)

  info <- as.list(Sys.info())

  if (info$sysname == "Darwin" && info$machine == "arm64") {
    base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
    name <- "Miniforge3-MacOSX-arm64.sh"
    url <- file.path(base, name)
  } else {
    version <- "3"
    arch <- reticulate:::miniconda_installer_arch(info)

    name <- if (.is_windows()) {
      sprintf("Miniconda%s-latest-Windows-%s.exe", version, arch)
    } else if (.is_osx()) {
      sprintf("Miniconda%s-latest-MacOSX-%s.sh", version, arch)
    } else if (.is_linux()) {
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

  conda <- reticulate:::miniconda_conda(miniconda_path)
  log_message("Miniconda installed at: {.file {miniconda_path}}")

  conda
}

#' Print environment information
#' @param conda Conda binary path
#' @param envname Environment name
env_info <- function(conda, envname) {
  envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]

  py_info <- utils::capture.output(reticulate::py_config())

  py_info_mesg <- c(
    cli::col_blue(
      "conda environment: "
    ),
    cli::col_grey(
      "  conda:          ", conda
    ),
    cli::col_grey(
      "  environment:    ", envs_dir, "/", get_envname()
    ),
    cli::col_blue(
      "python config: "
    ),
    cli::col_grey(
      "  ", py_info
    )
  )
  invisible(lapply(py_info_mesg, packageStartupMessage))
}

#' @title Python environment requirements
#'
#' @description
#' The function returns a list of requirements including the required Python version
#' and a list of packages with their corresponding versions.
#'
#' @param version A character vector specifying the version of the environment.
#' Default is "3.10-1".
#' @return A list of requirements for the specified version.
#'
#' @export
#' @examples
#' env_requirements("3.10-1")
env_requirements <- function(version = "3.10-1") {
  version <- match.arg(
    version,
    choices = c("3.8-1", "3.8-2", "3.9-1", "3.10-1", "3.11-1")
  )
  package_install_methods <- c(
    "leidenalg" = "conda",
    "tbb" = "conda",
    "python-igraph" = "conda",
    "matplotlib" = "pip",
    "numba" = "pip",
    "llvmlite" = "pip",
    "numpy" = "pip",
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
    "scvi-tools" = "pip",
    "cellrank" = "pip"
  )

  package_versions <- c(
    "leidenalg" = "leidenalg==0.10.2",
    "tbb" = "tbb==2022.2.0",
    "python-igraph" = "python-igraph==0.11.9",
    "matplotlib" = "matplotlib==3.10.3",
    "numba" = "numba==0.59.1",
    "llvmlite" = "llvmlite==0.42.0",
    "numpy" = "numpy==1.26.4",
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
    "cellrank" = "cellrank==2.0.7"
  )

  requirements <- list(
    python = version,
    packages = package_versions,
    install_methods = package_install_methods
  )

  return(requirements)
}

#' Show all the python packages in the environment
#'
#' @inheritParams check_python
#' @export
installed_python_pkgs <- function(
    envname = NULL,
    conda = "auto") {
  envname <- get_envname(envname)

  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }

  if (is.null(conda)) {
    log_message("Conda not found", message_type = "error")
  }

  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    log_message(
      "Cannot find the conda environment: {.file {envname}}",
      message_type = "error"
    )
  }

  log_message(
    "Retrieving package list for environment: {.file {envname}}"
  )

  tryCatch(
    {
      all_installed <- reticulate:::conda_list_packages(
        conda = conda,
        envname = envname,
        no_pip = FALSE
      )
      log_message("Found {.val {nrow(all_installed)}} packages installed")
      return(all_installed)
    },
    error = function(e) {
      log_message(
        "Failed to retrieve package list: {.val {e$message}}",
        message_type = "error"
      )
    }
  )
}

#' Check if the python package exists in the environment
#'
#' @inheritParams check_python
#' @export
exist_python_pkgs <- function(
    packages,
    envname = NULL,
    conda = "auto") {
  envname <- get_envname(envname)

  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }

  if (is.null(conda)) {
    log_message("Conda not found", message_type = "error")
  }

  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    log_message(
      "Cannot find the conda environment: {.file {envname}}",
      message_type = "error"
    )
  }

  log_message(
    "Checking {.val {length(packages)}} packages in environment: {.file {envname}}"
  )

  all_installed <- tryCatch(
    {
      installed_python_pkgs(envname = envname, conda = conda)
    },
    error = function(e) {
      log_message(
        "Failed to get installed packages: {.val {e$message}}",
        message_type = "warning"
      )
    }
  )

  packages_installed <- stats::setNames(
    rep(FALSE, length(packages)), packages
  )

  for (i in seq_along(packages)) {
    pkg <- packages[i]

    if (grepl("==", pkg)) {
      pkg_info <- strsplit(pkg, split = "==")[[1]]
      pkg_name <- names(pkg) %||% pkg_info[1]
      pkg_version <- pkg_info[2]
    } else if (grepl("git+", pkg)) {
      pkg_info <- strsplit(pkg, "/")[[1]]
      pkg_name <- names(pkg) %||% pkg_info[length(pkg_info)]
      pkg_version <- NA
    } else {
      pkg_name <- names(pkg) %||% pkg
      pkg_version <- NA
    }

    if (pkg_name %in% all_installed$package) {
      if (!is.na(pkg_version)) {
        installed_version <- all_installed$version[all_installed$package == pkg_name]
        packages_installed[pkg] <- installed_version == pkg_version
        if (packages_installed[pkg]) {
          log_message(
            "{.pkg {pkg_name}} {.pkg {pkg_version}}",
            message_type = "success"
          )
        } else {
          log_message(
            "{.pkg {pkg_name}} found but version mismatch: installed {.pkg {installed_version}}, required {.pkg {pkg_version}}",
            message_type = "warning"
          )
        }
      } else {
        packages_installed[pkg] <- TRUE
        installed_version <- all_installed$version[all_installed$package == pkg_name]
        log_message(
          "{.pkg {pkg_name}} version: {.pkg {installed_version}}",
          message_type = "success"
        )
      }
    } else {
      packages_installed[pkg] <- FALSE
      log_message(
        "{.pkg {pkg_name}} not found",
        message_type = "warning"
      )
    }
  }

  return(packages_installed)
}

#' Check if a conda environment exists
#'
#' @param envs_dir Directories in which conda environments are located.
#' @inheritParams check_python
env_exist <- function(
    conda = "auto",
    envname = NULL,
    envs_dir = NULL) {
  envname <- get_envname(envname)

  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }

  if (!is.null(conda)) {
    if (is.null(envs_dir)) {
      conda_info <- tryCatch(
        {
          reticulate:::conda_info(conda = conda)
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

      envs_dir <- conda_info$envs_dirs[1]
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

#' Find an appropriate conda binary
#'
#' @export
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
    conda_exist <- reticulate:::miniconda_exists(
      miniconda_path
    ) &&
      reticulate:::miniconda_test(miniconda_path)
    if (isTRUE(conda_exist)) {
      conda <- reticulate:::miniconda_conda(miniconda_path)
    } else {
      conda <- NULL
    }
  }
  return(conda)
}

#' Enhanced conda installation
#'
#' @inheritParams reticulate::conda_install
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
  reticulate:::check_forbidden_install("Python packages")

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
  envname <- reticulate:::condaenv_resolve(envname)

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
    log_message("Environment doesn't exist, creating: {.file {envname}}")
    reticulate::conda_create(
      envname = envname,
      packages = python_package %||% "python",
      forge = forge,
      channel = channel,
      conda = conda
    )
    python <- conda_python(envname = envname, conda = conda)
    log_message("Environment created with Python: {.pkg {python}}")
  }

  if (!is.null(python_version)) {
    log_message("Updating Python to version: {.pkg {python_version}}")
    args <- reticulate:::conda_args("install", envname, python_package)
    status <- reticulate:::system2t(conda, shQuote(args))
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
    log_message("Installing packages via {.pkg pip}...")
    result <- tryCatch(
      {
        reticulate:::pip_install(
          python = python,
          packages = packages,
          pip_options = pip_options,
          ignore_installed = pip_ignore_installed,
          conda = conda,
          envname = envname
        )
      },
      error = function(e) {
        log_message(
          "{.pkg pip} installation failed: {.val {e$message}}",
          message_type = "error"
        )
      }
    )

    if (!is.null(result)) {
      log_message(
        "{.pkg pip} installation completed",
        message_type = "success"
      )
    }
    return(result)
  }

  log_message("Installing packages via {.pkg conda}...")
  args <- reticulate:::conda_args("install", envname)

  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }

  for (ch in channels) {
    args <- c(args, "-c", ch)
    log_message("Using channel: {.val {ch}}")
  }

  args <- c(args, python_package, packages)

  log_message("Installing {.val {length(packages)}} packages...")

  result <- reticulate:::system2t(conda, shQuote(args))

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

#' Find the path to Python associated with a conda environment
#'
#' @inheritParams reticulate::conda_python
conda_python <- function(
    envname = NULL,
    conda = "auto",
    all = FALSE) {
  envname <- get_envname(envname)
  envname <- reticulate:::python_environment_resolve(envname)
  if (grepl("[/\\\\]", envname)) {
    suffix <- if (.is_windows()) "python.exe" else "bin/python"
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
  conda_envs <- conda_envs[
    grep(
      normalizePath(
        reticulate:::conda_info(conda = conda)$envs_dirs[1],
        mustWork = FALSE
      ),
      x = normalizePath(conda_envs$python, mustWork = FALSE),
      fixed = TRUE
    ), ,
    drop = FALSE
  ]
  env <- conda_envs[conda_envs$name == envname, , drop = FALSE]
  if (nrow(env) == 0) {
    log_message(
      "{.val {envname}} environment not found",
      message_type = "error"
    )
  }
  python <- if (all) env$python else env$python[[1L]]
  return(normalizePath(as.character(python), mustWork = FALSE))
}

#' Remove a conda environment
#'
#' @param envname The name of the conda environment to remove.
#' If NULL, uses the default scop environment name.
#' @param conda The path to a conda executable. Use "auto" to allow
#' reticulate to automatically find an appropriate conda binary.
#' @param force Whether to force removal without confirmation.
#' Default is FALSE.
#'
#' @details This function removes a conda environment completely.
#' This action cannot be undone, so use with caution.
#' If the environment is currently active, it will be deactivated first.
#'
#' @return Invisibly returns TRUE if successful, FALSE otherwise.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove the default scop environment
#' RemoveEnv()
#'
#' # Remove a specific environment
#' RemoveEnv("my_old_env")
#'
#' # Force removal without confirmation
#' RemoveEnv("my_old_env", force = TRUE)
#' }
RemoveEnv <- function(
    envname = NULL,
    conda = "auto",
    force = FALSE) {
  envname <- get_envname(envname)

  log_message("Removing conda environment: {.val {envname}}")

  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }

  if (is.null(conda)) {
    log_message("Conda not found", message_type = "error")
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

  envs_dir <- reticulate:::conda_info(conda)$envs_dirs[1]
  env_path <- file.path(envs_dir, envname)

  if (!force) {
    log_message("Environment path: {.file {env_path}}")
    log_message("This will permanently delete the environment and all its packages.")

    if (interactive()) {
      response <- readline(
        paste0(
          "Are you sure you want to remove environment: ", envname, "? (y/N): "
        )
      )
      if (!tolower(response) %in% c("y", "yes")) {
        log_message("Environment removal cancelled.")
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
#' @param conda The path to a conda executable.
#' Use `"auto"` to allow reticulate to automatically find an appropriate conda binary.
#'
#' @return A data frame of conda environments.
#' @export
ListEnv <- function(conda = "auto") {
  reticulate::conda_list(conda = conda)
}

.is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

.is_osx <- function() {
  Sys.info()[["sysname"]] == "Darwin"
}

.is_linux <- function() {
  identical(tolower(Sys.info()[["sysname"]]), "linux")
}
