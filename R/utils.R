#' Prepare the scop Python environment by installing the required dependencies and setting up the environment.
#'
#' @param miniconda_repo  Repositories for miniconda. Default is \code{https://repo.anaconda.com/miniconda}
#' @param force Whether to force a new environment to be created. If \code{TRUE}, the existing environment will be recreated. Default is \code{FALSE}
#' @param version A character vector specifying the version of the environment (default is "3.8-1").
#' @inheritParams check_python
#' @details This function prepares the scop Python environment by checking if conda is installed, creating a new conda environment if needed, installing the required packages, and setting up the Python environment for use with scop.
#' In order to create the environment, this function requires the path to the conda binary. If \code{conda} is set to \code{"auto"}, it will attempt to automatically find the conda binary.
#' If a conda environment with the specified name already exists and \code{force} is set to \code{FALSE}, the function will use the existing environment. If \code{force} set to \code{TRUE}, the existing environment will be recreated. Note that recreating the environment will remove any existing data in the environment.
#' The function also checks if the package versions in the environment meet the requirements specified by the \code{version} parameter. The default is \code{3.8-1}.
#'
#' @export
PrepareEnv <- function(
    conda = "auto",
    miniconda_repo = "https://repo.anaconda.com/miniconda",
    envname = NULL,
    version = "3.10-1",
    force = FALSE,
    ...) {
  log_message("=== Preparing scop Python Environment ===")

  envname <- get_envname(envname)
  log_message("Environment name: ", envname)

  requirements <- env_requirements(version = version)
  python_version <- requirements[["python"]]
  packages <- requirements[["packages"]]

  log_message("Python version: ", python_version)
  log_message("Number of packages to install: ", length(packages))

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
    log_message("Conda not found, will install miniconda")
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
      log_message("Using existing environment: ", paste0(envs_dir, "/", envname))
    }
  }

  if (isTRUE(env)) {
    python_path <- conda_python(conda = conda, envname = envname)
    log_message("Python path: ", python_path)
  } else {
    force <- TRUE

    if (is.null(conda)) {
      log_message("Installing miniconda...")
      conda <- install_miniconda_enhanced(miniconda_repo)
    }

    if (python_version < numeric_version("3.8.0") || python_version >= numeric_version("3.12.0")) {
      log_message(
        "scop currently supports Python versions 3.8-3.11. Requested: ", python_version,
        message_type = "error"
      )
    }

    log_message("Creating conda environment with Python ", python_version, "...")

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
        "Unable to find scop environment under the expected path: ", env_path, "\n",
        "conda: ", conda, "\n",
        "scop python: ", python_path,
        message_type = "error"
      )
    } else {
      log_message("Environment created successfully: ", env_path)
    }
  }

  log_message("Checking and installing required packages...")
  check_python(
    packages = packages,
    envname = envname,
    conda = conda,
    force = force,
    ...
  )

  log_message("Setting up Python environment...")
  Sys.unsetenv("RETICULATE_PYTHON")

  options(reticulate.keras.backend = "tensorflow")
  options(reticulate.miniconda.enabled = FALSE)

  python_path <- conda_python(conda = conda, envname = envname)
  reticulate::use_python(python_path, required = TRUE)

  env_info(conda, envname)

  initialize_python_modules()

  log_message("=== scop Python Environment Ready ===")
}

#' Enhanced miniconda installation
#' @param miniconda_repo Repository URL for miniconda
install_miniconda_enhanced <- function(miniconda_repo) {
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

    name <- if (reticulate:::is_windows()) {
      sprintf("Miniconda%s-latest-Windows-%s.exe", version, arch)
    } else if (reticulate:::is_osx()) {
      sprintf("Miniconda%s-latest-MacOSX-%s.sh", version, arch)
    } else if (reticulate:::is_linux()) {
      sprintf("Miniconda%s-latest-Linux-%s.sh", version, arch)
    } else {
      log_message(
        "Unsupported platform: ", info$sysname,
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

  log_message("Downloading and installing miniconda from: ", url)
  reticulate::install_miniconda(
    path = miniconda_path,
    force = TRUE,
    update = FALSE
  )

  conda <- reticulate:::miniconda_conda(miniconda_path)
  log_message("Miniconda installed at: ", miniconda_path)

  return(conda)
}

#' Print environment information
#' @param conda Conda binary path
#' @param envname Environment name
env_info <- function(conda, envname) {
  envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]

  pyinfo <- utils::capture.output(reticulate::py_config())

  pyinfo_mesg <- c(
    "==================== scop conda environment ====================",
    paste0("conda: ", conda),
    paste0("environment: ", paste0(envs_dir, "/", envname)),
    paste0("python: ", conda_python(conda = conda, envname = envname)),
    "==================== scop python config ====================",
    pyinfo,
    "=============================================================="
  )

  invisible(lapply(pyinfo_mesg, packageStartupMessage))
}

initialize_python_modules <- function() {
  log_message("Initializing Python modules...")

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
    },
    error = function(e) {
      log_message(
        "Warning: Could not set Python environment variables",
        message_type = "warning"
      )
    }
  )

  tryCatch(
    {
      reticulate::py_run_string("import matplotlib")
      reticulate::py_run_string("matplotlib.use('Agg')")
      reticulate::py_run_string("import matplotlib.pyplot as plt")
      log_message("matplotlib initialized", message_type = "success")
    },
    error = function(e) {
      log_message(
        "Warning: matplotlib initialization failed: ", e$message,
        message_type = "warning"
      )
    }
  )

  Sys.sleep(0.1)

  tryCatch(
    {
      reticulate::py_run_string("
import numba
numba.config.THREADING_LAYER = 'safe'
")
      reticulate::py_run_string("import scanpy as sc")
      reticulate::py_run_string("sc.settings.verbosity = 1")
      log_message("scanpy initialized", message_type = "success")
    },
    error = function(e) {
      log_message(
        "Warning: scanpy initialization failed: ", e$message,
        message_type = "warning"
      )
    }
  )

  tryCatch(
    {
      reticulate::py_run_string("import scvelo as scv")
      reticulate::py_run_string("scv.settings.verbosity = 1")
      log_message("scvelo initialized", message_type = "success")
    },
    error = function(e) {
      log_message("Note: scvelo not available (this is optional)")
    }
  )
}

#' env_requirements function
#'
#' This function provides the scop python environment requirements for a specific version.
#'
#' @param version A character vector specifying the version of the environment.
#' Default is "3.10-1".
#' @return A list of requirements for the specified version.
#' @details The function returns a list of requirements including the required Python version
#' and a list of packages with their corresponding versions.
#'
#' @export
#' @examples
#' \dontrun{
#' env_requirements("3.10-1")
#' }
env_requirements <- function(version = "3.10-1") {
  version <- match.arg(
    version,
    choices = c("3.8-1", "3.8-2", "3.9-1", "3.10-1", "3.11-1")
  )

  requirements <- list(
    python = version,
    packages = c(
      "leidenalg" = "leidenalg",
      "matplotlib" = "matplotlib>=3.5,<3.11",
      "numba" = "numba>=0.59,<0.60.0",
      "llvmlite" = "llvmlite>=0.42,<0.43.0",
      "numpy" = "numpy>=1.24,<1.27.0",
      "palantir" = "palantir",
      "pandas" = "pandas>=2.0,<2.1",
      "python-igraph" = "python-igraph",
      "scanpy" = "scanpy>=1.9,<1.12",
      "scikit-learn" = "scikit-learn",
      "scipy" = "scipy>=1.10",
      "scvelo" = "scvelo",
      "wot" = "wot",
      "trimap" = "trimap",
      "pacmap" = "pacmap",
      "phate" = "phate",
      "bbknn" = "bbknn",
      "scanorama" = "scanorama",
      "scvi-tools" = "scvi-tools",
      "tbb" = "tbb"
    )
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
      "Cannot find the conda environment: ", envname,
      message_type = "error"
    )
  }

  log_message("Retrieving package list for environment: ", envname)

  tryCatch(
    {
      all_installed <- reticulate:::conda_list_packages(
        conda = conda,
        envname = envname,
        no_pip = FALSE
      )
      log_message("Found ", nrow(all_installed), " packages installed")
      return(all_installed)
    },
    error = function(e) {
      log_message("Failed to retrieve package list: ", e$message, message_type = "error")
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
      "Cannot find the conda environment: ", envname,
      message_type = "error"
    )
  }

  log_message("Checking ", length(packages), " packages in environment: ", envname)

  all_installed <- tryCatch(
    {
      installed_python_pkgs(envname = envname, conda = conda)
    },
    error = function(e) {
      log_message(
        "Failed to get installed packages: ", e$message,
        message_type = "warning"
      )
    }
  )

  packages_installed <- setNames(rep(FALSE, length(packages)), packages)
  packages_checked <- 0

  for (i in seq_along(packages)) {
    pkg <- packages[i]
    packages_checked <- packages_checked + 1

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
          log_message(pkg_name, " ", pkg_version, message_type = "success")
        } else {
          log_message(
            "! ", pkg_name, " found but version mismatch: installed=",
            installed_version, ", required=", pkg_version
          )
        }
      } else {
        packages_installed[pkg] <- TRUE
        installed_version <- all_installed$version[all_installed$package == pkg_name]
        log_message(pkg_name, " ", installed_version, message_type = "success")
      }
    } else {
      packages_installed[pkg] <- FALSE
      log_message(pkg_name, " not found", message_type = "warning")
    }

    if (length(packages) > 10 && packages_checked %% 5 == 0) {
      log_message(
        "Progress: ", packages_checked, "/", length(packages), " packages checked"
      )
    }
  }

  n_installed <- sum(packages_installed)
  n_total <- length(packages_installed)
  log_message(
    "Package check complete: ", n_installed, "/", n_total, " packages available"
  )

  return(packages_installed)
}

#' Check if a conda environment exists
#'
#' @param envs_dir Directories in which conda environments are located.
#' @inheritParams check_python
env_exist <- function(conda = "auto", envname = NULL, envs_dir = NULL) {
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
          log_message("Failed to get conda info: ", e$message, message_type = "warning")
          return(NULL)
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
        log_message("Environment directory exists but appears invalid: ", env_path,
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
    envname <- getOption("scop_env_name", default = "scop_env")
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
      log_message("packages argument is required", message_type = "error")
    }
  }

  conda <- reticulate::conda_binary(conda)
  envname <- reticulate:::condaenv_resolve(envname)

  log_message(
    "Installing ", length(packages), " packages into environment: ", envname
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
    log_message("Environment doesn't exist, creating: ", envname)
    reticulate::conda_create(
      envname = envname,
      packages = python_package %||% "python",
      forge = forge,
      channel = channel,
      conda = conda
    )
    python <- conda_python(envname = envname, conda = conda)
    log_message("Environment created with Python: ", python)
  }

  if (!is.null(python_version)) {
    log_message("Updating Python to version: ", python_version)
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
    log_message("Installing packages via pip...")
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
        log_message("Pip installation failed: ", e$message, message_type = "error")
      }
    )

    if (!is.null(result)) {
      log_message("Pip installation completed", message_type = "success")
    }
    return(result)
  }

  log_message("Installing packages via conda...")
  args <- reticulate:::conda_args("install", envname)

  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }

  for (ch in channels) {
    args <- c(args, "-c", ch)
    log_message("Using channel: ", ch)
  }

  args <- c(args, python_package, packages)

  if (length(packages) <= 10) {
    log_message("Installing packages: ", paste(packages, collapse = ", "))
  } else {
    log_message(
      "Installing ", length(packages), " packages (",
      paste(packages[1:3], collapse = ", "), "...)"
    )
  }

  result <- reticulate:::system2t(conda, shQuote(args))

  if (result != 0L) {
    log_message(
      "Conda installation failed with error code: ", result,
      message_type = "warning"
    )
  } else {
    log_message(
      "Conda installation completed successfully",
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
  envname <- reticulate:::condaenv_resolve(envname)
  if (grepl("[/\\\\]", envname)) {
    suffix <- if (reticulate:::is_windows()) "python.exe" else "bin/python"
    path <- file.path(envname, suffix)
    if (file.exists(path)) {
      return(path)
    }
    fmt <- "no conda environment exists at path '%s'"
    log_message(
      sprintf(fmt, envname),
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
      "conda environment \"", envname, "\" not found",
      message_type = "error"
    )
  }
  python <- if (all) env$python else env$python[[1L]]
  return(normalizePath(as.character(python), mustWork = FALSE))
}

run_python <- function(
    command,
    envir = .GlobalEnv) {
  tryCatch(
    expr = {
      eval(
        {
          reticulate::py_run_string(command)
        },
        envir = envir
      )
    },
    error = function(error) {
      log_message(error, message_type = "error")
      log_message(
        "Failed to run \"", command, "\". Please check manually.",
        message_type = "error"
      )
    }
  )
}

#' Try to evaluate an expression a set number of times before failing
#'
#' The function is used as a fail-safe if your R code sometimes works and sometimes
#' doesn't, usually because it depends on a resource that may be temporarily
#' unavailable. It tries to evaluate the expression `max_tries` times. If all the
#' attempts fail, it throws an error; if not, the evaluated expression is returned.
#'
#' @param expr The expression to be evaluated.
#' @param max_tries The maximum number of attempts to evaluate the expression before giving up. Default is set to 5.
#' @param error_message a string, additional custom error message you would like to be displayed when an error occurs.
#' @param retry_message a string, a message displayed when a new try to evaluate the expression would be attempted.
#'
#' @return This function returns the evaluated expression if successful,
#' otherwise it throws an error if all attempts are unsuccessful.
#' @export
#'
#' @examples
#' f <- function() {
#'   value <- runif(1, min = 0, max = 1)
#'   if (value > 0.5) {
#'     log_message("value is larger than 0.5")
#'     return(value)
#'   } else {
#'     log_message(
#'       "value is smaller than 0.5",
#'       message_type = "error"
#'     )
#'   }
#' }
#' f_evaluated <- try_get(expr = f())
#' print(f_evaluated)
try_get <- function(
    expr,
    max_tries = 5,
    error_message = "",
    retry_message = "Retrying...") {
  out <- simpleError("start")
  ntry <- 0
  while (inherits(out, "error")) {
    ntry <- ntry + 1
    out <- tryCatch(
      expr = eval.parent(substitute(expr)),
      error = function(error) {
        log_message(error)
        log_message("")
        log_message(error_message)
        Sys.sleep(1)
        return(error)
      }
    )
    if (inherits(out, "error") && ntry >= max_tries) {
      log_message(
        out,
        message_type = "error"
      )
    } else {
      if (!inherits(out, "error")) {
        break
      } else {
        log_message(retry_message)
      }
    }
  }
  return(out)
}

#' Download File from the Internet
#'
#' @md
#' @inheritParams utils::download.file
#' @param methods Methods to be used for downloading files.
#' The default is to try different download methods in turn until the download is successfully completed.
#' @param max_tries Number of tries for each download method.
#' @param ... Other arguments passed to [utils::download.file]
#'
#' @export
download <- function(
    url,
    destfile,
    methods = c("auto", "wget", "libcurl", "curl", "wininet", "internal"),
    quiet = FALSE,
    ...,
    max_tries = 2) {
  if (missing(url) || missing(destfile)) {
    log_message(
      "'url' and 'destfile' must be both provided.",
      message_type = "error"
    )
  }
  ntry <- 0
  status <- NULL
  while (is.null(status)) {
    for (method in methods) {
      status <- tryCatch(
        expr = {
          suppressWarnings(
            utils::download.file(
              url = url,
              destfile = destfile,
              method = method,
              quiet = quiet,
              ...
            )
          )
          status <- 1
        },
        error = function(error) {
          log_message(
            error,
            message_type = "warning"
          )
          log_message(
            "Cannot download from the url: ", url,
            message_type = "warning"
          )
          log_message(
            "Failed to download using \"", method, "\". Retry...",
            message_type = "warning"
          )
          Sys.sleep(1)
          return(NULL)
        }
      )
      if (!is.null(status)) {
        break
      }
    }
    ntry <- ntry + 1
    if (is.null(status) && ntry >= max_tries) {
      log_message(
        "Download failed.",
        message_type = "error"
      )
    }
  }
  return(invisible(NULL))
}

kegg_get <- function(url) {
  temp <- tempfile()
  on.exit(unlink(temp))
  download(url = url, destfile = temp)
  content <- as.data.frame(
    do.call(
      rbind,
      strsplit(readLines(temp), split = "\t")
    )
  )
  return(content)
}

rescale <- function(
    x,
    from = range(x, na.rm = TRUE, finite = TRUE),
    to = c(0, 1)) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  } else {
    return((x - from[1]) / diff(from) * diff(to) + to[1])
  }
}

zero_range <- function(
    x,
    tol = 1000 * .Machine$double.eps) {
  if (length(x) == 1) {
    return(TRUE)
  }
  if (length(x) != 2) {
    log_message(
      "x must be length 1 or 2",
      message_type = "error"
    )
  }
  if (any(is.na(x))) {
    return(NA)
  }
  if (x[1] == x[2]) {
    return(TRUE)
  }
  if (all(is.infinite(x))) {
    return(FALSE)
  }
  m <- min(abs(x))
  if (m == 0) {
    return(FALSE)
  }
  abs((x[1] - x[2]) / m) < tol
}

col2hex <- function(cname) {
  colMat <- grDevices::col2rgb(cname)
  grDevices::rgb(
    red = colMat[1, ] / 255,
    green = colMat[2, ] / 255,
    blue = colMat[3, ] / 255
  )
}

#' Invoke a function with a list of arguments
#' @param .fn A function, or function name as a string.
#' @param .args A list of arguments.
#' @param ... Other arguments passed to the function.
#' @param .env Environment in which to evaluate the call. This will be most useful if .fn is a string, or the function has side-effects.
#'
#' @export
invoke_fun <- function(
    .fn,
    .args = list(),
    ...,
    .env = rlang::caller_env()) {
  args <- c(.args, list(...))
  .bury <- c(".fn", "")
  if (rlang::is_null(.bury) || !length(args)) {
    if (rlang::is_scalar_character(.fn)) {
      .fn <- rlang::env_get(.env, .fn, inherit = TRUE)
    }
    call <- rlang::call2(.fn, !!!args)
    return(.External2(rlang:::ffi_eval, call, .env))
  }
  if (!rlang::is_character(.bury, 2L)) {
    rlang::abort("`.bury` must be a character vector of length 2")
  }
  arg_prefix <- .bury[[2]]
  fn_nm <- .bury[[1]]
  buried_nms <- paste0(arg_prefix, seq_along(args))
  buried_args <- rlang::set_names(args, buried_nms)
  .env <- rlang::env(.env, !!!buried_args)
  args <- rlang::set_names(buried_nms, names(args))
  args <- rlang::syms(args)
  if (rlang::is_function(.fn)) {
    rlang::env_bind(.env, `:=`(!!fn_nm, .fn))
    .fn <- fn_nm
  }
  call <- rlang::call2(.fn, !!!args)
  .External2(rlang:::ffi_eval, call, .env)
}

#' Implement similar functions to the \code{unnest} function in the tidyr package
#' @param data A data frame.
#' @param cols Columns to unnest.
#' @param keep_empty By default, you get one row of output for each element of the list your unchopping/unnesting.
#' This means that if there's a size-0 element (like \code{NULL} or an empty data frame),
#' that entire row will be dropped from the output.
#' If you want to preserve all rows,
#' use \code{keep_empty = TRUE} to replace size-0 elements with a single row of missing values.
#' @export
unnest <- function(
    data,
    cols,
    keep_empty = FALSE) {
  if (nrow(data) == 0 || length(cols) == 0) {
    return(data)
  }
  for (col in cols) {
    col_expand <- unlist(data[[col]])
    expand_times <- sapply(data[[col]], length)
    if (isTRUE(keep_empty)) {
      data[[col]][expand_times == 0] <- NA
      col_expand <- unlist(data[[col]])
      expand_times[expand_times == 0] <- 1
    }
    data <- data[rep(seq_len(nrow(data)), times = expand_times), ]
    data[, col] <- col_expand
  }
  rownames(data) <- NULL
  return(data)
}

#' Capitalizes the characters
#' Making the first letter uppercase
#'
#' @param x A vector of character strings to be capitalized.
#' @param force_tolower Whether to force the remaining letters to be lowercase.
#' @export
#'
#' @examples
#' x <- c(
#'   "dna methylation",
#'   "rRNA processing",
#'   "post-Transcriptional gene silencing"
#' )
#' capitalize(x)
capitalize <- function(x, force_tolower = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  if (!inherits(x, "character")) {
    log_message(
      "x must be the type of character.",
      message_type = "error"
    )
  }
  if (isTRUE(force_tolower)) {
    x <- paste(
      toupper(substr(x, 1, 1)),
      tolower(substr(x, 2, nchar(x))),
      sep = ""
    )
  } else {
    first_word <- sapply(strsplit(x, "\\s|-"), function(s) s[1])
    index <- which(first_word == tolower(first_word))
    x[index] <- paste(
      toupper(substr(x[index], 1, 1)),
      substr(x[index], 2, nchar(x[index])),
      sep = ""
    )
  }
  return(x)
}

str_wrap <- function(x, width = 80) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "factor")) {
    x <- as.character(x)
  }
  x_wrap <- unlist(
    lapply(
      x,
      function(i) {
        paste0(strwrap(i, width = width), collapse = "\n")
      }
    )
  )
  return(x_wrap)
}

select_cells <- function(obj, celltypes, group.by) {
  metadata <- obj@meta.data
  cells_c <- c()
  for (celltype in celltypes) {
    cells_c <- c(
      cells_c,
      rownames(metadata[metadata[[group.by]] == celltype, ])
    )
  }
  return(cells_c)
}
