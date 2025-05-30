#' This function prepares the scop Python environment by installing the required dependencies and setting up the environment.
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
    version = "3.8-1",
    force = FALSE,
    ...) {
  envname <- get_envname(envname)

  requirements <- env_requirements(version = version)
  python_version <- requirements[["python"]]
  packages <- requirements[["packages"]]

  if (!is.null(conda)) {
    if (identical(conda, "auto")) {
      conda <- find_conda()
    } else {
      options(reticulate.conda_binary = conda)
      conda <- find_conda()
    }
  }

  if (is.null(conda)) {
    env <- FALSE
  } else {
    envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    env <- env_exist(
      conda = conda,
      envname = envname,
      envs_dir = envs_dir
    )
    if (isTRUE(force) && isTRUE(env)) {
      unlink(paste0(envs_dir, "/", envname), recursive = TRUE)
      env <- FALSE
    }
  }

  if (isTRUE(env)) {
    python_path <- conda_python(conda = conda, envname = envname)
  } else {
    force <- TRUE
    if (is.null(conda)) {
      message("Conda not found. Installing miniconda...")
      options(timeout = 360)
      version <- "3"
      info <- as.list(Sys.info())
      if (info$sysname == "Darwin" && info$machine == "arm64") {
        base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
        name <- "Miniforge3-MacOSX-arm64.sh"
        return(file.path(base, name))
      }
      base <- miniconda_repo
      info <- as.list(Sys.info())
      arch <- reticulate:::miniconda_installer_arch(info)
      version <- as.character(version)
      name <- if (reticulate:::is_windows()) {
        sprintf(
          "Miniconda%s-latest-Windows-%s.exe",
          version,
          arch
        )
      } else if (reticulate:::is_osx()) {
        sprintf(
          "Miniconda%s-latest-MacOSX-%s.sh",
          version,
          arch
        )
      } else if (reticulate:::is_linux()) {
        sprintf("Miniconda%s-latest-Linux-%s.sh", version, arch)
      } else {
        BBmisc::stopf(
          "unsupported platform %s",
          shQuote(Sys.info()[["sysname"]])
        )
      }
      url <- file.path(base, name)
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
      unlink(miniconda_path, recursive = TRUE)
      reticulate::install_miniconda(
        path = miniconda_path,
        force = TRUE,
        update = FALSE
      )
      conda <- reticulate:::miniconda_conda(miniconda_path)
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    }
    if (
      python_version < numeric_version("3.7.0") ||
        python_version >= numeric_version("3.10.0")
    ) {
      stop("scop currently only support python version 3.7-3.9!")
    }
    python_path <- reticulate::conda_create(
      conda = conda,
      envname = envname,
      python_version = python_version,
      packages = "pytables"
    )
    env_path <- paste0(envs_dir, "/", envname)
    env <- file.exists(env_path)
    if (isFALSE(env)) {
      print(reticulate:::conda_info(conda = conda))
      print(reticulate::conda_list(conda = conda))
      stop(
        "Unable to find scop environment under the expected path: ",
        env_path,
        "\n",
        "conda: ",
        conda,
        "\n",
        "scop python: ",
        python_path,
        "\n"
      )
    }
  }

  check_python(
    packages = packages,
    envname = envname,
    conda = conda,
    force = force,
    ...
  )

  Sys.unsetenv("RETICULATE_PYTHON")
  python_path <- conda_python(conda = conda, envname = envname)
  reticulate::use_python(python_path, required = TRUE)

  pyinfo <- utils::capture.output(reticulate::py_config())
  pyinfo_mesg <- c(
    "====================== scop conda environment ======================",
    paste0("conda: ", conda),
    paste0("environment: ", paste0(envs_dir, "/", get_envname())),
    "======================== scop python config ========================",
    pyinfo,
    "==================================================================="
  )
  invisible(
    lapply(pyinfo_mesg, packageStartupMessage)
  )
  invisible(
    run_python(
      command = "import matplotlib",
      envir = .GlobalEnv
    )
  )
  if (!interactive()) {
    invisible(
      run_python(
        command = "matplotlib.use('pdf')",
        envir = .GlobalEnv
      )
    )
  }
  invisible(
    run_python(
      command = "import matplotlib.pyplot as plt",
      envir = .GlobalEnv
    )
  )
  invisible(
    run_python(
      command = "import scanpy", envir = .GlobalEnv
    )
  )
}

#' env_requirements function
#'
#' This function provides the scop python environment requirements for a specific version.
#'
#' @param version A character vector specifying the version of the environment (default is "3.8-1").
#' @return A list of requirements for the specified version.
#' @details The function returns a list of requirements including the required Python version
#'          and a list of packages with their corresponding versions.
#' @examples
#' # Get requirements for version "3.8-1"
#' env_requirements("3.8-1")
#'
#' @export
env_requirements <- function(version = "3.8-1") {
  version <- match.arg(
    version,
    choices = c("3.8-1", "3.8-2", "3.9-1", "3.10-1", "3.11-1")
  )
  requirements <- switch(version,
    "3.8-1" = list(
      python = "3.8",
      packages = c(
        "leidenalg" = "leidenalg==0.10.1",
        "matplotlib" = "matplotlib==3.6.3",
        "numba" = "numba==0.55.2",
        "numpy" = "numpy==1.21.6",
        "palantir" = "palantir==1.0.1",
        "pandas" = "pandas==1.3.5",
        "python-igraph" = "python-igraph==0.10.2",
        "scanpy" = "scanpy==1.9.5",
        "scikit-learn" = "scikit-learn==1.3.2",
        "scipy" = "scipy==1.10.1",
        "scvelo" = "scvelo==0.2.5",
        "wot" = "wot==1.0.8.post2",
        "trimap" = "trimap==1.1.4",
        "pacmap" = "pacmap==0.7.0",
        "phate" = "phate==1.0.11",
        "bbknn" = "bbknn==1.6.0",
        "scanorama" = "scanorama==1.7.4",
        "scvi-tools" = "scvi-tools==0.20.3"
      )
    ),
    "3.8-2" = list(
      python = "3.8",
      packages = c(
        "leidenalg" = "leidenalg==0.10.1",
        "matplotlib" = "matplotlib==3.7.3",
        "numba" = "numba==0.58.1",
        "numpy" = "numpy==1.24.4",
        "palantir" = "palantir==1.3.0",
        "pandas" = "pandas==1.5.3",
        "python-igraph" = "python-igraph==0.10.8",
        "scanpy" = "scanpy==1.9.5",
        "scikit-learn" = "scikit-learn==1.3.2",
        "scipy" = "scipy==1.10.1",
        "scvelo" = "scvelo==0.2.5",
        "wot" = "wot==1.0.8.post2",
        "trimap" = "trimap==1.1.4",
        "pacmap" = "pacmap==0.7.0",
        "phate" = "phate==1.0.11",
        "bbknn" = "bbknn==1.6.0",
        "scanorama" = "scanorama==1.7.4",
        "scvi-tools" = "scvi-tools==0.20.3"
      )
    ),
    "3.9-1" = list(
      python = "3.9",
      packages = c(
        "leidenalg" = "leidenalg==0.10.1",
        "matplotlib" = "matplotlib==3.8.0",
        "numba" = "numba==0.58.1",
        "numpy" = "numpy==1.25.2",
        "palantir" = "palantir==1.3.0",
        "pandas" = "pandas==1.5.3",
        "python-igraph" = "python-igraph==0.10.8",
        "scanpy" = "scanpy==1.9.5",
        "scikit-learn" = "scikit-learn==1.3.2",
        "scipy" = "scipy==1.11.3",
        "scvelo" = "scvelo==0.2.5",
        "wot" = "wot==1.0.8.post2",
        "trimap" = "trimap==1.1.4",
        "pacmap" = "pacmap==0.7.0",
        "phate" = "phate==1.0.11",
        "bbknn" = "bbknn==1.6.0",
        "scanorama" = "scanorama==1.7.4",
        "scvi-tools" = "scvi-tools==0.20.3"
      )
    ),
    "3.10-1" = list(
      python = "3.10",
      packages = c(
        "leidenalg" = "leidenalg==0.10.1",
        "matplotlib" = "matplotlib==3.8.0",
        "numba" = "numba==0.58.1",
        "numpy" = "numpy==1.25.2",
        "palantir" = "palantir==1.3.0",
        "pandas" = "pandas==1.5.3",
        "python-igraph" = "python-igraph==0.10.8",
        "scanpy" = "scanpy==1.9.5",
        "scikit-learn" = "scikit-learn==1.3.2",
        "scipy" = "scipy==1.11.3",
        "scvelo" = "scvelo==0.2.5",
        "wot" = "wot==1.0.8.post2",
        "trimap" = "trimap==1.1.4",
        "pacmap" = "pacmap==0.7.0",
        "phate" = "phate==1.0.11",
        "bbknn" = "bbknn==1.6.0",
        "scanorama" = "scanorama==1.7.4",
        "scvi-tools" = "scvi-tools==0.20.3"
      )
    ),
    "3.11-1" = list(
      python = "3.10",
      packages = c(
        "leidenalg" = "leidenalg==0.10.1",
        "matplotlib" = "matplotlib==3.8.0",
        "numba" = "numba==0.58.1",
        "numpy" = "numpy==1.25.2",
        "palantir" = "palantir==1.3.0",
        "pandas" = "pandas==1.5.3",
        "python-igraph" = "python-igraph==0.10.8",
        "scanpy" = "scanpy==1.9.5",
        "scikit-learn" = "scikit-learn==1.3.2",
        "scipy" = "scipy==1.11.3",
        "scvelo" = "scvelo==0.2.5",
        "wot" = "wot==1.0.8.post2",
        "trimap" = "trimap==1.1.4",
        "pacmap" = "pacmap==0.7.0",
        "phate" = "phate==1.0.11",
        "bbknn" = "bbknn==1.6.0",
        "scanorama" = "scanorama==1.7.4",
        "scvi-tools" = "scvi-tools==0.20.3"
      )
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
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    stop("Can not find the conda environment: ", envname)
  }
  all_installed <- reticulate:::conda_list_packages(
    conda = conda,
    envname = envname,
    no_pip = FALSE
  )
  return(all_installed)
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
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    stop("Can not find the conda environment: ", envname)
  }
  all_installed <- installed_python_pkgs(envname = envname, conda = conda)
  packages_installed <- NULL
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
        packages_installed[pkg] <- all_installed$version[
          all_installed$package == pkg_name
        ] ==
          pkg_version
      } else {
        packages_installed[pkg] <- TRUE
      }
    } else {
      packages_installed[pkg] <- FALSE
    }
  }
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
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
    }
    exist <- file.exists(paste0(envs_dir, "/", envname))
  } else {
    exist <- FALSE
  }
  return(exist)
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

#' Installs a list of packages into a specified conda environment
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
      BBmisc::stopf(fmt, deparse1(substitute(envname)), call. = FALSE)
    } else {
      packages
    }
  }
  conda <- reticulate::conda_binary(conda)
  envname <- reticulate:::condaenv_resolve(envname)
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
    reticulate::conda_create(
      envname = envname,
      packages = python_package %||% "python",
      forge = forge,
      channel = channel,
      conda = conda
    )
    python <- conda_python(envname = envname, conda = conda)
  }
  if (!is.null(python_version)) {
    args <- reticulate:::conda_args("install", envname, python_package)
    status <- reticulate:::system2t(conda, shQuote(args))
    if (status != 0L) {
      fmt <- "installation of '%s' into environment '%s' failed [error code %i]"
      msg <- sprintf(fmt, python_package, envname, status)
      stop(msg, call. = FALSE)
    }
  }
  if (pip) {
    result <- reticulate:::pip_install(
      python = python,
      packages = packages,
      pip_options = pip_options,
      ignore_installed = pip_ignore_installed,
      conda = conda,
      envname = envname
    )
    return(result)
  }
  args <- reticulate:::conda_args("install", envname)
  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }
  for (ch in channels) args <- c(args, "-c", ch)
  args <- c(args, python_package, packages)
  result <- reticulate:::system2t(conda, shQuote(args))
  if (result != 0L) {
    fmt <- "one or more Python packages failed to install [error code %i]"
    BBmisc::stopf(fmt, result)
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
    stop(sprintf(fmt, envname))
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
    stop("conda environment \"", envname, "\" not found")
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
      message(error)
      stop("Failed to run \"", command, "\". Please check manually.")
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
#'     message("value is larger than 0.5")
#'     return(value)
#'   } else {
#'     stop("value is smaller than 0.5")
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
        message(error)
        message("")
        message(error_message)
        Sys.sleep(1)
        return(error)
      }
    )
    if (inherits(out, "error") && ntry >= max_tries) {
      stop(out)
    } else {
      if (!inherits(out, "error")) {
        break
      } else {
        message(retry_message)
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
    stop("'url' and 'destfile' must be both provided.")
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
          message(error)
          message("Cannot download from the url: ", url)
          message("Failed to download using \"", method, "\". Retry...\n")
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
      stop("Download failed.")
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
    stop("x must be length 1 or 2")
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
    stop("x must be the type of character.")
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
