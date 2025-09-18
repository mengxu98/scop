#' @title Check and install python packages
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param packages A character vector, indicating package names which should be installed or removed.
#' Use `"<package>==<version>"` to request the installation of a specific version of a package.
#' @param envname The name of a conda environment.
#' @param conda The path to a conda executable. Use `"auto"` to allow scop to automatically find an appropriate conda binary.
#' @param force Whether to force package installation. Default is FALSE.
#' @param pip Whether to use pip for package installation.
#' Default is TRUE, packages are installed from the active conda channels.
#' @param pip_options An optional character vector of additional command line arguments to be passed to `pip`.
#' Only relevant when `pip = TRUE`.
#' @param ... Other arguments passed to [reticulate::conda_install]
#'
#' @export
#'
#' @examples
#' check_python(
#'   packages = c("numpy", "pandas")
#' )
#'
#' \dontrun{
#' check_python(
#'   packages = "numpy==1.26.4",
#'   envname = "scop_env",
#'   pip_options = "-i https://pypi.tuna.tsinghua.edu.cn/simple"
#' )
#' }
check_python <- function(
    packages,
    envname = NULL,
    conda = "auto",
    force = FALSE,
    pip = TRUE,
    pip_options = character(),
    verbose = TRUE,
    ...) {
  envname <- get_envname(envname)
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    options(reticulate.conda_binary = conda)
    conda <- find_conda()
  }
  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    log_message(
      "{.arg envname}: {.val {envname}} python environment does not exist. Create it with {.fn PrepareEnv}",
      message_type = "warning",
      verbose = verbose
    )
    PrepareEnv()
  }

  if (isTRUE(force)) {
    pkg_installed <- stats::setNames(
      rep(FALSE, length(packages)),
      packages
    )
    pip_options <- c(pip_options, "--force-reinstall")
  } else {
    pkg_installed <- exist_python_pkgs(
      packages = packages,
      envname = envname,
      conda = conda
    )
  }

  if (sum(!pkg_installed) > 0) {
    pkgs_to_install <- names(pkg_installed)[!pkg_installed]
    log_message(
      "Try to install: {.pkg {pkgs_to_install}}",
      verbose = verbose
    )
    if (isTRUE(pip)) {
      pkgs_to_install <- c("pip", pkgs_to_install)
    }
    tryCatch(
      expr = {
        conda_install(
          conda = conda,
          packages = pkgs_to_install,
          envname = envname,
          pip = pip,
          pip_options = pip_options,
          ...
        )
      },
      error = identity
    )

    pkg_installed <- exist_python_pkgs(
      packages = packages,
      envname = envname,
      conda = conda
    )
  }

  if (sum(!pkg_installed) > 0) {
    failed_pkgs <- names(pkg_installed)[!pkg_installed]
    log_message(
      "Failed to install: {.pkg {failed_pkgs}} into the environment {.file {envname}}. Please install manually",
      message_type = "warning",
      verbose = verbose
    )
  } else {
    return(invisible(NULL))
  }
}

#' @title Remove Python packages from conda environment
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param packages A character vector of package names to remove.
#' @param envname The name of the conda environment.
#' If `NULL`, uses the default scop environment name.
#' @param conda The path to a conda executable.
#' Use `"auto"` to allow reticulate to automatically find an appropriate conda binary.
#' @param pip Whether to use pip for package removal.
#' Default is `FALSE` (use conda).
#' @param force Whether to force removal without confirmation.
#' Default is `FALSE`.
#'
#' @return Invisibly value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove a single package using conda
#' remove_python("numpy")
#'
#' # Remove multiple packages using conda
#' remove_python(c("numpy", "pandas"))
#'
#' # Remove packages using pip
#' remove_python("numpy", pip = TRUE)
#'
#' # Force removal without confirmation
#' remove_python("numpy", force = TRUE)
#'
#' # Remove packages from a specific environment
#' remove_python("numpy", envname = "env")
#' }
remove_python <- function(
    packages,
    envname = NULL,
    conda = "auto",
    pip = FALSE,
    force = FALSE,
    verbose = TRUE) {
  envname <- get_envname(envname)

  log_message(
    "Removing {.pkg {packages}} from environment: {.file {envname}}",
    verbose = verbose
  )

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
      "Cannot find the conda environment: {.file {envname}}",
      message_type = "error"
    )
    return(invisible(FALSE))
  }

  if (!force) {
    if (interactive()) {
      response <- readline(
        paste0(
          "Are you sure you want to remove these packages from environment '",
          envname, "'? (y/N): "
        )
      )
      if (!tolower(response) %in% c("y", "yes")) {
        log_message(
          "{.pkg {packages}} removal cancelled",
          message_type = "warning",
          verbose = verbose
        )
        return(invisible(FALSE))
      }
    } else {
      log_message(
        "Automatically remove {.pkg {packages}} in non-interactive mode",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  python <- tryCatch(
    {
      conda_python(envname = envname, conda = conda)
    },
    error = function(e) {
      log_message(
        "Failed to get Python path: {.val {e$message}}",
        message_type = "error"
      )
      NULL
    }
  )

  if (is.null(python)) {
    return(invisible(FALSE))
  }

  if (pip) {
    log_message(
      "Removing {.pkg {packages}} via {.pkg pip}...",
      verbose = verbose
    )

    result <- tryCatch(
      {
        for (pkg in packages) {
          log_message(
            "Removing {.pkg {pkg}}",
            verbose = verbose
          )

          args <- c(
            "-m", "pip", "uninstall", "-y", pkg
          )

          status <- reticulate:::system2t(python, shQuote(args))

          if (status != 0L) {
            log_message(
              "Failed to remove {.pkg {pkg}} via {.pkg pip} [error code {.val {status}}]",
              message_type = "warning",
              verbose = verbose
            )
          } else {
            log_message(
              "Package {.pkg {pkg}} removed successfully via {.pkg pip}",
              message_type = "success",
              verbose = verbose
            )
          }
        }
        TRUE
      },
      error = function(e) {
        log_message(
          "Pip removal failed: {.val {e$message}}",
          message_type = "error"
        )
        FALSE
      }
    )
  } else {
    log_message(
      "Removing {.pkg {packages}} via {.pkg conda}...",
      verbose = verbose
    )

    result <- tryCatch(
      {
        args <- reticulate:::conda_args("remove", envname)
        args <- c(args, packages)

        status <- reticulate:::system2t(conda, shQuote(args))

        if (status != 0L) {
          log_message(
            "{.pkg {packages}} removal failed via {.pkg conda} with error code: {.val {status}}",
            message_type = "warning",
            verbose = verbose
          )
          FALSE
        } else {
          log_message(
            "{.pkg {packages}} removed successfully via {.pkg conda}",
            message_type = "success",
            verbose = verbose
          )
          TRUE
        }
      },
      error = function(e) {
        log_message(
          "Conda removal failed: {.val {e$message}}",
          message_type = "warning",
          verbose = verbose
        )
        FALSE
      }
    )

    if (!result && !pip) {
      log_message(
        "{.pkg {packages}} removal failed via {.pkg conda}, trying {.pkg pip} as fallback...",
        message_type = "warning",
        verbose = verbose
      )

      result <- tryCatch(
        {
          for (pkg in packages) {
            log_message(
              "Removing {.pkg {pkg}}",
              verbose = verbose
            )

            args <- c(
              "-m", "pip", "uninstall", "-y", pkg
            )

            status <- reticulate:::system2t(python, shQuote(args))

            if (status != 0L) {
              log_message(
                "Failed to remove {.pkg {pkg}} via {.pkg pip} [error code {.val {status}}]",
                message_type = "warning",
                verbose = verbose
              )
            } else {
              log_message(
                "{.pkg {pkg}} removed successfully via {.pkg pip}",
                message_type = "success",
                verbose = verbose
              )
            }
          }
          TRUE
        },
        error = function(e) {
          log_message(
            "{.pkg {packages}} removal failed via {.pkg pip} as fallback: {.val {e$message}}",
            message_type = "error",
            verbose = verbose
          )
          FALSE
        }
      )
    }
  }

  if (result) {
    log_message(
      "{.pkg {packages}} removal completed successfully",
      message_type = "success",
      verbose = verbose
    )
  } else {
    log_message(
      "{.pkg {packages}} removal failed",
      message_type = "warning",
      verbose = verbose
    )
  }

  return(invisible(result))
}

#' @title Check and install R packages
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param packages Package to be installed.
#' Package source can be CRAN, Bioconductor or Github.
#' By default, the package name is extracted according to the `packages` parameter.
#' @param lib The location of the library directories where to install the packages.
#' @param force Whether to force the installation of packages.
#' Default is `FALSE`.
#'
#' @export
#' @examples
#' check_r(c("Seurat", "reticulate"))
check_r <- function(
    packages,
    lib = .libPaths()[1],
    force = FALSE,
    verbose = TRUE) {
  status_list <- list()
  for (pkg in packages) {
    version <- NULL
    if (grepl("/", pkg)) {
      pkg_name <- strsplit(pkg, split = "/|@|==", perl = TRUE)[[1]][[2]]
    } else {
      pkg_info <- strsplit(pkg, split = "@|==", perl = TRUE)[[1]]
      pkg_name <- pkg_info[[1]]
      if (length(pkg_info) > 1) {
        version <- pkg_info[[2]]
      }
    }
    dest <- gsub("@.*|==.*|>=.*", "", pkg)

    check_pkg <- .check_pkg_status(pkg_name, lib = lib)

    force_update <- FALSE
    if (check_pkg && !is.null(version)) {
      current_version <- utils::packageVersion(pkg_name)
      force_update <- current_version < package_version(version)
    }
    force_update <- force_update || isTRUE(force)

    if (!check_pkg || force_update) {
      log_message(
        "Installing: {.pkg {pkg_name}}...",
        verbose = verbose
      )
      status_list[[pkg]] <- FALSE
      tryCatch(
        expr = {
          old_lib_paths <- .libPaths()
          .libPaths(lib)
          pak::pak(dest)
          .libPaths(old_lib_paths)
        },
        error = function(e) {
          status_list[[pkg]] <- FALSE
          log_message(
            "Failed to install: {.pkg {pkg_name}}. Error: {.val {e$message}}",
            message_type = "warning",
            verbose = verbose
          )
        }
      )
      status_list[[pkg]] <- .check_pkg_status(pkg_name, lib = lib)
    } else {
      status_list[[pkg]] <- TRUE
    }
  }

  success <- sapply(status_list, isTRUE)
  failed <- names(status_list)[!success]

  if (length(failed) > 0) {
    log_message(
      "Failed to install: {.pkg {failed}}. Please install manually",
      message_type = "warning",
      verbose = verbose
    )
  } else {
    log_message(
      "{.pkg {packages}} installed successfully",
      verbose = verbose
    )
  }

  return(invisible(status_list))
}

#' @title Check and remove R packages
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param packages Package to be removed.
#' @param lib The location of the library directories where to remove the packages.
#'
#' @export
remove_r <- function(
    packages,
    lib = .libPaths()[1],
    verbose = TRUE) {
  status_list <- list()
  for (pkg in packages) {
    pkg_installed <- .check_pkg_status(pkg, lib = lib)

    if (pkg_installed) {
      log_message(
        "Removing: {.pkg {pkg}}...",
        verbose = verbose
      )
      status_list[[pkg]] <- FALSE
      tryCatch(
        expr = {
          old_lib_paths <- .libPaths()
          .libPaths(lib)
          pak::pkg_remove(pkg)
          .libPaths(old_lib_paths)
        },
        error = function(e) {
          log_message(
            "Warning during removal: {.pkg {pkg}}. Error: {.val {e$message}}",
            message_type = "warning",
            verbose = verbose
          )
        }
      )
      status_list[[pkg]] <- !.check_pkg_status(pkg, lib = lib)
    } else {
      log_message(
        "{.pkg {pkg}} is not installed, skipping removal",
        verbose = verbose
      )
      status_list[[pkg]] <- TRUE
    }
  }

  success <- sapply(status_list, isTRUE)
  failed <- names(status_list)[!success]

  if (length(failed) > 0) {
    log_message(
      "Failed to remove: {.pkg {failed}}. Please remove manually",
      message_type = "warning",
      verbose = verbose
    )
  } else {
    log_message(
      "{.pkg {packages}} removed successfully",
      verbose = verbose
    )
  }

  return(invisible(status_list))
}

.check_pkg_status <- function(pkg, lib = .libPaths()[1]) {
  installed_pkgs <- utils::installed.packages(lib.loc = lib)
  installed_pkgs <- installed_pkgs[, "Package"]
  pkg_exists <- pkg %in% installed_pkgs

  if (isFALSE(pkg_exists)) {
    return(FALSE)
  }

  TRUE
}
