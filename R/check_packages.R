#' @title Check and install python packages
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams PrepareEnv
#' @param packages A character vector, indicating package names which should be installed or removed.
#' Use `"<package>==<version>"` to request the installation of a specific version of a package.
#' @param force Whether to force package installation.
#' Default is `FALSE`.
#' @param pip Whether to use pip for package installation.
#' Default is `TRUE`, packages are installed from the active conda channels.
#' @param pip_options An optional character vector of additional command line arguments to be passed to `pip`.
#' Only relevant when `pip = TRUE`.
#' @param ... Other arguments passed to [reticulate::conda_install]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' check_python(
#'   packages = c("numpy", "pandas")
#' )
#'
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
      "{.arg envname}: {.val {envname}} python environment does not exist",
      message_type = "warning",
      verbose = verbose
    )
    PrepareEnv(envname = envname)
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

#' Check if the python package exists in the environment
#'
#' @inheritParams PrepareEnv
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
    "Checking {.val {length(packages)}} package{?s} in environment: {.file {envname}}"
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
  system2t <- get_namespace_fun("reticulate", "system2t")
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

          status <- system2t(python, shQuote(args))

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
        conda_args <- get_namespace_fun("reticulate", "conda_args")
        args <- conda_args("remove", envname)
        args <- c(args, packages)

        status <- system2t(conda, shQuote(args))

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

            status <- system2t(python, shQuote(args))

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
#' Package source can be *CRAN*, *Bioconductor* or *Github*.
#' By default, the package name is extracted according to the `packages` parameter.
#' @param lib The location of the library directories where to install the packages.
#' @param dependencies Whether to install dependencies of the packages.
#' Default is `TRUE`.
#' @param force Whether to force the installation of packages.
#' Default is `FALSE`.
#'
#' @return Package installation status.
#'
#' @export
#' @examples
#' check_r(c("ggplot2", "dplyr"))
check_r <- function(
    packages,
    lib = .libPaths()[1],
    dependencies = TRUE,
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
    check_pkg <- check_pkg_status(
      pkg_name,
      version = version,
      lib = lib
    )

    force_update <- FALSE
    if (check_pkg && !is.null(version)) {
      current_version <- utils::packageVersion(pkg_name)
      force_update <- current_version < package_version(version)
    }
    force_update <- force_update || isTRUE(force)

    if (!check_pkg || force_update) {
      log_message(
        "Installing: {.pkg {pkg_name}}...",
        message_type = "running",
        verbose = verbose
      )
      status_list[[pkg]] <- FALSE
      tryCatch(
        expr = {
          old_lib_paths <- .libPaths()
          .libPaths(lib)
          if (isTRUE(verbose)) {
            pak::pak(
              pkg,
              lib = lib,
              dependencies = dependencies
            )
          } else {
            invisible(
              suppressMessages(
                pak::pak(
                  pkg,
                  lib = lib,
                  dependencies = dependencies
                )
              )
            )
          }
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
      status_list[[pkg]] <- check_pkg_status(
        pkg_name,
        version = version,
        lib = lib
      )
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
      message_type = "success",
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
    pkg_installed <- check_pkg_status(pkg, lib = lib)

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
      status_list[[pkg]] <- !check_pkg_status(pkg, lib = lib)
    } else {
      log_message(
        "{.pkg {pkg}} is not installed, skipping removal",
        message_type = "warning",
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
      message_type = "success",
      verbose = verbose
    )
  }

  return(invisible(status_list))
}

check_pkg_status <- function(pkg, version = NULL, lib = .libPaths()[1]) {
  installed_pkgs_info <- utils::installed.packages(lib.loc = lib)
  installed_pkgs <- installed_pkgs_info[, "Package"]
  installed_pkgs_version <- installed_pkgs_info[, "Version"]
  pkg_exists <- pkg %in% installed_pkgs
  if (is.null(version)) {
    version_match <- TRUE
  } else {
    version_match <- installed_pkgs_version[installed_pkgs == pkg] == version
  }

  if (isFALSE(pkg_exists) || isFALSE(version_match)) {
    return(FALSE)
  }

  TRUE
}
