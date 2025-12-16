#' @title Check and install python packages
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams PrepareEnv
#' @param packages A character vector of package names to check and install.
#' Use `"<package>==<version>"` to request a specific version.
#' @param force Whether to force package reinstallation.
#' Default is `FALSE`.
#' @param pip Whether to use `pip`/`uv` (`TRUE`) or `conda` (`FALSE`) for installation.
#' Default is `TRUE`. When `TRUE`, uv is used as the primary installer with pip as fallback.
#' @param pip_options Additional command line arguments to be passed to `uv`/`pip` when `pip = TRUE`.
#' @param ... Other arguments to be passed to [conda_install()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PrepareEnv()
#'
#' # Then check/install packages
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
  conda <- resolve_conda(conda)

  env <- env_exist(conda = conda, envname = envname)
  if (isFALSE(env)) {
    log_message(
      "Python environment {.val {envname}} not found. Run {.fn PrepareEnv} first",
      message_type = "error"
    )
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

exist_python_pkgs <- function(
    packages,
    envname = NULL,
    conda = "auto") {
  envname <- get_envname(envname)
  conda <- resolve_conda(conda)

  if (!ensure_conda(conda)) {
    return(invisible(FALSE))
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

  requirements <- env_requirements()
  pkg_name_mapping <- requirements$package_aliases

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

    check_pkg_names <- pkg_name
    if (length(pkg_name_mapping) > 0) {
      if (pkg_name %in% names(pkg_name_mapping)) {
        actual_name <- pkg_name_mapping[[pkg_name]]
        check_pkg_names <- c(pkg_name, actual_name)
      }
      if (pkg_name %in% pkg_name_mapping) {
        logical_name <- names(pkg_name_mapping)[pkg_name_mapping == pkg_name]
        check_pkg_names <- unique(c(check_pkg_names, logical_name))
      }
    }

    found_pkg_name <- NULL
    if (length(check_pkg_names) > 1) {
      for (chk_name in check_pkg_names) {
        if (chk_name %in% all_installed$package) {
          found_pkg_name <- chk_name
          break
        }
      }
    } else {
      if (check_pkg_names %in% all_installed$package) {
        found_pkg_name <- check_pkg_names
      }
    }

    if (!is.null(found_pkg_name)) {
      if (!is.na(pkg_version)) {
        installed_version <- all_installed$version[all_installed$package == found_pkg_name]

        get_major_version <- function(version_str) {
          tryCatch(
            {
              ver <- package_version(version_str)
              as.numeric(ver[1, 1])
            },
            error = function(e) {
              match_result <- regmatches(version_str, regexpr("^([0-9]+)", version_str))
              if (length(match_result) > 0) {
                as.numeric(match_result)
              } else {
                NA
              }
            }
          )
        }

        version_match <- tryCatch(
          {
            installed_ver <- package_version(installed_version)
            required_ver <- package_version(pkg_version)
            installed_ver == required_ver
          },
          error = function(e) {
            normalize_version_string <- function(v) {
              v <- gsub("(\\.0+)+$", "", v)
              v
            }
            installed_norm <- normalize_version_string(installed_version)
            required_norm <- normalize_version_string(pkg_version)
            installed_norm == required_norm || installed_version == pkg_version
          }
        )

        if (!version_match) {
          installed_major <- get_major_version(installed_version)
          required_major <- get_major_version(pkg_version)

          if (!is.na(installed_major) && !is.na(required_major) &&
            installed_major == required_major) {
            version_match <- TRUE
            log_message(
              "{.pkg {pkg_name}} compatible (major version {.pkg {required_major}}): installed {.pkg {installed_version}}, required {.pkg {pkg_version}}",
              message_type = "warning"
            )
          } else {
            log_message(
              "{.pkg {pkg_name}} found but version mismatch: installed {.pkg {installed_version}}, required {.pkg {pkg_version}}",
              message_type = "warning"
            )
          }
        } else {
          log_message(
            "{.pkg {pkg_name}} {.pkg {pkg_version}}",
            message_type = "success"
          )
        }

        packages_installed[pkg] <- version_match
      } else {
        packages_installed[pkg] <- TRUE
        installed_version <- all_installed$version[all_installed$package == found_pkg_name]
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

installed_python_pkgs <- function(
    envname = NULL,
    conda = "auto") {
  envname <- get_envname(envname)
  conda <- resolve_conda(conda)

  if (!ensure_conda(conda)) {
    return(invisible(NULL))
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
      all_installed <- get_namespace_fun(
        "reticulate", "conda_list_packages"
      )(
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

#' @title Remove Python packages from conda environment
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams PrepareEnv
#' @param packages A character vector of package names to remove.
#' @param pip Whether to use pip for package removal.
#' Default is `FALSE` (use conda).
#' @param force Whether to force removal without confirmation.
#' Default is `FALSE`.
#'
#' @return Invisibly returns.
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

  conda <- resolve_conda(conda)
  if (!ensure_conda(conda)) {
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
    result <- remove_via_pip_uv(packages, python, envname, conda, verbose)
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
        "{.pkg {packages}} removal failed via {.pkg conda}, trying {.pkg uv} as fallback...",
        message_type = "warning",
        verbose = verbose
      )

      result <- tryCatch(
        {
          uv <- find_uv(
            python = python, envname = envname, conda = conda, auto_install = TRUE
          )

          if (is.null(uv)) {
            log_message(
              "{.pkg uv} not found and installation failed, falling back to {.pkg pip}",
              message_type = "warning"
            )
            for (pkg in packages) {
              log_message(
                "Removing {.pkg {pkg}} via {.pkg pip}...",
                verbose = verbose
              )

              args <- c(
                "-m", "pip", "uninstall", "-y", pkg
              )

              status <- system2t(python, shQuote(args))

              if (status != 0L) {
                log_message(
                  "{.pkg {pkg}} removal failed via {.pkg pip} [error code {.val {status}}]",
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
          } else {
            for (pkg in packages) {
              log_message(
                "Removing {.pkg {pkg}}",
                verbose = verbose
              )

              args <- c("pip", "uninstall", "--python", python, "-y", pkg)

              if (uv == "python -m uv") {
                args <- c("-m", "uv", args)
                status <- system2t(python, shQuote(args))
              } else {
                status <- system2t(uv, shQuote(args))
              }

              if (status != 0L) {
                log_message(
                  "Failed to remove {.pkg {pkg}} via {.pkg uv} [error code {.val {status}}], trying pip",
                  message_type = "warning",
                  verbose = verbose
                )
                args <- c("-m", "pip", "uninstall", "-y", pkg)
                status <- system2t(python, shQuote(args))
                if (status == 0L) {
                  log_message(
                    "{.pkg {pkg}} removed successfully via {.pkg pip}",
                    message_type = "success",
                    verbose = verbose
                  )
                }
              } else {
                log_message(
                  "{.pkg {pkg}} removed successfully via {.pkg uv}",
                  message_type = "success",
                  verbose = verbose
                )
              }
            }
          }
          TRUE
        },
        error = function(e) {
          log_message(
            "{.pkg {packages}} removal failed via {.pkg uv} as fallback: {.val {e$message}}",
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

remove_via_pip_uv <- function(
    packages, python, envname, conda, verbose = TRUE) {
  system2t <- get_namespace_fun("reticulate", "system2t")

  uv <- find_uv(python = python, envname = envname, conda = conda, auto_install = TRUE)

  cmd_args <- c("pip", "uninstall", "-y")

  if (!is.null(uv)) {
    log_message(
      "Removing {.pkg {packages}} via {.pkg uv}...",
      verbose = verbose
    )

    if (uv == "python -m uv") {
      runner <- python
      base_args <- c("-m", "uv", cmd_args)
    } else {
      runner <- uv
      base_args <- c(cmd_args, "--python", python)
    }

    final_args <- c(base_args, packages)
    status <- system2t(runner, shQuote(final_args))

    if (status == 0L) {
      log_message(
        "{.pkg {packages}} removed successfully via {.pkg uv}",
        message_type = "success",
        verbose = verbose
      )
      return(TRUE)
    } else {
      log_message(
        "{.pkg uv} removal failed [error code {.val {status}}], falling back to pip",
        message_type = "warning",
        verbose = verbose
      )
    }
  } else {
    log_message(
      "uv not found, falling back to pip",
      message_type = "warning", verbose = verbose
    )
  }

  log_message(
    "Removing {.pkg {packages}} via {.pkg pip}...",
    verbose = verbose
  )
  args <- c("-m", "pip", "uninstall", "-y", packages)
  status <- system2t(python, shQuote(args))

  if (status == 0L) {
    log_message(
      "{.pkg {packages}} removed successfully via {.pkg pip}",
      message_type = "success",
      verbose = verbose
    )
    return(TRUE)
  } else {
    log_message(
      "Failed to remove {.pkg {packages}} via {.pkg pip} [error code {.val {status}}]",
      message_type = "warning",
      verbose = verbose
    )
    return(FALSE)
  }
}
