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
#' @param ... Other arguments to be passed to other functions.
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
  ...
) {
  packages <- resolve_requested_python_packages(packages)
  envname <- get_envname(envname)
  conda <- resolve_conda(conda)
  pip_options <- normalize_cli_args(pip_options)

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
      conda = conda,
      verbose = verbose
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
      conda = conda,
      verbose = verbose
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
  conda = "auto",
  verbose = TRUE
) {
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
    "Checking {.val {length(packages)}} package{?s} in environment: {.file {envname}}",
    verbose = verbose
  )

  all_installed <- tryCatch(
    {
      installed_python_pkgs(
        envname = envname,
        conda = conda,
        verbose = verbose
      )
    },
    error = function(e) {
      log_message(
        "Failed to get installed packages: {.val {e$message}}",
        message_type = "warning",
        verbose = verbose
      )
    }
  )

  packages_installed <- stats::setNames(
    rep(FALSE, length(packages)),
    packages
  )
  if (is.null(all_installed) || is.null(all_installed$package)) {
    return(packages_installed)
  }

  requirements <- env_requirements(
    include_optional = TRUE
  )
  pkg_name_mapping <- requirements$package_aliases

  for (i in seq_along(packages)) {
    pkg <- packages[i]
    requirement <- parse_python_requirement(pkg, fallback_name = names(pkg))
    pkg_name <- requirement$name
    pkg_version <- requirement$version
    pkg_operator <- requirement$operator

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
        installed_version <- unique(all_installed$version[
          all_installed$package == found_pkg_name
        ])
        version_match <- version_satisfies_requirement(
          installed_version = installed_version,
          required_version = pkg_version,
          operator = pkg_operator
        )

        if (version_match) {
        log_message(
          "{.pkg {pkg_name}} {.pkg {pkg_operator}} {.pkg {pkg_version}}",
          message_type = "success",
          verbose = verbose
        )
      } else {
        log_message(
          "{.pkg {pkg_name}} found but version mismatch: installed {.pkg {paste(installed_version, collapse = ', ')}}, required {.pkg {pkg_operator}} {.pkg {pkg_version}}",
          message_type = "warning",
          verbose = verbose
        )
      }

        packages_installed[pkg] <- version_match
      } else {
        packages_installed[pkg] <- TRUE
        installed_version <- unique(all_installed$version[
          all_installed$package == found_pkg_name
        ])
        log_message(
          "{.pkg {pkg_name}} version: {.pkg {paste(installed_version, collapse = ', ')}}",
          message_type = "success",
          verbose = verbose
        )
      }
    } else {
      packages_installed[pkg] <- FALSE
      log_message(
        "{.pkg {pkg_name}} not found",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  return(packages_installed)
}

resolve_requested_python_packages <- function(packages) {
  requirements <- env_requirements(
    include_optional = TRUE
  )
  package_versions <- requirements$packages

  resolved <- vapply(
    packages,
    function(pkg) {
      requirement <- parse_python_requirement(pkg)
      if (
        !is.na(requirement$operator) ||
          !requirement$name %in% names(package_versions)
      ) {
        return(pkg)
      }
      unname(package_versions[[requirement$name]])
    },
    character(1),
    USE.NAMES = FALSE
  )

  stats::setNames(resolved, names(packages))
}

parse_python_requirement <- function(pkg, fallback_name = NULL) {
  if (grepl("^git\\+", pkg)) {
    pkg_name <- fallback_name %||%
      basename(strsplit(pkg, "@", fixed = TRUE)[[1]][1])
    pkg_name <- sub("\\.git$", "", pkg_name)
    return(list(
      name = pkg_name,
      operator = NA_character_,
      version = NA_character_
    ))
  }

  match <- regexec("^([^<>=!~]+?)\\s*([<>=!~]{1,2})\\s*(.+)$", pkg, perl = TRUE)
  parts <- regmatches(pkg, match)[[1]]
  if (length(parts) == 4) {
    return(list(
      name = fallback_name %||% trimws(parts[2]),
      operator = trimws(parts[3]),
      version = trimws(parts[4])
    ))
  }

  list(
    name = fallback_name %||% pkg,
    operator = NA_character_,
    version = NA_character_
  )
}

normalize_package_version <- function(version) {
  version <- trimws(as.character(version))
  sub("^v", "", version)
}

compare_versions_safe <- function(installed_version, required_version) {
  installed_version <- normalize_package_version(installed_version)
  required_version <- normalize_package_version(required_version)

  vapply(
    installed_version,
    function(version) {
      tryCatch(
        suppressWarnings(
          as.integer(utils::compareVersion(version, required_version))
        ),
        error = function(e) {
          if (identical(version, required_version)) {
            0L
          } else {
            NA_integer_
          }
        }
      )
    },
    integer(1)
  )
}

version_satisfies_requirement <- function(
  installed_version,
  required_version,
  operator
) {
  if (is.na(required_version) || is.na(operator) || !nzchar(required_version)) {
    return(TRUE)
  }

  comparison <- compare_versions_safe(installed_version, required_version)
  if (all(is.na(comparison))) {
    return(FALSE)
  }

  any(
    switch(operator,
      "==" = comparison == 0L,
      "=" = comparison == 0L,
      ">=" = comparison >= 0L,
      "<=" = comparison <= 0L,
      ">" = comparison > 0L,
      "<" = comparison < 0L,
      "!=" = comparison != 0L,
      "~=" = comparison >= 0L,
      FALSE
    ),
    na.rm = TRUE
  )
}

installed_python_pkgs <- function(
  envname = NULL,
  conda = "auto",
  verbose = TRUE
) {
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
    "Retrieving package list for environment: {.file {envname}}",
    verbose = verbose
  )

  tryCatch(
    {
      all_installed <- get_namespace_fun(
        "reticulate",
        "conda_list_packages"
      )(
        conda = conda,
        envname = envname,
        no_pip = FALSE
      )
      log_message(
        "Found {.val {nrow(all_installed)}} packages installed",
        verbose = verbose
      )
      return(all_installed)
    },
    error = function(e) {
      log_message(
        "Failed to retrieve package list: {.val {e$message}}",
        message_type = "error",
        verbose = verbose
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
  verbose = TRUE
) {
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
          envname,
          "'? (y/N): "
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

  python <- resolve_python_executable(
    envname = envname,
    conda = conda,
    error_message = "Failed to get Python path"
  )

  if (is.null(python)) {
    return(invisible(FALSE))
  }

  remove_many_python_packages <- function(packages, uv = NULL) {
    statuses <- vapply(
      packages,
      remove_single_python_package,
      logical(1),
      python = python,
      uv = uv,
      verbose = verbose
    )
    all(statuses)
  }

  if (pip) {
    result <- remove_pkg(
      packages,
      python,
      envname,
      conda,
      verbose
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
        "{.pkg {packages}} removal failed via {.pkg conda}, trying {.pkg uv} as fallback...",
        message_type = "warning",
        verbose = verbose
      )

      result <- tryCatch(
        {
          uv <- find_uv(
            python = python,
            envname = envname,
            conda = conda,
            auto_install = TRUE
          )

          if (is.null(uv)) {
            log_message(
              "{.pkg uv} not found and installation failed, falling back to {.pkg pip}",
              message_type = "warning"
            )
            remove_many_python_packages(packages, uv = NULL)
          } else {
            remove_many_python_packages(packages, uv = uv)
          }
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

remove_pkg <- function(
  packages,
  python,
  envname,
  conda,
  verbose
) {
  uv <- find_uv(
    python = python,
    envname = envname,
    conda = conda,
    auto_install = TRUE
  )

  if (!is.null(uv)) {
    log_message(
      "Removing {.pkg {packages}} via {.pkg uv}...",
      verbose = verbose
    )

    statuses <- vapply(
      packages,
      remove_single_python_package,
      logical(1),
      python = python,
      uv = uv,
      verbose = verbose
    )

    if (all(statuses)) {
      log_message(
        "{.pkg {packages}} removed successfully via {.pkg uv}",
        message_type = "success",
        verbose = verbose
      )
      return(TRUE)
    } else {
      log_message(
        "{.pkg uv} removal failed, falling back to {.pkg pip}...",
        message_type = "warning",
        verbose = verbose
      )
    }
  } else {
    log_message(
      "{.pkg uv} not found, falling back to {.pkg pip}...",
      message_type = "warning",
      verbose = verbose
    )
  }

  log_message(
    "Removing {.pkg {packages}} via {.pkg pip}...",
    verbose = verbose
  )
  statuses <- vapply(
    packages,
    remove_single_python_package,
    logical(1),
    python = python,
    uv = NULL,
    verbose = verbose
  )

  if (all(statuses)) {
    log_message(
      "{.pkg {packages}} removed successfully via {.pkg pip}",
      message_type = "success",
      verbose = verbose
    )
    return(TRUE)
  } else {
    log_message(
      "Failed to remove {.pkg {packages}} via {.pkg pip}",
      message_type = "warning",
      verbose = verbose
    )
    return(FALSE)
  }
}

remove_single_python_package <- function(
  pkg,
  python,
  uv = NULL,
  verbose = TRUE
) {
  if (!is.null(uv)) {
    log_message(
      "Removing {.pkg {pkg}} via {.pkg uv}...",
      verbose = verbose
    )

    args <- c("pip", "uninstall")
    if (uv != "python -m uv") {
      args <- c(args, "--python", python)
    }
    args <- c(args, pkg)

    status <- run_uv_command(uv, python, args)
    if (status == 0L) {
      log_message(
        "{.pkg {pkg}} removed successfully via {.pkg uv}",
        message_type = "success",
        verbose = verbose
      )
      return(TRUE)
    }

    log_message(
      "Failed to remove {.pkg {pkg}} via {.pkg uv} [error code {.val {status}}], trying {.pkg pip} as fallback...",
      message_type = "warning",
      verbose = verbose
    )
  }

  log_message(
    "Removing {.pkg {pkg}} via {.pkg pip}...",
    verbose = verbose
  )
  status <- run_pip_command(python, c("uninstall", "-y", pkg))

  if (status == 0L) {
    log_message(
      "{.pkg {pkg}} removed successfully via {.pkg pip}",
      message_type = "success",
      verbose = verbose
    )
    return(TRUE)
  }

  log_message(
    "{.pkg {pkg}} removal failed via {.pkg pip} [error code {.val {status}}]",
    message_type = "warning",
    verbose = verbose
  )
  FALSE
}
