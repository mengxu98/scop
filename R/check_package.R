#' @title Check and install python packages
#'
#' @md
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
      message_type = "warning"
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
      "Try to install: {.pkg {pkgs_to_install}}"
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
      message_type = "warning"
    )
  } else {
    return(invisible(NULL))
  }
}

#' @title Check and install R packages
#'
#' @md
#' @param packages Package to be installed.
#' Package source can be CRAN, Bioconductor or Github.
#' By default, the package name is extracted according to the `packages` parameter.
#' @param install_methods Functions used to install R packages.
#' @param lib The location of the library directories where to install the packages.
#' @param force Whether to force the installation of packages.
#' Default is FALSE.
#'
#' @export
check_r <- function(
    packages,
    install_methods = c(
      "pak::pak",
      "BiocManager::install",
      "install.packages",
      "devtools::install_github"
    ),
    lib = .libPaths()[1],
    force = FALSE) {
  status_list <- list()
  for (pkg in packages) {
    version <- NULL
    if (grepl("/", pkg)) {
      pkg_name <- strsplit(pkg, split = "/|@|==", perl = TRUE)[[1]][[2]]
    } else {
      pkg_name <- strsplit(pkg, split = "@|==", perl = TRUE)[[1]][[1]]
      if (length(strsplit(pkg, split = "@|==", perl = TRUE)[[1]]) > 1) {
        version <- strsplit(pkg, split = "@|==", perl = TRUE)[[1]][[2]]
      }
    }
    dest <- gsub("@.*|==.*|>=.*", "", pkg)
    if (is.null(version)) {
      force_update <- isTRUE(force)
    } else {
      force_update <- isTRUE(
        utils::packageVersion(pkg_name) < package_version(version)
      ) ||
        isTRUE(force)
    }
    if (
      !suppressPackageStartupMessages(
        requireNamespace(
          pkg_name,
          quietly = TRUE
        )
      ) ||
        isTRUE(force_update)
    ) {
      log_message(
        "Installing package: {.pkg {pkg_name}}..."
      )
      status_list[[pkg]] <- FALSE
      i <- 1
      while (isFALSE(status_list[[pkg]])) {
        tryCatch(
          expr = {
            if (grepl("pak::pak", install_methods[i])) {
              if (!requireNamespace("pak", quietly = TRUE)) {
                utils::install.packages("pak", lib = lib)
              }
              if (!requireNamespace("withr", quietly = TRUE)) {
                utils::install.packages("withr", lib = lib)
              }
              eval(
                str2lang(
                  paste0(
                    "withr::with_libpaths(new = \"",
                    lib,
                    "\", ",
                    install_methods[i],
                    "(\"",
                    dest,
                    "\", ask = FALSE))"
                  )
                )
              )
            } else if (grepl("BiocManager", install_methods[i])) {
              if (!requireNamespace("BiocManager", quietly = TRUE)) {
                utils::install.packages("BiocManager", lib = lib)
              }
              eval(
                str2lang(
                  paste0(
                    install_methods[i],
                    "(\"",
                    dest,
                    "\", lib=\"",
                    lib,
                    "\", update = FALSE, upgrade = \"never\", ask = FALSE, force = TRUE)"
                  )
                )
              )
            } else if (grepl("devtools", install_methods[i])) {
              if (!requireNamespace("devtools", quietly = TRUE)) {
                utils::install.packages("devtools", lib = lib)
              }
              if (!requireNamespace("withr", quietly = TRUE)) {
                utils::install.packages("withr", lib = lib)
              }
              eval(
                str2lang(
                  paste0(
                    "withr::with_libpaths(new = \"",
                    lib,
                    "\", ",
                    install_methods[i],
                    "(\"",
                    dest,
                    "\", upgrade = \"never\", force = TRUE))"
                  )
                )
              )
            } else {
              eval(
                str2lang(
                  paste0(
                    install_methods[i],
                    "(\"",
                    dest,
                    "\", lib=\"",
                    lib,
                    "\", force = TRUE)"
                  )
                )
              )
            }
          },
          error = function(e) {
            status_list[[pkg]] <- FALSE
          }
        )
        status_list[[pkg]] <- requireNamespace(pkg_name, quietly = TRUE)
        i <- i + 1
        if (i > length(install_methods)) {
          break
        }
      }
    } else {
      status_list[[pkg]] <- TRUE
    }
  }
  out <- sapply(status_list, isTRUE)
  out <- out[!out]
  if (length(out) > 0) {
    log_message(
      "Failed to install: {.pkg {names(out)}}. Please install manually",
      message_type = "error"
    )
  }
}
