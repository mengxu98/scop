#' Check and install python packages
#'
#' @param packages A character vector, indicating package names which should be installed or removed. Use \code{⁠<package>==<version>}⁠ to request the installation of a specific version of a package.
#' @param envname The name of a conda environment.
#' @param conda The path to a conda executable. Use \code{"auto"} to allow scop to automatically find an appropriate conda binary.
#' @param force Whether to force package installation. Default is \code{FALSE}.
#' @param pip Whether to use pip for package installation. By default, packages are installed from the active conda channels.
#' @param pip_options An optional character vector of additional command line arguments to be passed to \code{pip}. Only relevant when \code{pip = TRUE}.
#' @param ... Other arguments passed to \code{\link[reticulate]{conda_install}}
#'
#' @export
#'
#' @examples
#' check_python(packages = c("bbknn", "scanorama"))
#' \dontrun{
#' check_python(
#'   packages = "scvi-tools==0.20.0",
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
    warning(
      envname,
      " python environment does not exist. Create it with the PrepareEnv function...",
      immediate. = TRUE
    )
    PrepareEnv()
  }
  if (isTRUE(force)) {
    pkg_installed <- stats::setNames(rep(FALSE, length(packages)), packages)
    pip_options <- c(pip_options, "--force-reinstall")
  } else {
    pkg_installed <- exist_Python_pkgs(
      packages = packages,
      envname = envname,
      conda = conda
    )
  }
  if (sum(!pkg_installed) > 0) {
    pkgs_to_install <- names(pkg_installed)[!pkg_installed]
    message("Try to install ", paste0(pkgs_to_install, collapse = ","), " ...")
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
  }

  pkg_installed <- exist_Python_pkgs(
    packages = packages,
    envname = envname,
    conda = conda
  )
  if (sum(!pkg_installed) > 0) {
    stop(
      "Failed to install the package(s): ",
      paste0(names(pkg_installed)[!pkg_installed], collapse = ","),
      " into the environment \"",
      envname,
      "\". Please install manually."
    )
  } else {
    return(invisible(NULL))
  }
}

#' Check and install R packages
#'
#' @param packages Package to be installed. Package source can be CRAN, Bioconductor or Github, e.g. scmap, quadbiolab/simspec.
#' By default, the package name is extracted according to the \code{packages} parameter.
#' @param install_methods Functions used to install R packages.
#' @param lib The location of the library directories where to install the packages.
#' @param force Whether to force the installation of packages. Default is \code{FALSE}.
#'
#' @export
check_r <- function(
    packages,
    install_methods = c(
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
      message("Install package: \"", pkg_name, "\" ...")
      status_list[[pkg]] <- FALSE
      i <- 1
      while (isFALSE(status_list[[pkg]])) {
        tryCatch(
          expr = {
            if (grepl("BiocManager", install_methods[i])) {
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
    stop(
      "Failed to install the package(s): ",
      paste0(names(out), collapse = ","),
      ". Please install manually."
    )
  }
}
