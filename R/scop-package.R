# -*- coding: utf-8 -*-

#' @title scop Single-Cell Omics analysis Pipeline
#'
#' @description
#' An end-to-end Single-Cell Omics analysis Pipeline designed to facilitate comprehensive analysis and exploration of single-cell omics data.
#' 
#' @author Meng xu (Maintainer), \email{mengxu98@qq.com}
#'
#' @source \url{https://github.com/mengxu98/scop}
#'
#' @md
#' @docType package
#' @name scop-package
"_PACKAGE"

#' @title scop logo
#'
#' @description
#' The scop logo, using ASCII or Unicode characters
#' Use [cli::ansi_strip()] to get rid of the colors.
#' @md
#' @param unicode Whether to use Unicode symbols. Default is `TRUE` on UTF-8 platforms.
#'
#' @references
#'  \url{https://github.com/tidyverse/tidyverse/blob/main/R/logo.R}
#'
#' @export
#' @examples
#' scop_logo()
scop_logo <- function(
    unicode = cli::is_utf8_output()) {
  logo <- c(
    "       0        1      2           3    4
       ______________  ____
      / ___/ ___/ __ ./ __ .
     (__  ) /__/ /_/ / /_/ /
    /____/.___/.____/ .___/
                   /_/
    5             6      7      8       9   "
  )

  hexa <- c("*", ".", "o", "*", ".", "*", ".", "o", ".", "*")
  if (unicode) {
    hexa <- c("*" = "\u2b22", "o" = "\u2b21", "." = ".")[hexa]
  }

  cols <- c(
    "red", "yellow", "green", "magenta", "cyan",
    "yellow", "green", "white", "magenta", "cyan"
  )

  col_hexa <- purrr::map2(
    hexa, cols, ~ cli::make_ansi_style(.y)(.x)
  )

  for (i in 0:9) {
    pat <- paste0("\\b", i, "\\b")
    logo <- sub(pat, col_hexa[[i + 1]], logo)
  }

  structure(cli::col_blue(logo), class = "logo")
}

#' @title print logo
#'
#' @param x Input infromation.
#' @param ... Other parameters.
#'
#' @method print logo
#'
#' @export
print.logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0(
    "-------------------------------------------------------
",
    cli::col_blue(" ", pkgname, " version ", version),
    "
   This message can be suppressed by:
     suppressPackageStartupMessages(library(scop))
-------------------------------------------------------"
  )

  packageStartupMessage(scop_logo())
  packageStartupMessage(msg)

  # options(future.globals.maxSize = Inf)
  # env <- FALSE
  # if (isTRUE(getOption("scop_env_init", default = TRUE))) {
  #   conda <- find_conda()
  #   if (!is.null(conda)) {
  #     envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
  #     env <- env_exist(conda = conda, envname = get_envname(), envs_dir = envs_dir)
  #     if (isFALSE(env)) {
  #       packageStartupMessage("scop python environment not found.")
  #     }
  #   } else {
  #     packageStartupMessage("Conda not found.")
  #   }
  #   if (isTRUE(env)) {
  #     status <- tryCatch(
  #       {
  #         Sys.unsetenv("RETICULATE_PYTHON")
  #         python_path <- conda_python(conda = conda)
  #         reticulate::use_python(python_path, required = TRUE)

  #         pyinfo <- utils::capture.output(reticulate::py_config())
  #         pyinfo_mesg <- c(
  #           "====================== scop conda environment ======================",
  #           paste0("conda:          ", conda),
  #           paste0("environment:    ", paste0(envs_dir, "/", get_envname())),
  #           "======================== scop python config ========================",
  #           pyinfo,
  #           "==================================================================="
  #         )
  #         invisible(lapply(pyinfo_mesg, packageStartupMessage))
  #         invisible(run_python(command = "import matplotlib", envir = .GlobalEnv))
  #         if (!interactive()) {
  #           invisible(run_python(command = "matplotlib.use('pdf')", envir = .GlobalEnv))
  #         }
  #         invisible(run_python(command = "import matplotlib.pyplot as plt", envir = .GlobalEnv))
  #         invisible(run_python(command = "import scanpy", envir = .GlobalEnv))
  #         packageStartupMessage("Conda path can be specified with the command `options(reticulate.conda_binary = \"/path/to/conda\")` before loading the package")
  #         packageStartupMessage("scop python environment can be disabled with the command `options(scop_env_init = FALSE)` before loading the package")
  #       },
  #       error = identity
  #     )
  #     if (inherits(status, "error")) {
  #       packageStartupMessage(status)
  #     }
  #   } else {
  #     packageStartupMessage("If you have already created an scop python environment using conda, you can specify the conda path by setting options(reticulate.conda_binary = \"/path/to/conda\", scop_env_name = \"scop_env\") before loading the package.")
  #   }
  # }
}
