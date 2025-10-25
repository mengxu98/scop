# -*- coding: utf-8 -*-

#' @title Single-Cell Omics analysis Pipeline
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
#' Use [cli::ansi_strip] to get rid of the colors.
#' @md
#' @param unicode Whether to use Unicode symbols on UTF-8 platforms.
#' Default is [cli::is_utf8_output].
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
    "          0          1        2             3     4
                     _____ _________  ____
                    / ___// ___/ __ ./ __ .
                   (__  )/ /__/ /_/ / /_/ /
                  /____/ .___/.____/ .___/
                                  /_/
      5               6      7        8          9"
  )

  hexa <- c("*", ".", "o", "*", ".", "*", ".", "o", ".", "*")
  if (unicode) {
    hexa <- c("*" = "\u2b22", "o" = "\u2b21", "." = ".")[hexa]
  }

  cols <- c(
    "red", "yellow", "green", "magenta", "cyan",
    "yellow", "green", "white", "magenta", "cyan"
  )

  col_hexa <- mapply(
    function(x, y) cli::make_ansi_style(y)(x),
    hexa, cols,
    SIMPLIFY = FALSE
  )

  for (i in 0:9) {
    pat <- paste0("\\b", i, "\\b")
    logo <- sub(pat, col_hexa[[i + 1]], logo)
  }

  structure(cli::col_blue(logo), class = "scop_logo")
}

#' @title print scop logo
#'
#' @param x Input infromation.
#' @param ... Other parameters.
#'
#' @method print scop_logo
#'
#' @export
print.scop_logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  verbose <- thisutils::get_verbose()
  if (isTRUE(verbose)) {
    version <- utils::packageDescription(
      pkgname,
      fields = "Version"
    )
    scop_env_init <- getOption("scop_env_init", default = FALSE)
    version <- utils::packageDescription(pkgname, fields = "Version")
    msg <- paste0(
      cli::col_grey(strrep("-", 60)),
      "\n",
      cli::col_blue(pkgname, " version ", version),
      "\n"
    )
    if (isFALSE(scop_env_init)) {
      msg <- paste0(
        msg,
        "\n",
        cli::col_grey("Python environment initialization is disabled"),
        "\n",
        cli::col_grey("To enable it, set: options(scop_env_init = TRUE)"),
        "\n"
      )
    }
    suppress_msg <- paste0(
      cli::col_grey("The message can be suppressed by: "),
      "\n",
      cli::col_grey("  suppressPackageStartupMessages(library(scop))"),
      "\n",
      cli::col_grey("  or options(log_message.verbose = FALSE)")
    )
    if (isFALSE(scop_env_init)) {
      msg <- paste0(
        msg,
        "\n",
        suppress_msg,
        "\n",
        cli::col_grey(strrep("-", 60))
      )
    }

    packageStartupMessage(scop_logo())
    packageStartupMessage(msg)

    if (isTRUE(scop_env_init)) {
      tryCatch(
        {
          conda <- find_conda()
          if (is.null(conda)) {
            packageStartupMessage(
              cli::col_grey(
                "Conda not found. Run: PrepareEnv() to create the environment"
              )
            )
            return(invisible(NULL))
          }
          envname <- get_envname()
          envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
          env <- env_exist(
            conda = conda,
            envname = envname,
            envs_dir = envs_dir
          )

          if (isFALSE(env)) {
            packageStartupMessage(
              cli::col_grey(
                "Python environment not found. Run: PrepareEnv() to create the environment"
              )
            )
            return(invisible(NULL))
          }
          set_python_env(conda = conda, envname = envname, verbose = FALSE)

          packageStartupMessage(
            cli::col_green("conda environment initialized successfully")
          )

          env_info(conda = conda, envname = envname)

          packageStartupMessage(
            "\n",
            cli::col_grey(
              "Configure conda path: options(reticulate.conda_binary = \"/path/to/conda\")"
            ),
            "\n",
            cli::col_grey(
              "Disable Python initialization information: options(scop_env_init = FALSE)"
            ),
            "\n\n",
            suppress_msg,
            "\n",
            cli::col_grey(strrep("-", 60))
          )
        },
        error = function(e) {
          packageStartupMessage(
            cli::col_grey(
              "Failed to initialize Python environment: ",
              e$message,
              "\n",
              "Run: PrepareEnv() to set up the environment, or disable: options(scop_env_init = FALSE)",
              "\n\n",
              suppress_msg,
              "\n",
              strrep("-", 60)
            )
          )
        }
      )
    }
  }
}
