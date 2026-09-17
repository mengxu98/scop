# -*- coding: utf-8 -*-

#' @title Single-cell and Spatial omics analysis pipeline
#'
#' @description
#' An R package for single-cell and spatial omics analysis, integration, visualization, and interactive exploration.
#'
#' @author Meng xu (Maintainer), \email{mengxu98@qq.com}
#'
#' @source \url{https://github.com/mengxu98/scop}
#'
#' @md
#' @docType package
#' @name scop-package
#' @useDynLib scop, .registration = TRUE
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
#' @rdname scop_logo
#' @export
#' @examples
#' scop_logo()
scop_logo <- function(
  unicode = cli::is_utf8_output()
) {
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
#' @param x Input logo object.
#' @param ... Additional arguments passed to [cat()].
#'
#' @rdname scop_logo
#' @method print scop_logo
#'
#' @export
print.scop_logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  options(scop_env_cache = NULL)

  if (
    isTRUE(thisutils::get_verbose()) &&
      isTRUE(getOption("scop_env_init", default = FALSE))
  ) {
    tryCatch(
      {
        conda <- find_conda()
        if (is.null(conda)) {
          packageStartupMessage(
            cli::col_grey(
              "Conda-compatible environment manager not found. Run: PrepareEnv() to create the environment"
            )
          )
          return(invisible(NULL))
        }
        envname <- get_envname()
        envs_dir <- get_conda_envs_dir(conda = conda)
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
        python_path <- conda_python(
          conda = conda,
          envname = envname
        )
        configure_python_runtime(python_path)

        packageStartupMessage(cli::col_green("Python environment initialized"))

        env_info(conda = conda, envname = envname)
      },
      error = function(e) {
        packageStartupMessage(
          cli::col_grey(
            "Failed to initialize Python environment: ",
            e$message
          )
        )
      }
    )
  }
}
