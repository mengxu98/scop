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

#' @title print scop logo
#'
#' @param x Input infromation.
#' @param ... Other parameters.
#'
#' @method print scop
#'
#' @export
print.scop <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  options(memory.limit = Inf)

  scop_env_init <- getOption("scop_env_init", default = FALSE)
  version <- utils::packageDescription(pkgname, fields = "Version")
  msg <- cli::col_blue(pkgname, " version ", version)
  if (!isTRUE(scop_env_init)) {
    msg <- paste0(
      msg,
      "\n",
      cli::col_grey("Python environment initialization is disabled."),
      "\n",
      cli::col_grey("To enable it, set: options(scop_env_init = TRUE)")
    )
  }
  msg <- paste0(
    msg,
    "\n",
    cli::col_grey("This message can be suppressed by: "),
    "\n",
    cli::col_grey("  suppressPackageStartupMessages(library(scop))")
  )
  msg <- paste0(
    strrep("-", 60),
    "\n",
    msg
  )
  if (!isTRUE(scop_env_init)) {
    msg <- paste0(
      msg,
      "\n",
      strrep("-", 60)
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
            cli::col_grey("Conda not found. Run: PrepareEnv() to create the environment")
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

        if (!isTRUE(env)) {
          packageStartupMessage(
            cli::col_grey(
              "Python environment not found. Run: PrepareEnv() to create the environment"
            )
          )
          return(invisible(NULL))
        }
        set_python_env(conda = conda, envname = envname)

        pyinfo <- utils::capture.output(reticulate::py_config())
        packageStartupMessage(
          cli::col_grey("Conda environment initialized successfully")
        )

        pyinfo_mesg <- c(
          cli::col_blue("conda environment: "),
          cli::col_grey(paste0("  conda:          ", conda)),
          cli::col_grey(paste0("  environment:    ", paste0(envs_dir, "/", get_envname()))),
          cli::col_blue("python config: "),
          cli::col_grey(paste0("  ", pyinfo))
        )
        invisible(lapply(pyinfo_mesg, packageStartupMessage))

        packageStartupMessage(
          "\n",
          cli::col_grey(
            "Configure conda path: options(reticulate.conda_binary = \"/path/to/conda\")"
          ),
          "\n",
          cli::col_grey(
            "Disable Python initialization information: options(scop_env_init = FALSE)"
          ),
          "\n",
          strrep("-", 60)
        )
      },
      error = function(e) {
        packageStartupMessage(
          cli::col_grey(
            "Failed to initialize Python environment: ",
            e$message,
            "\n",
            "Run: PrepareEnv() to set up the environment, or disable: options(scop_env_init = FALSE)",
            "\n",
            strrep("-", 60)
          )
        )
      }
    )
  }
}
