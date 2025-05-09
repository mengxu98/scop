.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  env <- FALSE
  if (isTRUE(getOption("scop_env_init", default = TRUE))) {
    conda <- find_conda()
    if (!is.null(conda)) {
      envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
      env <- env_exist(conda = conda, envname = get_envname(), envs_dir = envs_dir)
      if (isFALSE(env)) {
        packageStartupMessage("scop python environment not found.")
      }
    } else {
      packageStartupMessage("Conda not found.")
    }
    if (isTRUE(env)) {
      status <- tryCatch(
        {
          Sys.unsetenv("RETICULATE_PYTHON")
          python_path <- conda_python(conda = conda)
          reticulate::use_python(python_path, required = TRUE)

          pyinfo <- utils::capture.output(reticulate::py_config())
          pyinfo_mesg <- c(
            "====================== scop conda environment ======================",
            paste0("conda:          ", conda),
            paste0("environment:    ", paste0(envs_dir, "/", get_envname())),
            "======================== scop python config ========================",
            pyinfo,
            "==================================================================="
          )
          invisible(lapply(pyinfo_mesg, packageStartupMessage))
          invisible(run_Python(command = "import matplotlib", envir = .GlobalEnv))
          if (!interactive()) {
            invisible(run_Python(command = "matplotlib.use('pdf')", envir = .GlobalEnv))
          }
          invisible(run_Python(command = "import matplotlib.pyplot as plt", envir = .GlobalEnv))
          invisible(run_Python(command = "import scanpy", envir = .GlobalEnv))
          packageStartupMessage("Conda path can be specified with the command `options(reticulate.conda_binary = \"/path/to/conda\")` before loading the package")
          packageStartupMessage("scop python environment can be disabled with the command `options(scop_env_init = FALSE)` before loading the package")
        },
        error = identity
      )
      if (inherits(status, "error")) {
        packageStartupMessage(status)
      }
    } else {
      packageStartupMessage("If you have already created an scop python environment using conda, you can specify the conda path by setting options(reticulate.conda_binary = \"/path/to/conda\", scop_env_name = \"scop_env\") before loading the package.")
    }
  }
}
