#' @title Prepare the python environment
#'
#' @description
#' Prepare the python environment by installing the required dependencies and setting up the environment.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param envname The name of the conda-compatible Python environment.
#' If `NULL`, the environment name will be set to `"scop_env"`.
#' Default is `NULL`.
#' @param conda The path or command name of a conda-compatible executable
#' (`conda`, `mamba`, or `micromamba`). Use `"auto"` to allow automatically
#' finding an appropriate environment manager. If `"micromamba"` is requested
#' and micromamba is not available on `PATH`, a package-managed micromamba is
#' downloaded automatically.
#' @param miniconda_repo Repository URL for miniconda.
#' Default is \url{https://repo.anaconda.com/miniconda}.
#' @param force Whether to force recreation of the environment.
#' If `TRUE`, the existing environment will be removed and recreated.
#' Default is `FALSE`.
#' @param version The Python version.
#' Default is `"3.10-1"` on macOS and Unix and `"3.11-1"` on Windows.
#' @param modules Optional Python dependency modules to install in addition to
#' the default scientific stack.
#' Supported values are `"scanpy"`, `"scvi"`, `"glue"`, `"scanorama"`, `"bbknn"`,
#' `"celltypist"`, `"cellphonedb"`, `"magic"`, `"scrublet"`,
#' `"sccoda"`, `"doubletdetection"`, `"doublet"`, `"palantir"`, `"scvelo"`,
#' `"cellrank"`, `"wot"`, `"phate"`, `"pacmap"`, `"trimap"`, `"multimap"`,
#' `"scomm"`, `"scenic"`, `"seacells"`, `"tage"`,
#' `"scmalignantfinder"`, `"secact"`, `"scpagwas"`, and
#' `"external_wrappers"`.
#' If `NULL` or omitted in [PrepareEnv()], the default environment is installed.
#' The default excludes `"sccoda"` and `"scomm"` because their TensorFlow stacks
#' are not compatible with the default JAX/scVI stack in the same environment;
#' request them explicitly for scCODA/scOMM workflows. `"scenic"` is also
#' excluded from the default environment and is prepared in `"scenic_env"` by
#' default because SCENIC requires an older Python/numpy stack. On Windows, the
#' default also excludes `"scvi"`, `"glue"`, and `"multimap"` because those
#' upstream stacks are more reliable when requested explicitly for
#' method-specific workflows.
#' @param pip_options Additional command line arguments to be passed to `uv`/`pip` when installing pip packages.
#' @param ... Additional arguments passed to package installation functions.
#'
#' @export
PrepareEnv <- function(
  envname = NULL,
  conda = "auto",
  miniconda_repo = "https://repo.anaconda.com/miniconda",
  version = if (is_windows()) "3.11-1" else "3.10-1",
  force = FALSE,
  modules = NULL,
  pip_options = character(),
  verbose = TRUE,
  ...
) {
  modules <- normalize_env_modules(modules = modules)
  if ("scenic" %in% modules) {
    scenic_allowed <- c("scenic", "regdiffusion")
    if (length(setdiff(modules, scenic_allowed)) > 0) {
      log_message(
        "{.arg modules = 'scenic'} must be prepared as a standalone environment, optionally with {.val regdiffusion}. Run {.code PrepareEnv(envname = 'scenic_env', modules = c('scenic', 'regdiffusion'))}.",
        message_type = "error"
      )
    }
    if (is.null(envname)) {
      envname <- "scenic_env"
    }
    log_message(
      "{.pkg SCENIC} pins {.pkg numpy} to {.val 1.23.5}. Prepare it in an isolated environment such as {.val {envname}}; installing it into an environment shared with {.pkg scanpy}, {.pkg scvi}, or {.pkg scvelo} may downgrade dependencies and break those workflows.",
      message_type = "warning"
    )
    if (!identical(version, "3.10-1")) {
      log_message(
        "{.pkg SCENIC} requires Python 3.10 in {.pkg scop}; using {.val 3.10-1} for {.arg modules = 'scenic'}.",
        message_type = "warning"
      )
      version <- "3.10-1"
    }
  }
  envname <- get_envname(envname)
  pip_options <- normalize_cli_args(pip_options)
  requirements <- env_requirements(
    version = version,
    modules = modules
  )
  python_version <- requirements[["python"]]

  conda_auto <- identical(conda, "auto")
  if (!is.null(conda)) {
    conda <- resolve_conda(conda)
  }

  cache_spec <- build_env_cache_spec(
    envname = envname,
    python_version = python_version,
    conda = conda,
    modules = modules,
    pip_options = pip_options,
    requirements = requirements
  )
  if (is_cached_env_valid(cache_spec) && isFALSE(force)) {
    python_cached <- getOption("scop_env_cache", default = NULL)[["python"]] %||%
      tryCatch(
        conda_python(envname = envname, conda = conda),
        error = function(...) NULL
      )
    if (!is.null(python_cached) && nzchar(python_cached) && file.exists(python_cached)) {
      assert_python_runtime_switchable(
        python_cached,
        restart_hint = python_runtime_restart_hint(envname = envname, modules = modules)
      )
      configure_python_runtime(python_cached)
    }
    remember_python_environment(envname = envname, conda = conda)
    if ("scomm" %in% modules) {
      ensure_scomm_runtime_support(
        envname = envname,
        conda = conda,
        pip_options = pip_options,
        keep_jax = "scvi" %in% modules
      )
    }
    ensure_external_wrapper_r_packages(
      modules = modules,
      verbose = verbose
    )
    log_message(
      "{cli::col_green('Python environment completed')}\n",
      "{cli::col_grey('Until next loading, the environment will be cached')}",
      message_type = "success"
    )
    return(invisible(NULL))
  }

  log_message(
    "Preparing python environment...",
    text_color = "blue",
    message_type = "running"
  )

  log_message(
    "Environment name: {.file {envname}} and python version: {.pkg {python_version}}"
  )

  if (isTRUE(conda_auto)) {
    log_message("Auto-detecting conda-compatible environment manager...")
  }
  if (!is.null(conda)) {
    log_message(
      "Using {.pkg {conda_manager_label(conda)}} executable: {.file {conda}}"
    )
  }

  if (is.null(conda)) {
    env <- FALSE
    log_message(
      "Conda-compatible environment manager not found, will install miniconda",
      message_type = "warning"
    )
  } else {
    envs_dir <- get_conda_envs_dir(conda = conda)
    env_path <- conda_env_path(
      envname = envname,
      conda = conda,
      envs_dir = envs_dir
    )
    if (!is.null(env_path)) {
      assert_python_runtime_switchable(
        conda_env_python_path(env_path),
        restart_hint = python_runtime_restart_hint(envname = envname, modules = modules)
      )
    }
    env <- env_exist(
      conda = conda,
      envname = envname,
      envs_dir = envs_dir
    )

    if (isTRUE(force) && isTRUE(env)) {
      log_message("Force recreating environment...")
      unlink(env_path, recursive = TRUE)
      env <- FALSE
    }

    if (isTRUE(env)) {
      log_message(
        "Using existing environment: {.file {env_path}}"
      )
    }
  }

  if (isFALSE(env)) {
    force <- TRUE

    if (is.null(conda)) {
      log_message("Installing miniconda...")
      conda <- install_miniconda2(miniconda_repo)
      log_message(
        "Using {.pkg {conda_manager_label(conda)}} executable: {.file {conda}}"
      )
      cache_spec <- build_env_cache_spec(
        envname = envname,
        python_version = python_version,
        conda = conda,
        modules = modules,
        pip_options = pip_options,
        requirements = requirements
      )
    }

    log_message(
      "Creating Python environment with python {.val {python_version}} using {.pkg {conda_manager_label(conda)}}..."
    )

    accept_conda_tos(conda = conda)

    base_packages <- c("pip", "setuptools<81", "wheel")

    python_path <- create_conda_env(
      conda = conda,
      envname = envname,
      python_version = python_version,
      packages = base_packages
    )

    python <- conda_python(envname = envname, conda = conda)

    env_path <- conda_env_path(envname = envname, conda = conda) %||%
      dirname(dirname(python_path))
    env <- file.exists(env_path)

    if (isFALSE(env)) {
      log_message(
        "Environment creation failed",
        message_type = "warning"
      )
      log_message(
        "Conda info: {.val {conda_info_json(conda = conda)}}",
        verbose = verbose
      )
      log_message(
        "Conda environments: {.val {reticulate::conda_list(conda = conda)}}",
        verbose = verbose
      )
      log_message(
        "Unable to find environment under the expected path: {.file {env_path}}\n",
        "manager: {.pkg {conda_manager_label(conda)}}\n",
        "executable: {.file {conda}}\n",
        "python: {.file {python_path}}",
        message_type = "error"
      )
    } else {
      log_message(
        "Environment created successfully: {.file {env_path}}",
        message_type = "success"
      )
    }
  }

  install_methods <- requirements[["install_methods"]]
  packages <- requirements[["packages"]]

  python <- conda_python(envname = envname, conda = conda)
  uv <- find_uv(
    python = python,
    envname = envname,
    conda = conda,
    auto_install = TRUE,
    pip_options = pip_options
  )

  conda_packages <- packages[install_methods == "conda"]
  pip_packages <- packages[install_methods == "pip"]

  log_message(
    "{.val {length(packages)}} package{?s} to install, {.val {length(conda_packages)}} conda packages and {.val {length(pip_packages)}} pip packages"
  )

  install_ok <- TRUE
  if (length(conda_packages) > 0) {
    log_message(
      "Installing {.val {length(conda_packages)}} {.pkg conda} packages"
    )
    install_ok <- isTRUE(check_python(
      packages = conda_packages,
      envname = envname,
      conda = conda,
      force = force,
      pip = FALSE,
      verbose = FALSE,
      ...
    )) && install_ok
  }

  if (length(pip_packages) > 0) {
    if (!is.null(uv)) {
      log_message(
        "Installing {.val {length(pip_packages)}} packages using {.pkg uv}"
      )
    } else {
      log_message(
        "Installing {.val {length(pip_packages)}} {.pkg pip} packages"
      )
    }
    install_ok <- isTRUE(check_python(
      packages = pip_packages,
      envname = envname,
      conda = conda,
      force = force,
      pip = TRUE,
      pip_options = pip_options,
      verbose = FALSE,
      ...
    )) && install_ok
  }

  if ("scvi" %in% modules) {
    ensure_windows_scvi_support(
      envname = envname,
      conda = conda,
      force = force,
      pip_options = pip_options
    )
  }

  if ("scomm" %in% modules) {
    ensure_scomm_runtime_support(
      envname = envname,
      conda = conda,
      pip_options = pip_options,
      keep_jax = "scvi" %in% modules
    )
  }
  ensure_external_wrapper_r_packages(
    modules = modules,
    verbose = verbose
  )

  assert_python_runtime_switchable(
    python,
    restart_hint = python_runtime_restart_hint(envname = envname, modules = modules)
  )
  configure_python_runtime(python)
  remember_python_environment(envname = envname, conda = conda)

  if (isTRUE(install_ok)) {
    log_message(
      "{cli::col_green('Python environment ready')}\n",
      "{cli::col_grey('Until next loading, the environment will be cached')}",
      message_type = "success"
    )
  } else {
    log_message(
      "Python environment setup completed with missing packages. Review the warnings above before using modules that require them.",
      message_type = "warning"
    )
  }

  env_info(
    conda = conda,
    envname = envname,
    verbose = isTRUE(getOption("log_message.verbose", TRUE))
  )

  set_env_cache(cache_spec, python)
}

normalize_cli_args <- function(args) {
  if (length(args) == 0) {
    return(character())
  }

  args <- unlist(args, use.names = FALSE)
  if (length(args) == 1 && grepl("\\s", args)) {
    args <- strsplit(args, "\\s+")[[1]]
  }

  args <- trimws(as.character(args))
  args[nzchar(args)]
}

supported_env_modules <- function() {
  c(
    "scanpy",
    "scvi",
    "glue",
    "scanorama",
    "bbknn",
    "celltypist",
    "cellphonedb",
    "magic",
    "scfea",
    "scrublet",
    "sccoda",
    "doubletdetection",
    "palantir",
    "scvelo",
    "cellrank",
    "wot",
    "phate",
    "pacmap",
    "trimap",
    "multimap",
    "scomm",
    "scenic",
    "regdiffusion",
    "scenicplus",
    "seacells",
    "tage",
    "scmalignantfinder",
    "secact",
    "scpagwas",
    "external_wrappers"
  )
}

default_env_modules <- function() {
  excluded <- c(
    "scfea",
    "sccoda",
    "scomm",
    "scenic",
    "regdiffusion",
    "scenicplus",
    "scmalignantfinder",
    "secact",
    "scpagwas",
    "external_wrappers"
  )
  if (is_windows()) {
    excluded <- c(excluded, "scvi", "glue", "multimap")
  }
  setdiff(supported_env_modules(), excluded)
}

optional_env_modules <- function() {
  c(
    "celltypist",
    "cellphonedb",
    "magic",
    "scrublet",
    "doubletdetection"
  )
}

env_module_dependencies <- function() {
  list(
    celltypist = "scanpy",
    cellphonedb = "scanpy",
    sccoda = "scanpy",
    palantir = "scanpy",
    scvelo = "scanpy",
    cellrank = c("scanpy", "scvelo"),
    wot = "scanpy",
    scvi = "scanpy",
    glue = "scanpy",
    scanorama = "scanpy",
    bbknn = "scanpy",
    multimap = "scanpy",
    scmalignantfinder = "scanpy"
  )
}

env_module_requirements <- function() {
  list(
    scanpy = scanpy_python_requirements(),
    scvi = scvi_python_requirements(),
    glue = glue_python_requirements(),
    scanorama = scanorama_python_requirements(),
    bbknn = bbknn_python_requirements(),
    celltypist = celltypist_python_requirements(),
    cellphonedb = cellphonedb_python_requirements(),
    magic = magic_python_requirements(),
    scfea = scfea_python_requirements(),
    scrublet = scrublet_python_requirements(),
    sccoda = sccoda_python_requirements(),
    doubletdetection = doubletdetection_python_requirements(),
    palantir = palantir_python_requirements(),
    scvelo = scvelo_python_requirements(),
    cellrank = cellrank_python_requirements(),
    wot = wot_python_requirements(),
    phate = phate_python_requirements(),
    pacmap = pacmap_python_requirements(),
    trimap = trimap_python_requirements(),
    multimap = multimap_python_requirements(),
    scomm = scomm_python_requirements(),
    scenic = scenic_python_requirements(),
    regdiffusion = regdiffusion_python_requirements(),
    scenicplus = scenicplus_python_requirements(),
    seacells = seacells_python_requirements(),
    tage = tage_python_requirements(),
    scmalignantfinder = scmalignantfinder_python_requirements()
  )
}

normalize_env_modules <- function(modules = NULL, include_optional = FALSE) {
  modules <- modules %||% character(0)
  modules <- unlist(modules, use.names = FALSE)
  modules <- trimws(as.character(modules))
  modules <- modules[nzchar(modules)]

  if (length(modules) == 0) {
    modules <- default_env_modules()
  }

  if (isTRUE(include_optional)) {
    modules <- c(modules, optional_env_modules())
  }

  if ("optional" %in% modules) {
    modules <- c(setdiff(modules, "optional"), optional_env_modules())
  }

  if ("doublet" %in% modules) {
    modules <- c(
      setdiff(modules, "doublet"),
      "scrublet",
      "doubletdetection"
    )
  }

  if ("external_wrappers" %in% modules) {
    modules <- c(
      setdiff(modules, "external_wrappers"),
      "scmalignantfinder",
      "secact",
      "scpagwas"
    )
  }

  module_dependencies <- env_module_dependencies()
  for (module in modules) {
    deps <- module_dependencies[[module]]
    if (!is.null(deps)) {
      modules <- c(modules, deps)
    }
  }

  modules <- unique(modules)
  invalid_modules <- setdiff(modules, c(supported_env_modules(), "optional", "doublet"))
  if (length(invalid_modules) > 0) {
    log_message(
      "{.arg modules} contains unsupported values: {.val {invalid_modules}}",
      message_type = "error"
    )
  }
  modules
}

ensure_external_wrapper_r_packages <- function(modules, verbose = TRUE) {
  repos <- c(
    secact = "data2intelligence/SecAct",
    scpagwas = "sulab-wmu/scPagwas"
  )
  modules <- intersect(names(repos), modules)
  if (length(modules) == 0) {
    return(invisible(TRUE))
  }

  ok <- TRUE
  for (module in modules) {
    repo <- unname(repos[[module]])
    installed <- isTRUE(check_r(repo, dependencies = NA, verbose = verbose))
    ok <- installed && ok
  }

  invisible(ok)
}

build_env_cache_spec <- function(
  envname,
  python_version,
  conda,
  modules,
  pip_options,
  requirements = NULL
) {
  conda_path <- if (is.null(conda) || identical(conda, "")) {
    NULL
  } else {
    normalizePath(conda, mustWork = FALSE)
  }

  list(
    envname = get_envname(envname),
    python_version = python_version,
    conda = conda_path,
    modules = normalize_env_modules(modules = modules),
    pip_options = normalize_cli_args(pip_options),
    requirements = requirements[c("packages", "install_methods")]
  )
}

clear_env_cache <- function() {
  options(scop_env_cache = NULL)
  invisible(NULL)
}

remember_python_environment <- function(envname, conda) {
  opts <- list(scop_envname = envname)
  if (!is.null(conda)) {
    opts[["scop_conda"]] <- conda
    opts[["reticulate.conda_binary"]] <- conda
  }
  do.call(options, opts)
  invisible(NULL)
}

set_env_cache <- function(spec, python = NULL) {
  cache <- spec
  if (!is.null(python) && nzchar(python) && file.exists(python)) {
    cache[["python"]] <- normalizePath(python, mustWork = FALSE)
  }
  options(scop_env_cache = cache)
  invisible(cache)
}

is_cached_env_valid <- function(spec) {
  cache <- getOption("scop_env_cache", default = NULL)
  required_fields <- c(
    "envname",
    "python_version",
    "conda",
    "modules",
    "pip_options",
    "requirements"
  )

  if (!is.list(cache) || !all(required_fields %in% names(cache))) {
    return(FALSE)
  }

  if (
    !identical(unname(cache[required_fields]), unname(spec[required_fields]))
  ) {
    return(FALSE)
  }

  if (is.null(spec[["conda"]]) || !file.exists(spec[["conda"]])) {
    return(FALSE)
  }

  if (!env_exist(conda = spec[["conda"]], envname = spec[["envname"]])) {
    return(FALSE)
  }

  python <- cache[["python"]] %||%
    tryCatch(
      conda_python(envname = spec[["envname"]], conda = spec[["conda"]]),
      error = function(e) NULL
    )

  !is.null(python) && file.exists(python)
}

prepend_path_var <- function(var, values) {
  values <- values[nzchar(values)]
  values <- values[file.exists(values)]
  if (length(values) == 0) {
    return(invisible(NULL))
  }

  values <- unique(normalizePath(values, mustWork = FALSE))
  current <- Sys.getenv(var, unset = "")
  current_values <- strsplit(current, .Platform$path.sep, fixed = TRUE)[[1]]
  current_values <- current_values[nzchar(current_values)]
  do.call(
    Sys.setenv,
    stats::setNames(
      list(paste(unique(c(values, current_values)), collapse = .Platform$path.sep)),
      var
    )
  )
  invisible(NULL)
}

find_conda_shared_library <- function(env_path, names) {
  library_dirs <- file.path(
    env_path,
    c("lib", "lib64", file.path("Library", "lib"), file.path("Library", "bin"))
  )
  candidates <- as.vector(outer(library_dirs, names, file.path))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    return(NULL)
  }

  normalizePath(candidates[[1]], mustWork = FALSE)
}

configure_linux_cpp_runtime <- function(env_path) {
  if (isFALSE(is_linux())) {
    return(invisible(FALSE))
  }

  libgcc <- find_conda_shared_library(env_path, c("libgcc_s.so.1", "libgcc_s.so"))
  libstdcpp <- find_conda_shared_library(env_path, c("libstdc++.so.6", "libstdc++.so"))
  libcrypto <- find_conda_shared_library(env_path, c("libcrypto.so.3", "libcrypto.so"))
  libssl <- find_conda_shared_library(env_path, c("libssl.so.3", "libssl.so"))
  runtime_libs <- c(libcrypto, libssl, libgcc, libstdcpp)
  runtime_libs <- runtime_libs[nzchar(runtime_libs)]
  if (length(runtime_libs) == 0) {
    return(invisible(FALSE))
  }

  prepend_path_var("LD_PRELOAD", runtime_libs)

  for (lib in runtime_libs) {
    tryCatch(
      dyn.load(lib, local = FALSE, now = TRUE),
      error = function(...) NULL
    )
  }

  invisible(TRUE)
}

normalize_python_runtime_path <- function(path) {
  path <- as.character(path %||% "")
  if (length(path) == 0 || !nzchar(path[[1]])) {
    return(NULL)
  }
  normalizePath(path[[1]], winslash = "/", mustWork = FALSE)
}

same_python_runtime <- function(x, y) {
  x <- normalize_python_runtime_path(x)
  y <- normalize_python_runtime_path(y)
  if (is.null(x) || is.null(y)) {
    return(FALSE)
  }
  if (identical(x, y)) {
    return(TRUE)
  }

  x_base <- tolower(basename(x))
  y_base <- tolower(basename(y))
  identical(dirname(x), dirname(y)) &&
    startsWith(x_base, "python") &&
    startsWith(y_base, "python")
}

active_python_runtime <- function() {
  if (isFALSE(reticulate::py_available(initialize = FALSE))) {
    return(NULL)
  }

  config <- tryCatch(
    reticulate::py_config(),
    error = function(...) NULL
  )
  normalize_python_runtime_path(config[["python"]])
}

assert_python_runtime_switchable <- function(python_path, restart_hint = NULL) {
  active_python <- active_python_runtime()
  if (is.null(active_python) || same_python_runtime(active_python, python_path)) {
    return(invisible(TRUE))
  }

  if (is.null(restart_hint)) {
    restart_hint <- "Restart R, then run {.code PrepareEnv(...)} and the downstream analysis with the same {.arg envname} and {.arg conda}."
  }

  log_message(
    "Python is already initialized with {.file {active_python}} and cannot be switched to {.file {python_path}} in the same R session. {restart_hint}",
    message_type = "error"
  )
}

scenic_runtime_restart_hint <- function(envname = "scenic_env") {
  paste0(
    "Restart R, then run ",
    "PrepareEnv(envname = \"", envname, "\", modules = \"scenic\") ",
    "before RunSCENIC()."
  )
}

python_runtime_restart_hint <- function(envname = "scop_env", modules = NULL) {
  if ("scenic" %in% modules) {
    return(scenic_runtime_restart_hint(envname = envname))
  }
  if ("cellphonedb" %in% modules) {
    return(paste0(
      "Restart R, then run ",
      "PrepareEnv(envname = \"", envname, "\", modules = \"cellphonedb\") ",
      "before RunCellphoneDB()."
    ))
  }
  module_code <- if (!is.null(modules) && length(modules) > 0L) {
    paste0(
      ", modules = ",
      paste(deparse(modules), collapse = "")
    )
  } else {
    ""
  }
  paste0(
    "Restart R, then run ",
    "PrepareEnv(envname = \"", envname, "\"", module_code, ") ",
    "before the downstream Python-backed analysis."
  )
}

configure_python_thread_env <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
  Sys.setenv(NUMEXPR_NUM_THREADS = "1")
  Sys.setenv(KMP_WARNINGS = "0")
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  Sys.setenv(NUMBA_NUM_THREADS = "1")
}

configure_python_runtime <- function(python_path) {
  python_path <- normalizePath(python_path, mustWork = FALSE)
  if (!nzchar(python_path) || !file.exists(python_path)) {
    return(invisible(FALSE))
  }

  configure_python_thread_env()
  assert_python_runtime_switchable(python_path)

  python_dir <- dirname(python_path)
  env_path <- dirname(python_dir)

  Sys.setenv(
    RETICULATE_PYTHON = python_path,
    PYTHONNOUSERSITE = "1",
    PIP_USER = "0",
    MPLBACKEND = "Agg"
  )
  Sys.unsetenv("PYTHONPATH")

  if (isFALSE(reticulate::py_available(initialize = FALSE))) {
    tryCatch(
      reticulate::use_python(python_path, required = TRUE),
      error = function(...) NULL
    )
  }

  prepend_path_var(
    "PATH",
    c(
      python_dir,
      file.path(env_path, "bin"),
      file.path(env_path, "Library", "bin"),
      file.path(env_path, "Scripts")
    )
  )

  if (is_linux()) {
    prepend_path_var(
      "LD_LIBRARY_PATH",
      c(
        file.path(env_path, "lib"),
        file.path(env_path, "lib64"),
        file.path(env_path, "Library", "lib")
      )
    )
    configure_linux_cpp_runtime(env_path)
  }

  if (is_osx()) {
    prepend_path_var(
      "DYLD_FALLBACK_LIBRARY_PATH",
      c(
        file.path(env_path, "lib"),
        file.path(env_path, "lib64")
      )
    )
  }

  invisible(TRUE)
}

set_python_env <- function(conda, envname, verbose = TRUE) {
  Sys.unsetenv("RETICULATE_PYTHON")
  options(reticulate.miniconda.enabled = FALSE)

  configure_python_thread_env()

  python_path <- conda_python(
    conda = conda,
    envname = envname
  )

  configure_python_runtime(python_path)

  reticulate::use_python(
    python_path,
    required = TRUE
  )

  tryCatch(
    {
      reticulate::py_run_string(
        "
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['KMP_WARNINGS'] = '0'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
os.environ['NUMBA_NUM_THREADS'] = '1'
try:
    import numba
    if hasattr(numba, 'set_num_threads'):
        try:
            numba.set_num_threads(1)
        except RuntimeError:
            pass
    if hasattr(numba, 'config'):
        numba.config.NUMBA_NUM_THREADS = 1
except ImportError:
    pass
except Exception:
    pass
"
      )
    },
    error = function(e) {}
  )
}

ensure_scomm_runtime_support <- function(
  envname,
  conda,
  pip_options = character(),
  keep_jax = FALSE,
  verbose = TRUE
) {
  installed <- tryCatch(
    installed_python_pkgs(envname = envname, conda = conda),
    error = function(...) NULL
  )
  if (is.null(installed) || is.null(installed$package)) {
    return(invisible(FALSE))
  }

  removable <- if (isTRUE(keep_jax)) {
    character()
  } else {
    intersect(c("jax", "jaxlib"), installed$package)
  }
  if (length(removable) > 0) {
    log_message(
      "Removing conflicting packages for {.pkg scOMM}: {.pkg {removable}}",
      verbose = verbose
    )
    remove_python(
      packages = removable,
      envname = envname,
      conda = conda,
      pip = TRUE,
      force = TRUE,
      verbose = verbose
    )
  }

  check_python(
    packages = if (isTRUE(keep_jax)) {
      "ml-dtypes>=0.3.2"
    } else {
      "ml-dtypes==0.3.2"
    },
    envname = envname,
    conda = conda,
    force = FALSE,
    pip = TRUE,
    pip_options = pip_options,
    verbose = verbose
  )

  invisible(TRUE)
}

install_miniconda2 <- function(
  miniconda_repo,
  timeout = 600
) {
  log_message("Installing miniconda...")
  options(timeout = timeout)

  info <- as.list(Sys.info())

  if (is_apple_silicon()) {
    base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
    name <- "Miniforge3-MacOSX-arm64.sh"
    url <- file.path(base, name)
  } else {
    version <- "3"
    arch <- get_namespace_fun(
      "reticulate",
      "miniconda_installer_arch"
    )(info)

    name <- if (is_windows()) {
      sprintf("Miniconda%s-latest-Windows-%s.exe", version, arch)
    } else if (is_osx()) {
      sprintf("Miniconda%s-latest-MacOSX-%s.sh", version, arch)
    } else if (is_linux()) {
      sprintf("Miniconda%s-latest-Linux-%s.sh", version, arch)
    } else {
      log_message(
        "Unsupported platform: {.val {info$sysname}}",
        message_type = "error"
      )
    }

    url <- file.path(miniconda_repo, name)
  }

  options(reticulate.miniconda.url = url)

  if (!is.na(Sys.getenv("USER", unset = NA))) {
    miniconda_path <- gsub(
      pattern = "\\$USER",
      replacement = Sys.getenv("USER"),
      reticulate::miniconda_path()
    )
  } else {
    miniconda_path <- reticulate::miniconda_path()
  }

  if (dir.exists(miniconda_path)) {
    log_message("Removing existing miniconda installation...")
    unlink(miniconda_path, recursive = TRUE)
  }

  log_message("Downloading and installing miniconda from: {.url {url}}")
  reticulate::install_miniconda(
    path = miniconda_path,
    force = TRUE,
    update = FALSE
  )

  conda <- get_namespace_fun(
    "reticulate",
    "miniconda_conda"
  )(miniconda_path)
  log_message("Miniconda installed at: {.file {miniconda_path}}")

  conda
}

micromamba_platform <- function() {
  machine <- tolower(Sys.info()[["machine"]])
  if (is_windows()) {
    if (machine %in% c("x86_64", "amd64")) {
      return("win-64")
    }
  } else if (is_osx()) {
    if (machine %in% c("arm64", "aarch64")) {
      return("osx-arm64")
    }
    if (machine %in% c("x86_64", "amd64")) {
      return("osx-64")
    }
  } else if (is_linux()) {
    if (machine %in% c("aarch64", "arm64")) {
      return("linux-aarch64")
    }
    if (machine %in% c("x86_64", "amd64")) {
      return("linux-64")
    }
  }

  log_message(
    "Unsupported platform for automatic micromamba installation: {.val {Sys.info()[['sysname']]}} {.val {machine}}",
    message_type = "error"
  )
}

managed_micromamba_dir <- function() {
  file.path(
    tools::R_user_dir("scop", which = "cache"),
    "micromamba",
    micromamba_platform()
  )
}

managed_micromamba_root <- function() {
  file.path(
    tools::R_user_dir("scop", which = "cache"),
    "micromamba-root"
  )
}

legacy_managed_micromamba_roots <- function() {
  file.path(
    tools::R_user_dir("scop", which = "data"),
    "micromamba-root"
  )
}

find_managed_micromamba <- function() {
  install_dir <- managed_micromamba_dir()
  candidates <- file.path(
    install_dir,
    c(
      if (is_windows()) "Library/bin/micromamba.exe" else "bin/micromamba",
      if (is_windows()) "micromamba.exe" else "micromamba"
    )
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    return(NULL)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

configure_managed_micromamba_root <- function() {
  root <- managed_micromamba_root()
  current <- Sys.getenv("MAMBA_ROOT_PREFIX", unset = "")
  current_norm <- normalize_conda_paths(current)
  managed_roots <- normalize_conda_paths(c(root, legacy_managed_micromamba_roots()))

  if (
    !nzchar(current) ||
      length(current_norm) == 0 ||
      current_norm %in% managed_roots ||
      grepl("\\s", current_norm)
  ) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(MAMBA_ROOT_PREFIX = root)
  }
  invisible(Sys.getenv("MAMBA_ROOT_PREFIX"))
}

install_micromamba <- function(timeout = 600) {
  existing <- find_managed_micromamba()
  if (!is.null(existing)) {
    configure_managed_micromamba_root()
    return(existing)
  }

  platform <- micromamba_platform()
  install_dir <- managed_micromamba_dir()
  archive <- file.path(tempdir(), sprintf("micromamba-%s.tar.bz2", platform))
  url <- sprintf(
    "https://micro.mamba.pm/api/micromamba/%s/latest",
    platform
  )

  log_message(
    "Installing micromamba for platform {.val {platform}}..."
  )
  log_message(
    "Downloading micromamba from: {.url {url}}"
  )

  dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)
  options(timeout = timeout)
  status <- tryCatch(
    utils::download.file(url, archive, mode = "wb", quiet = FALSE),
    error = function(e) e
  )
  if (inherits(status, "error") || !identical(status, 0L)) {
    msg <- if (inherits(status, "error")) status$message else status
    log_message(
      "Failed to download micromamba: {.val {msg}}",
      message_type = "error"
    )
  }

  unlink(install_dir, recursive = TRUE, force = TRUE)
  dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(archive, exdir = install_dir)

  micromamba <- find_managed_micromamba()
  if (is.null(micromamba)) {
    log_message(
      "Micromamba installation completed but executable was not found under {.file {install_dir}}",
      message_type = "error"
    )
  }

  Sys.chmod(micromamba, mode = "0755")
  configure_managed_micromamba_root()
  log_message(
    "Micromamba installed at: {.file {micromamba}}"
  )

  micromamba
}

#' @title Print environment information
#' @inheritParams PrepareEnv
#' @param verbose Whether to print environment information.
#' @md
env_info <- function(conda, envname, verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  manager <- conda_manager_label(conda = conda)
  env_path <- conda_env_path(envname = envname, conda = conda)
  envs_dir <- if (!is.null(env_path)) {
    dirname(env_path)
  } else {
    get_conda_envs_dir(conda = conda)
  }

  python_path <- tryCatch(
    conda_python(envname = envname, conda = conda),
    error = function(...) NULL
  )
  env_display <- env_path %||% if (!is.null(envs_dir)) {
    file.path(envs_dir, get_envname(envname))
  } else {
    get_envname(envname)
  }
  python_path <- python_path %||% if (!is.null(env_display) && nzchar(env_display)) {
    file.path(
      env_display,
      if (is_windows()) "python.exe" else "bin/python"
    )
  } else {
    NULL
  }
  py_version <- if (!is.null(python_path) && nzchar(python_path) && file.exists(python_path)) {
    out <- tryCatch(
      suppressWarnings(system2(python_path, "--version", stdout = TRUE, stderr = TRUE)),
      error = function(...) character()
    )
    out[[1]] %||% "unknown"
  } else {
    "unknown"
  }

  py_info_mesg <- c(
    cli::col_blue(
      "Environment manager config:"
    ),
    cli::col_grey(
      " Manager:       ",
      manager
    ),
    cli::col_grey(
      " Executable:    ",
      conda
    ),
    cli::col_grey(
      " Environment:   ",
      env_display
    ),
    cli::col_blue(
      "Python config:"
    ),
    cli::col_grey(
      " python:        ",
      python_path
    ),
    cli::col_grey(
      " version:       ",
      py_version
    )
  )
  invisible(lapply(py_info_mesg, packageStartupMessage))
}

#' @title Python environment requirements
#'
#' @description
#' The function returns a list of requirements including the required Python version,
#' package versions, and package name aliases for platform-specific packages.
#' All packages will be installed using uv as the primary tool.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param version The Python version of the environment.
#' Default is `"3.10-1"`.
#' @param include_optional Whether to include optional Python dependencies.
#' @param modules Optional requirement modules to include. Supported values are
#' `"scanpy"`, `"scvi"`, `"scanorama"`, `"bbknn"`, `"celltypist"`,
#' `"cellphonedb"`, `"magic"`, `"scrublet"`, `"doubletdetection"`,
#' `"sccoda"`, `"doublet"`, `"palantir"`, `"scvelo"`, `"cellrank"`, `"wot"`,
#' `"phate"`, `"pacmap"`, `"trimap"`, `"multimap"`,
#' `"scomm"`, `"scenic"`, `"seacells"`, `"tage"`,
#' `"scmalignantfinder"`, `"secact"`, `"scpagwas"`, and
#' `"external_wrappers"`. If `NULL`, the default environment is returned. The default
#' excludes `"sccoda"`, `"scomm"`, and `"scenic"` because these workflows
#' require dependency stacks that should be prepared explicitly. The
#' `"scenic"` module is standalone and always uses Python `"3.10-1"`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{python}{Python version string}
#'   \item{packages}{Named vector of package version specifications}
#'   \item{package_aliases}{Named list mapping logical package names to actual installed names}
#' }
#'
#' @export
#' @examples
#' env_requirements("3.10-1")
env_requirements <- function(
  version = "3.10-1",
  include_optional = FALSE,
  modules = NULL,
  verbose = TRUE
) {
  version <- match.arg(
    version,
    choices = c("3.10-1", "3.11-1")
  )
  modules <- normalize_env_modules(
    modules = modules,
    include_optional = include_optional
  )
  if ("scenic" %in% modules) {
    scenic_allowed <- c("scenic", "regdiffusion")
    if (length(setdiff(modules, scenic_allowed)) > 0) {
      log_message(
        "{.arg modules = 'scenic'} must be used as a standalone environment module, optionally with {.val regdiffusion}.",
        message_type = "error"
      )
    }
    version <- "3.10-1"
  }

  base_requirements <- if ("scenic" %in% modules) {
    scenic_core_python_requirements()
  } else {
    core_python_requirements()
  }
  package_install_methods <- base_requirements$install_methods
  package_versions <- base_requirements$packages
  package_aliases <- base_requirements$package_aliases

  module_requirements <- env_module_requirements()

  for (module in modules) {
    req_i <- module_requirements[[module]]
    if (is.null(req_i)) {
      next
    }
    package_install_methods <- c(package_install_methods, req_i$install_methods)
    package_versions <- c(package_versions, req_i$packages)
    package_aliases <- c(package_aliases, req_i$package_aliases)
  }

  if (all(c("scvi", "scomm") %in% modules) && "ml_dtypes" %in% names(package_versions)) {
    package_versions[["ml_dtypes"]] <- "ml-dtypes>=0.3.2"
  }
  if ("scenic" %in% modules && "numpy" %in% names(package_versions)) {
    package_versions[["numpy"]] <- "numpy==1.23.5"
  }

  requirements <- list(
    python = version,
    packages = package_versions,
    install_methods = package_install_methods,
    package_aliases = package_aliases
  )

  return(requirements)
}

core_python_requirements <- function() {
  tbb_install_method <- if (is_apple_silicon()) {
    "pip"
  } else {
    "conda"
  }
  list(
    packages = c(
      "setuptools" = "setuptools<81",
      "leidenalg" = "leidenalg==0.10.2",
      "tbb" = if (is_apple_silicon()) {
        "pxr-tbb==2022.2.0.0"
      } else {
        "tbb==2022.2.0"
      },
      "python-igraph" = "python-igraph==0.11.9",
      "matplotlib" = "matplotlib==3.10.8",
      "numba" = "numba==0.59.1",
      "llvmlite" = "llvmlite==0.42.0",
      "numpy" = "numpy==1.26.4",
      "packaging" = "packaging>=24.0",
      "pandas" = "pandas==2.2.0",
      "scikit-learn" = "scikit-learn==1.7.0",
      "scipy" = "scipy==1.15.3",
      if (is_apple_silicon()) {
        c("llvm-openmp" = "llvm-openmp>=17")
      } else {
        NULL
      }
    ),
    install_methods = c(
      "setuptools" = "pip",
      "leidenalg" = "conda",
      "tbb" = tbb_install_method,
      "python-igraph" = if (is_apple_silicon()) {
        "pip"
      } else {
        "conda"
      },
      "matplotlib" = "pip",
      "numba" = "pip",
      "llvmlite" = "pip",
      "numpy" = "pip",
      "packaging" = "pip",
      "pandas" = "pip",
      "scikit-learn" = "pip",
      "scipy" = "pip",
      if (is_apple_silicon()) {
        c("llvm-openmp" = "conda")
      } else {
        NULL
      }
    ),
    package_aliases = c(
      list("python-igraph" = "igraph"),
      if (is_apple_silicon()) {
        list("tbb" = "pxr-tbb")
      } else {
        list()
      }
    )
  )
}

scenic_core_python_requirements <- function() {
  list(
    packages = c(
      "setuptools" = "setuptools<81"
    ),
    install_methods = c(
      "setuptools" = "pip"
    ),
    package_aliases = list()
  )
}

scenic_backend_package <- function() {
  paste0("py", "scenic")
}

scenic_backend_requirement <- function(version = "0.12.1") {
  paste0(scenic_backend_package(), "==", version)
}

scanpy_python_requirements <- function() {
  list(
    packages = c(
      "scanpy" = "scanpy==1.11.3",
      "loompy" = "loompy"
    ),
    install_methods = c(
      "scanpy" = "pip",
      "loompy" = "pip"
    ),
    package_aliases = list()
  )
}

scvi_python_requirements <- function() {
  scvi_install_method <- if (is_windows()) "pip" else "conda"
  list(
    packages = c(
      "scvi-tools" = "scvi-tools==1.2.1",
      "jax" = "jax[cpu]==0.4.38"
    ),
    install_methods = c(
      "scvi-tools" = scvi_install_method,
      "jax" = "pip"
    ),
    package_aliases = list()
  )
}

glue_python_requirements <- function() {
  list(
    packages = c(
      "scglue" = "scglue==0.4.0",
      "bedtools" = "bedtools",
      "zlib" = "zlib"
    ),
    install_methods = c(
      "scglue" = "pip",
      "bedtools" = "conda",
      "zlib" = "conda"
    ),
    package_aliases = list()
  )
}

scanorama_python_requirements <- function() {
  list(
    packages = c(
      "scanorama" = "scanorama==1.7.4"
    ),
    install_methods = c(
      "scanorama" = "pip"
    ),
    package_aliases = list()
  )
}

bbknn_python_requirements <- function() {
  list(
    packages = c(
      "bbknn" = "bbknn==1.6.0"
    ),
    install_methods = c(
      "bbknn" = "pip"
    ),
    package_aliases = list()
  )
}

celltypist_python_requirements <- function() {
  list(
    packages = c(
      "celltypist" = "celltypist==1.7.1"
    ),
    install_methods = c(
      "celltypist" = "pip"
    ),
    package_aliases = list()
  )
}

cellphonedb_python_requirements <- function() {
  list(
    packages = c(
      "cellphonedb" = "cellphonedb==5.0.1"
    ),
    install_methods = c(
      "cellphonedb" = "pip"
    ),
    package_aliases = list()
  )
}

magic_python_requirements <- function() {
  list(
    packages = c(
      "magic-impute" = "magic-impute==3.0.0"
    ),
    install_methods = c(
      "magic-impute" = "pip"
    ),
    package_aliases = list()
  )
}

scfea_python_requirements <- function() {
  list(
    packages = c(
      "torch" = "torch",
      "numpy" = "numpy",
      "pandas" = "pandas",
      "tqdm" = "tqdm"
    ),
    install_methods = c(
      "torch" = "pip",
      "numpy" = "pip",
      "pandas" = "pip",
      "tqdm" = "pip"
    ),
    package_aliases = list()
  )
}

tage_python_requirements <- function() {
  list(
    packages = c(
      "joblib" = "joblib",
      "pandas" = "pandas",
      "scikit-learn" = "scikit-learn"
    ),
    install_methods = c(
      "joblib" = "pip",
      "pandas" = "pip",
      "scikit-learn" = "pip"
    ),
    package_aliases = list()
  )
}

scrublet_python_requirements <- function() {
  list(
    packages = c(
      "scrublet" = "scrublet==0.2.3"
    ),
    install_methods = c(
      "scrublet" = "pip"
    ),
    package_aliases = list()
  )
}

sccoda_python_requirements <- function() {
  list(
    packages = c(
      "tensorflow" = "tensorflow==2.16.2",
      "tensorflow-probability" = "tensorflow-probability==0.24.0",
      "tf-keras" = "tf-keras==2.16.0",
      "sccoda" = "sccoda==0.1.9"
    ),
    install_methods = c(
      "tensorflow" = "pip",
      "tensorflow-probability" = "pip",
      "tf-keras" = "pip",
      "sccoda" = "pip"
    ),
    package_aliases = list(
      "sccoda" = "scCODA",
      "tf-keras" = "tf_keras"
    )
  )
}

doubletdetection_python_requirements <- function() {
  list(
    packages = c(
      "doubletdetection" = "doubletdetection==4.3.0.post1",
      "louvain" = "louvain==0.8.2"
    ),
    install_methods = c(
      "doubletdetection" = "pip",
      "louvain" = "pip"
    ),
    package_aliases = list()
  )
}

palantir_python_requirements <- function() {
  list(
    packages = c(
      "palantir" = "palantir==1.4.1"
    ),
    install_methods = c(
      "palantir" = "pip"
    ),
    package_aliases = list()
  )
}

scvelo_python_requirements <- function() {
  list(
    packages = c(
      "scvelo" = "scvelo==0.3.3"
    ),
    install_methods = c(
      "scvelo" = "pip"
    ),
    package_aliases = list()
  )
}

seacells_python_requirements <- function() {
  list(
    packages = c(
      "SEACells" = "SEACells"
    ),
    install_methods = c(
      "SEACells" = "pip"
    ),
    package_aliases = list()
  )
}

scmalignantfinder_python_requirements <- function() {
  packages <- c(
    "scMalignantFinder" = "git+https://github.com/Jonyyqn/scMalignantFinder.git",
    "xgboost" = "xgboost"
  )
  install_methods <- c(
    "scMalignantFinder" = "pip",
    "xgboost" = "pip"
  )
  if (is_osx()) {
    packages <- c(packages, "libcxx" = "libcxx")
    install_methods <- c(install_methods, "libcxx" = "conda")
  }

  list(
    packages = packages,
    install_methods = install_methods,
    package_aliases = list()
  )
}

scenic_python_requirements <- function() {
  scenic_backend <- scenic_backend_package()
  scenic_package <- stats::setNames(scenic_backend_requirement(), scenic_backend)
  scenic_install_method <- stats::setNames("pip", scenic_backend)
  list(
    packages = c(
      scenic_package,
      "arboreto" = "arboreto==0.1.6",
      "ctxcore" = "ctxcore==0.2.0",
      "numpy" = "numpy==1.23.5",
      "dask" = "dask==2024.2.1",
      "distributed" = "distributed==2024.2.1",
      "pyarrow" = "pyarrow"
    ),
    install_methods = c(
      scenic_install_method,
      "arboreto" = "pip",
      "ctxcore" = "pip",
      "numpy" = "pip",
      "dask" = "pip",
      "distributed" = "pip",
      "pyarrow" = "pip"
    ),
    package_aliases = list()
  )
}

regdiffusion_python_requirements <- function() {
  list(
    packages = c(
      "regdiffusion" = "regdiffusion"
    ),
    install_methods = c(
      "regdiffusion" = "pip"
    ),
    package_aliases = list()
  )
}

scenicplus_python_requirements <- function() {
  list(
    packages = c(
      "scenicplus" = "scenicplus @ git+https://github.com/aertslab/scenicplus.git"
    ),
    install_methods = c(
      "scenicplus" = "pip"
    ),
    package_aliases = list()
  )
}

cellrank_python_requirements <- function() {
  list(
    packages = c(
      "cellrank" = "cellrank==2.0.7"
    ),
    install_methods = c(
      "cellrank" = "pip"
    ),
    package_aliases = list()
  )
}

wot_python_requirements <- function() {
  list(
    packages = c(
      "wot" = "wot==1.0.8.post2"
    ),
    install_methods = c(
      "wot" = "pip"
    ),
    package_aliases = list()
  )
}

phate_python_requirements <- function() {
  list(
    packages = c(
      "phate" = "phate==1.0.11"
    ),
    install_methods = c(
      "phate" = "pip"
    ),
    package_aliases = list()
  )
}

pacmap_python_requirements <- function() {
  list(
    packages = c(
      "pacmap" = "pacmap==0.8.0"
    ),
    install_methods = c(
      "pacmap" = "pip"
    ),
    package_aliases = list()
  )
}

trimap_python_requirements <- function() {
  list(
    packages = c(
      "trimap" = "trimap==1.1.4"
    ),
    install_methods = c(
      "trimap" = "pip"
    ),
    package_aliases = list()
  )
}

multimap_python_requirements <- function() {
  list(
    packages = c(
      "multimap" = "git+https://github.com/Teichlab/MultiMAP.git"
    ),
    install_methods = c(
      "multimap" = "pip"
    ),
    package_aliases = list(
      "multimap" = "MultiMAP",
      "MultiMAP" = "multimap"
    )
  )
}

scomm_python_requirements <- function() {
  list(
    packages = c(
      "tensorflow" = "tensorflow==2.16.2",
      "keras" = "keras==3.3.3",
      "tf_keras" = "tf_keras==2.16.0",
      "ml_dtypes" = "ml-dtypes==0.3.2"
    ),
    install_methods = c(
      "tensorflow" = "pip",
      "keras" = "pip",
      "tf_keras" = "pip",
      "ml_dtypes" = "pip"
    ),
    package_aliases = list(
      "tf_keras" = "tf-keras",
      "tf-keras" = "tf_keras",
      "ml_dtypes" = "ml-dtypes",
      "ml-dtypes" = "ml_dtypes"
    )
  )
}

env_exist <- function(
  conda = "auto",
  envname = NULL,
  envs_dir = NULL
) {
  envname <- get_envname(envname)
  conda <- resolve_conda(conda)

  if (!is.null(conda)) {
    env_paths <- conda_env_paths(
      envname = envname,
      conda = conda,
      envs_dir = envs_dir
    )

    if (length(env_paths) > 0) {
      existing_paths <- env_paths[file.exists(env_paths)]
      valid_paths <- existing_paths[valid_conda_env_paths(existing_paths)]
      if (length(valid_paths) > 0) {
        return(TRUE)
      }
      if (length(existing_paths) > 0) {
        log_message(
          "Environment directory exists but appears invalid: {.file {existing_paths[1]}}",
          message_type = "warning"
        )
        return(FALSE)
      }
    }

    envs_dirs <- envs_dir %||% get_conda_envs_dirs(conda = conda)
    envs <- tryCatch(
      reticulate::conda_list(conda = conda),
      error = function(...) NULL
    )
    envs <- filter_conda_env_table(
      conda_envs = envs,
      envs_dirs = envs_dirs,
      root_prefix = conda_root_prefix(conda = conda)
    )
    if (
      is.data.frame(envs) &&
        nrow(envs) > 0 &&
        any(envs$name == envname & file.exists(envs$python))
    ) {
      return(TRUE)
    }

    return(FALSE)
  } else {
    return(FALSE)
  }
}

get_envname <- function(envname = NULL) {
  if (is.character(envname)) {
    envname <- envname
  } else {
    envname <- getOption("scop_envname", default = "scop_env")
  }
  return(envname)
}

normalize_conda_paths <- function(paths, must_exist = FALSE) {
  paths <- unlist(paths, use.names = FALSE)
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0) {
    return(character())
  }
  paths <- path.expand(paths)
  if (isTRUE(must_exist)) {
    paths <- paths[dir.exists(paths)]
  }
  if (length(paths) == 0) {
    return(character())
  }
  unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
}

normalize_conda_paths_vector <- function(paths) {
  paths <- as.character(paths)
  invalid <- is.na(paths) | !nzchar(paths)
  paths <- path.expand(paths)
  paths[!invalid] <- normalizePath(
    paths[!invalid],
    winslash = "/",
    mustWork = FALSE
  )
  paths[invalid] <- NA_character_
  paths
}

conda_env_python_suffix <- function() {
  if (is_windows()) "python.exe" else "bin/python"
}

conda_env_python_path <- function(env_path) {
  file.path(env_path, conda_env_python_suffix())
}

valid_conda_env_paths <- function(env_paths) {
  env_paths <- normalize_conda_paths_vector(env_paths)
  if (length(env_paths) == 0) {
    return(logical())
  }

  !is.na(env_paths) &
    dir.exists(file.path(env_paths, "conda-meta")) &
    file.exists(conda_env_python_path(env_paths))
}

conda_env_path_in_dirs <- function(env_paths, envs_dirs, root_prefix = NULL) {
  env_paths <- normalize_conda_paths_vector(env_paths)
  envs_dirs <- normalize_conda_paths(envs_dirs)
  root_prefix <- normalize_conda_paths(root_prefix)

  if (length(env_paths) == 0) {
    return(logical())
  }

  keep <- rep(FALSE, length(env_paths))
  if (length(envs_dirs) > 0) {
    keep <- keep | dirname(env_paths) %in% envs_dirs
  }
  if (length(root_prefix) > 0) {
    keep <- keep | env_paths %in% root_prefix
  }

  keep
}

filter_conda_env_table <- function(conda_envs, envs_dirs = NULL, root_prefix = NULL) {
  if (!is.data.frame(conda_envs) || nrow(conda_envs) == 0) {
    return(conda_envs)
  }
  if (!all(c("name", "python") %in% names(conda_envs))) {
    return(conda_envs[0, , drop = FALSE])
  }

  env_paths <- normalize_conda_paths_vector(dirname(dirname(conda_envs$python)))
  keep <- file.exists(conda_envs$python) & valid_conda_env_paths(env_paths)

  scoped <- conda_env_path_in_dirs(
    env_paths = env_paths,
    envs_dirs = envs_dirs,
    root_prefix = root_prefix
  )
  if (any(scoped)) {
    keep <- keep & scoped
  }

  conda_envs[keep, , drop = FALSE]
}

parse_conda_json_output <- function(output) {
  if (is.null(output) || length(output) == 0) {
    return(NULL)
  }
  output <- output[nzchar(output)]
  if (length(output) == 0) {
    return(NULL)
  }
  text <- paste(output, collapse = "\n")
  json_start <- regexpr("[\\[{]", text, perl = TRUE)[1]
  if (is.na(json_start) || json_start < 1) {
    return(NULL)
  }
  text <- substring(text, json_start)

  tryCatch(
    get_namespace_fun("jsonlite", "fromJSON")(text, simplifyVector = TRUE),
    error = function(...) NULL
  )
}

run_conda_json <- function(conda, args) {
  output <- tryCatch(
    suppressWarnings(system2(conda, args, stdout = TRUE, stderr = FALSE)),
    error = function(...) NULL
  )
  status <- attr(output, "status") %||% 0L
  if (!identical(status, 0L)) {
    return(NULL)
  }
  parse_conda_json_output(output)
}

conda_info_json <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (is.null(conda)) {
    return(NULL)
  }

  info <- run_conda_json(conda, c("info", "--json"))
  if (!is.null(info)) {
    return(info)
  }

  info <- tryCatch(
    get_namespace_fun("reticulate", "conda_info")(conda = conda),
    error = function(...) NULL
  )
  if (!is.null(info)) {
    return(info)
  }

  NULL
}

conda_env_list_json <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (is.null(conda)) {
    return(NULL)
  }
  run_conda_json(conda, c("env", "list", "--json"))
}

conda_manager_type <- function(conda = "auto") {
  if (identical(conda, "auto")) {
    conda <- resolve_conda(conda)
  }
  if (is.null(conda)) {
    return("conda")
  }
  exe <- tolower(basename(as.character(conda[[1]])))
  if (startsWith(exe, "micromamba")) {
    return("micromamba")
  }
  if (startsWith(exe, "mamba")) {
    return("mamba")
  }
  "conda"
}

conda_manager_label <- function(conda = "auto") {
  conda_manager_type(conda = conda)
}

conda_root_prefix <- function(conda = "auto", info = NULL) {
  conda <- resolve_conda(conda)
  if (is.null(conda)) {
    return(NULL)
  }

  info <- info %||% conda_info_json(conda = conda)
  root <- NULL
  if (!is.null(info)) {
    root <- info[["root_prefix"]] %||%
      info[["root prefix"]] %||%
      info[["base environment"]] %||%
      info[["base_prefix"]]
  }

  if (is.null(root) && identical(conda_manager_type(conda), "micromamba")) {
    root_env <- Sys.getenv("MAMBA_ROOT_PREFIX", unset = "")
    if (nzchar(root_env)) {
      root <- root_env
    }
  }

  if (is.null(root) && file.exists(conda)) {
    root <- dirname(dirname(conda))
  }

  root <- normalize_conda_paths(root)
  if (length(root) == 0) {
    return(NULL)
  }
  root[[1]]
}

get_conda_envs_dirs <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (is.null(conda)) {
    return(NULL)
  }

  info <- conda_info_json(conda = conda)
  dirs <- character()
  if (!is.null(info)) {
    dirs <- c(
      dirs,
      info[["envs directories"]],
      info[["envs_dirs"]]
    )
  }

  root <- conda_root_prefix(conda = conda, info = info)
  if (!is.null(root)) {
    dirs <- c(dirs, file.path(root, "envs"))
  }

  dirs <- normalize_conda_paths(dirs)
  if (length(dirs) == 0) {
    return(NULL)
  }

  existing_dirs <- dirs[dir.exists(dirs)]
  if (length(existing_dirs) > 0) {
    return(existing_dirs)
  }
  dirs
}

get_conda_envs_dir <- function(conda = "auto") {
  envs_dirs <- get_conda_envs_dirs(conda = conda)
  envs_dirs[[1]] %||% NULL
}

conda_env_paths <- function(
  envname = NULL,
  conda = "auto",
  envs_dir = NULL
) {
  envname <- get_envname(envname)
  if (grepl("[/\\\\]", envname)) {
    return(normalize_conda_paths(envname))
  }

  envs_dirs <- envs_dir %||% get_conda_envs_dirs(conda = conda)
  if (is.null(envs_dirs) || length(envs_dirs) == 0) {
    return(character())
  }

  normalize_conda_paths(file.path(envs_dirs, envname))
}

conda_env_path <- function(
  envname = NULL,
  conda = "auto",
  envs_dir = NULL
) {
  paths <- conda_env_paths(
    envname = envname,
    conda = conda,
    envs_dir = envs_dir
  )
  if (length(paths) == 0) {
    return(NULL)
  }
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  paths[[1]]
}

infer_conda_envs_dir <- function(conda = "auto") {
  envs_dir <- get_conda_envs_dir(conda = conda)
  if (is.null(envs_dir)) {
    return(NULL)
  }
  normalizePath(envs_dir, winslash = "/", mustWork = FALSE)
}

find_conda <- function() {
  for (candidate in c(
    getOption("scop_conda", default = NULL),
    getOption("reticulate.conda_binary", default = NULL),
    Sys.getenv("RETICULATE_CONDA", unset = NA)
  )) {
    resolved <- resolve_conda_executable(
      candidate,
      error_if_missing = FALSE
    )
    if (!is.null(resolved)) {
      return(resolved)
    }
  }

  conda <- tryCatch(
    reticulate::conda_binary(conda = "auto"),
    error = identity
  )
  conda_exist <- !inherits(conda, "error")
  if (isTRUE(conda_exist)) {
    return(normalizePath(conda, winslash = "/", mustWork = FALSE))
  }

  if (isFALSE(conda_exist)) {
    if (!is.na(Sys.getenv("USER", unset = NA))) {
      miniconda_path <- gsub(
        pattern = "\\$USER",
        replacement = Sys.getenv("USER"),
        reticulate::miniconda_path()
      )
    } else {
      miniconda_path <- reticulate::miniconda_path()
    }
    conda_exist <- get_namespace_fun(
      "reticulate",
      "miniconda_exists"
    )(
      miniconda_path
    ) &&
      get_namespace_fun(
        "reticulate",
        "miniconda_test"
      )(miniconda_path)
    if (isTRUE(conda_exist)) {
      conda <- get_namespace_fun(
        "reticulate",
        "miniconda_conda"
      )(miniconda_path)
    } else {
      conda <- NULL
    }
  }

  if (!is.null(conda)) {
    return(normalizePath(conda, winslash = "/", mustWork = FALSE))
  }

  for (candidate in c("micromamba", "mamba", "conda")) {
    resolved <- resolve_conda_executable(
      candidate,
      error_if_missing = FALSE
    )
    if (!is.null(resolved)) {
      return(resolved)
    }
  }

  return(conda)
}

resolve_conda_executable <- function(
  conda,
  error_if_missing = TRUE,
  install_if_missing = FALSE
) {
  if (is.null(conda) || length(conda) == 0 || is.na(conda[[1]])) {
    return(NULL)
  }

  input <- path.expand(as.character(conda[[1]]))
  if (!nzchar(input)) {
    return(NULL)
  }

  candidate <- input
  if (!grepl("[/\\\\]", candidate)) {
    found <- Sys.which(candidate)
    if (nzchar(found)) {
      candidate <- found
    } else if (
      identical(candidate, "micromamba") &&
        isTRUE(install_if_missing)
    ) {
      candidate <- install_micromamba()
    } else if (identical(candidate, "micromamba")) {
      managed <- find_managed_micromamba()
      if (!is.null(managed)) {
        candidate <- managed
        configure_managed_micromamba_root()
      }
    }
  }

  resolved <- tryCatch(
    reticulate::conda_binary(candidate),
    error = identity
  )
  if (inherits(resolved, "error")) {
    if (isTRUE(error_if_missing)) {
      log_message(
        "Unable to find conda-compatible executable {.val {input}}: {.val {resolved$message}}",
        message_type = "error"
      )
    }
    return(NULL)
  }

  resolved <- normalizePath(resolved, winslash = "/", mustWork = FALSE)
  if (identical(conda_manager_type(resolved), "micromamba")) {
    managed <- find_managed_micromamba()
    if (!is.null(managed) && identical(resolved, managed)) {
      configure_managed_micromamba_root()
    }
  }

  resolved
}

resolve_conda <- function(conda = "auto") {
  if (identical(conda, "auto")) {
    conda <- find_conda()
  } else {
    conda_name <- as.character(conda[[1]])
    conda <- resolve_conda_executable(
      conda,
      error_if_missing = TRUE,
      install_if_missing = identical(conda_name, "micromamba")
    )
  }
  return(conda)
}

ensure_conda <- function(conda, error_if_missing = TRUE) {
  if (is.null(conda)) {
    if (error_if_missing) {
      log_message(
        "Conda-compatible environment manager not found",
        message_type = "error"
      )
    }
    return(FALSE)
  }
  return(TRUE)
}

get_conda_version <- function(conda) {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = FALSE)) {
    return(NULL)
  }

  system2t <- get_namespace_fun("reticulate", "system2t")
  version_output <- tryCatch(
    {
      system2t(conda, "--version", stdout = TRUE, stderr = TRUE)
    },
    error = function(e) NULL
  )

  if (is.null(version_output) || length(version_output) == 0) {
    return(NULL)
  }

  version_line <- paste(version_output, collapse = " ")
  version_match <- regexec(
    "\\b(?:conda|mamba|micromamba)\\s+([0-9]+(?:\\.[0-9]+){1,})\\b",
    version_line,
    perl = TRUE
  )
  match_parts <- regmatches(version_line, version_match)[[1]]
  if (length(match_parts) < 2) {
    return(NULL)
  }

  tryCatch(
    {
      utils::packageVersion(match_parts[2])
    },
    error = function(e) NULL
  )
}

conda_supports_tos <- function(conda) {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = FALSE)) {
    return(FALSE)
  }
  if (!identical(conda_manager_type(conda), "conda")) {
    return(FALSE)
  }

  system2t <- get_namespace_fun("reticulate", "system2t")
  status <- tryCatch(
    {
      system2t(
        conda,
        c("tos", "--help"),
        stdout = FALSE,
        stderr = FALSE
      )
    },
    error = function(e) NULL
  )

  identical(status, 0L)
}

install_uv <- function(
  python = NULL,
  envname = NULL,
  conda = "auto",
  pip_options = character()
) {
  log_message("Attempting to install uv...")
  pip_options <- normalize_cli_args(pip_options)

  python <- resolve_python_executable(
    python = python,
    envname = envname,
    conda = conda,
    error_message = "Python not found, cannot install uv"
  )
  if (is.null(python)) {
    return(NULL)
  }

  system2t <- get_namespace_fun("reticulate", "system2t")

  log_message("Installing {.pkg uv} using {.pkg pip}...")

  install_cmd <- c("-m", "pip", "install", "--upgrade")
  if (is.null(envname) && !is.null(python)) {
    if (!grepl("conda|envs", python)) {
      install_cmd <- c(install_cmd, "--user")
    }
  }
  if (length(pip_options) > 0) {
    install_cmd <- c(install_cmd, pip_options)
  }
  install_cmd <- c(install_cmd, "uv")
  install_status <- system2t(python, shQuote(install_cmd))

  if (install_status == 0L) {
    log_message(
      "{.pkg uv} installed successfully via {.pkg pip}",
      message_type = "success"
    )

    python_dir <- dirname(python)
    uv_path <- file.path(python_dir, "uv")
    if (is_windows()) {
      uv_path <- paste0(uv_path, ".exe")
    }

    if (file.exists(uv_path)) {
      return(uv_path)
    }

    uv_system <- Sys.which("uv")
    if (uv_system != "" && file.exists(uv_system)) {
      return(uv_system)
    }

    test_cmd <- c("-m", "uv", "--version")
    test_status <- system2t(python, shQuote(test_cmd))
    if (test_status == 0L) {
      return("python -m uv")
    }

    log_message(
      "{.pkg uv} installed but not found in expected locations",
      message_type = "warning"
    )
    return(NULL)
  } else {
    log_message(
      "Failed to install {.pkg uv} [error code {.val {install_status}}]",
      message_type = "warning"
    )
    return(NULL)
  }
}

find_uv <- function(
  python = NULL,
  envname = NULL,
  conda = "auto",
  auto_install = TRUE,
  pip_options = character()
) {
  pip_options <- normalize_cli_args(pip_options)
  uv <- Sys.which("uv")
  if (uv != "" && file.exists(uv)) {
    return(uv)
  }

  if (is.null(python)) {
    python <- resolve_python_executable(
      envname = envname,
      conda = conda,
      error_if_missing = FALSE
    )
  }

  if (!is.null(python) && file.exists(python)) {
    python_dir <- dirname(python)
    uv_path <- file.path(python_dir, "uv")
    if (is_windows()) {
      uv_path <- paste0(uv_path, ".exe")
    }
    if (file.exists(uv_path)) {
      return(uv_path)
    }

    system2t <- get_namespace_fun("reticulate", "system2t")
    test_cmd <- c("-m", "uv", "--version")
    test_status <- system2t(python, shQuote(test_cmd))
    if (test_status == 0L) {
      return("python -m uv")
    }
  }

  if (isTRUE(auto_install)) {
    uv <- install_uv(
      python = python,
      envname = envname,
      conda = conda,
      pip_options = pip_options
    )
    if (!is.null(uv)) {
      return(uv)
    }
  }

  return(NULL)
}

create_conda_env <- function(
  conda,
  envname,
  python_version = NULL,
  packages = character(),
  forge = TRUE,
  channel = character()
) {
  manager <- conda_manager_type(conda)
  if (identical(manager, "conda")) {
    return(reticulate::conda_create(
      conda = conda,
      envname = envname,
      python_version = python_version,
      packages = packages,
      forge = forge,
      channel = channel
    ))
  }

  envname <- get_namespace_fun(
    "reticulate",
    "condaenv_resolve"
  )(envname)
  conda_args <- get_namespace_fun(
    "reticulate",
    "conda_args"
  )
  system2t <- get_namespace_fun(
    "reticulate",
    "system2t"
  )

  python_package <- if (is.null(python_version)) {
    "python"
  } else {
    sprintf("python=%s", python_version)
  }

  args <- conda_args("create", envname, c(python_package, packages))
  args <- c(args, "--quiet")

  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  } else {
    character()
  }
  for (ch in channels) {
    args <- c(args, "-c", ch)
  }

  result <- system2t(conda, shQuote(args))
  if (result != 0L) {
    log_message(
      "Error creating Python environment {.file {envname}} with {.pkg {manager}} [exit code {.val {result}}]",
      message_type = "error"
    )
  }

  conda_python(envname = envname, conda = conda)
}

conda_install <- function(
  envname = NULL,
  packages,
  forge = TRUE,
  channel = character(),
  pip = FALSE,
  pip_options = character(),
  pip_ignore_installed = FALSE,
  conda = "auto",
  python_version = NULL,
  ...
) {
  envname <- get_envname(envname)
  get_namespace_fun(
    "reticulate",
    "check_forbidden_install"
  )("Python packages")
  conda_args <- get_namespace_fun(
    "reticulate",
    "conda_args"
  )
  system2t <- get_namespace_fun(
    "reticulate",
    "system2t"
  )

  if (missing(packages)) {
    if (!is.null(envname)) {
      fmt <- paste(
        "argument \"packages\" is missing, with no default",
        "- did you mean 'conda_install(<envname>, %1$s)'?",
        "- use 'py_install(%1$s)' to install into the active Python environment",
        sep = "\n"
      )
      log_message(
        sprintf(fmt, deparse1(substitute(envname))),
        message_type = "error"
      )
    } else {
      log_message(
        "{.arg packages} argument is required",
        message_type = "error"
      )
    }
  }

  conda <- resolve_conda(conda)
  if (!ensure_conda(conda)) {
    return(invisible(packages))
  }
  manager <- conda_manager_label(conda)
  envname <- get_namespace_fun(
    "reticulate",
    "condaenv_resolve"
  )(envname)

  log_message(
    "Installing {.val {length(packages)}} packages into environment: {.file {envname}} using {.pkg {manager}}"
  )

  python_package <- if (is.null(python_version)) {
    NULL
  } else if (grepl("[><=]", python_version)) {
    paste0("python", python_version)
  } else {
    sprintf("python=%s", python_version)
  }

  python <- tryCatch(
    conda_python(envname = envname, conda = conda),
    error = identity
  )

  if (inherits(python, "error") || !file.exists(python)) {
    log_message(
      "Python environment does not exist, creating: {.file {envname}}"
    )
    accept_conda_tos(conda = conda)
    create_conda_env(
      conda = conda,
      envname = envname,
      python_version = python_version,
      packages = character(),
      forge = forge,
      channel = channel
    )
    python <- conda_python(envname = envname, conda = conda)
    log_message(
      "Environment created with python: {.pkg {python}}"
    )
  }

  if (!is.null(python_version)) {
    log_message(
      "Updating python version to: {.pkg {python_version}}"
    )
    args <- conda_args("install", envname, python_package)
    status <- system2t(conda, shQuote(args))
    if (status != 0L) {
      fmt <- "installation of '%s' into environment '%s' failed [error code %i]"
      msg <- sprintf(fmt, python_package, envname, status)
      log_message(msg, message_type = "error")
    } else {
      log_message(
        "Python version updated successfully",
        message_type = "success"
      )
    }
  }

  if (pip) {
    install_python_packages(
      packages = packages,
      python = python,
      envname = envname,
      conda = conda,
      pip_options = pip_options,
      pip_ignore_installed = pip_ignore_installed
    )

    return(invisible(packages))
  }

  log_message("Installing packages via {.pkg {manager}}...")
  args <- conda_args("install", envname)

  channels <- if (length(channel)) {
    channel
  } else if (forge) {
    "conda-forge"
  }

  for (ch in channels) {
    args <- c(args, "-c", ch)
    log_message("Using channel: {.val {ch}}")
  }

  conda_packages <- if (is.character(packages)) {
    gsub("==", "=", packages, fixed = TRUE)
  } else {
    packages
  }

  args <- c(args, python_package, conda_packages)

  log_message("Installing {.val {length(packages)}} packages...")

  result <- system2t(conda, shQuote(args))

  if (result != 0L) {
    log_message(
      "{.pkg {manager}} installation failed with error code: {.val {result}}",
      message_type = "warning"
    )
  } else {
    log_message(
      "{.pkg {manager}} installation completed successfully",
      message_type = "success"
    )
  }

  invisible(packages)
}

resolve_python_executable <- function(
  python = NULL,
  envname = NULL,
  conda = "auto",
  error_if_missing = TRUE,
  error_message = NULL
) {
  if (!is.null(python)) {
    python <- normalizePath(python, mustWork = FALSE)
    if (file.exists(python)) {
      return(python)
    }
  }

  if (!is.null(envname)) {
    python <- tryCatch(
      conda_python(envname = envname, conda = conda),
      error = function(e) {
        if (error_if_missing) {
          log_message(
            "Failed to get Python path for environment {.file {envname}}: {.val {e$message}}",
            message_type = "warning"
          )
        }
        NULL
      }
    )
    if (!is.null(python) && file.exists(python)) {
      return(normalizePath(python, mustWork = FALSE))
    }
  }

  for (candidate in c(Sys.which("python3"), Sys.which("python"))) {
    if (nzchar(candidate) && file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }

  if (error_if_missing) {
    log_message(
      error_message %||% "Python executable not found",
      message_type = "error"
    )
  }

  NULL
}

run_uv_command <- function(uv, python, args) {
  system2t <- get_namespace_fun("reticulate", "system2t")

  if (uv == "python -m uv") {
    return(system2t(python, shQuote(c("-m", "uv", args))))
  }

  system2t(uv, shQuote(args))
}

run_pip_command <- function(python, args) {
  system2t <- get_namespace_fun("reticulate", "system2t")
  system2t(python, shQuote(c("-m", "pip", args)))
}

install_python_packages <- function(
  packages,
  python,
  envname,
  conda,
  pip_options = character(),
  pip_ignore_installed = FALSE
) {
  pip_options <- normalize_cli_args(pip_options)
  uv <- find_uv(
    python = python,
    envname = envname,
    conda = conda,
    auto_install = TRUE,
    pip_options = pip_options
  )

  if (!is.null(uv)) {
    log_message("Installing packages via {.pkg uv}...")

    args <- c("pip", "install")
    if (uv != "python -m uv") {
      args <- c(args, "--python", python)
    }
    if (length(pip_options) > 0) {
      args <- c(args, pip_options)
    }
    if (pip_ignore_installed) {
      args <- c(args, "--reinstall")
    }
    args <- c(args, packages)

    status <- run_uv_command(uv, python, args)
    if (status == 0L) {
      log_message(
        "{.pkg uv} installation completed",
        message_type = "success"
      )
      return(invisible(TRUE))
    }

    log_message(
      "{.pkg uv} installation failed [error code {.val {status}}], using {.pkg pip}",
      message_type = "warning"
    )
  } else {
    log_message(
      "{.pkg uv} not available, using {.pkg pip}...",
      message_type = "warning"
    )
  }

  get_namespace_fun("reticulate", "pip_install")(
    python = python,
    packages = packages,
    pip_options = pip_options,
    ignore_installed = pip_ignore_installed,
    conda = conda,
    envname = envname
  )

  invisible(TRUE)
}

ensure_windows_scvi_support <- function(
  envname,
  conda,
  force = FALSE,
  pip_options = character()
) {
  if (!is_windows()) {
    return(invisible(FALSE))
  }

  jax_requirement <- c("jax" = "jax[cpu]==0.4.38")

  if (isTRUE(exist_python_pkgs(
    jax_requirement,
    envname = envname,
    conda = conda
  ))) {
    return(invisible(TRUE))
  }
  pip_options <- normalize_cli_args(pip_options)
  status <- isTRUE(check_python(
    packages = jax_requirement,
    envname = envname,
    conda = conda,
    force = FALSE,
    pip = TRUE,
    pip_options = pip_options,
    verbose = FALSE
  ))

  if (isTRUE(status)) {
    log_message(
      "{.pkg jax[cpu]} verified for Windows {.pkg scvi-tools} support",
      message_type = "success"
    )
  } else {
    log_message(
      "Failed to verify {.pkg jax[cpu]} for Windows {.pkg scvi-tools} support",
      message_type = "warning"
    )
  }

  invisible(status)
}

conda_python <- function(
  envname = NULL,
  conda = "auto",
  all = FALSE
) {
  envname <- get_envname(envname)
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda)) {
    return(NULL)
  }
  envname <- get_namespace_fun(
    "reticulate",
    "python_environment_resolve"
  )(envname)
  if (grepl("[/\\\\]", envname)) {
    path <- conda_env_python_path(envname)
    if (file.exists(path)) {
      return(path)
    }
    log_message(
      "No conda-compatible Python environment exists at path {.file {envname}}",
      message_type = "error"
    )
  }
  envs_dirs <- get_conda_envs_dirs(conda = conda)
  for (env_path in conda_env_paths(
    envname = envname,
    conda = conda,
    envs_dir = envs_dirs
  )) {
    python_path <- conda_env_python_path(env_path)

    if (file.exists(python_path) && valid_conda_env_paths(env_path)) {
      return(normalizePath(python_path, mustWork = FALSE))
    }
  }

  conda_envs <- tryCatch(
    reticulate::conda_list(conda = conda),
    error = function(...) data.frame(name = character(), python = character())
  )
  conda_envs <- filter_conda_env_table(
    conda_envs = conda_envs,
    envs_dirs = envs_dirs,
    root_prefix = conda_root_prefix(conda = conda)
  )
  env <- conda_envs[conda_envs$name == envname, , drop = FALSE]
  if (nrow(env) == 0) {
    for (env_path in conda_env_paths(
      envname = envname,
      conda = conda,
      envs_dir = envs_dirs
    )) {
      python_path <- conda_env_python_path(env_path)

      if (file.exists(python_path)) {
        return(normalizePath(python_path, mustWork = FALSE))
      }
    }

    log_message(
      "{.val {envname}} environment not found",
      message_type = "error"
    )
  }
  python <- if (all) env$python else env$python[[1L]]
  return(normalizePath(as.character(python), mustWork = FALSE))
}

#' @title Remove a conda-compatible Python environment
#'
#' @md
#' @inheritParams PrepareEnv
#' @param force Whether to force removal without confirmation.
#' Default is `FALSE`.
#'
#' @return Invisibly returns `TRUE` if successful, `FALSE` otherwise.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove default environment
#' RemoveEnv()
#'
#' # Remove a specific environment
#' RemoveEnv("my_old_env")
#'
#' # Removal without confirmation
#' RemoveEnv("my_old_env", force = TRUE)
#' }
RemoveEnv <- function(
  envname = NULL,
  conda = "auto",
  force = FALSE,
  verbose = TRUE
) {
  envname <- get_envname(envname)

  conda <- resolve_conda(conda)
  if (!ensure_conda(conda)) {
    return(invisible(FALSE))
  }
  manager <- conda_manager_label(conda)

  log_message(
    "Removing environment: {.file {envname}} using {.pkg {manager}}",
    verbose = verbose
  )

  env_exists <- env_exist(envname = envname, conda = conda)
  if (isFALSE(env_exists)) {
    log_message(
      "{.file {envname}} environment does not exist",
      message_type = "warning",
      verbose = verbose
    )
    return(invisible(FALSE))
  }

  env_path <- conda_env_path(envname = envname, conda = conda)

  if (!force) {
    log_message(
      "Environment path: {.file {env_path}}",
      verbose = verbose
    )
    log_message(
      "This will permanently delete the environment and all its packages",
      verbose = verbose
    )

    if (interactive()) {
      response <- readline(
        paste0(
          "Are you sure you want to remove environment: ",
          envname,
          "? (y/N): "
        )
      )
      if (!tolower(response) %in% c("y", "yes")) {
        log_message(
          "Environment removal cancelled",
          message_type = "warning",
          verbose = verbose
        )
        return(invisible(FALSE))
      }
    } else {
      log_message(
        "Use {.arg force = TRUE} to remove environment in non-interactive mode",
        message_type = "warning",
        verbose = verbose
      )
      return(invisible(FALSE))
    }
  }

  log_message("Removing {.file {envname}} environment...", verbose = verbose)

  result <- tryCatch(
    {
      reticulate::conda_remove(
        envname = envname,
        packages = NULL,
        conda = conda
      )
      log_message(
        "Environment removed successfully using {.fn reticulate::conda_remove}",
        message_type = "success",
        verbose = verbose
      )
      TRUE
    },
    error = function(e) {
      log_message(
        "{.fn reticulate::conda_remove} failed for {.pkg {manager}}: {.val {e$message}}",
        message_type = "warning",
        verbose = verbose
      )
      FALSE
    }
  )

  if (!result) {
    log_message(
      "Attempting direct removal of environment directory: {.file {env_path}}",
      verbose = verbose
    )

    result <- tryCatch(
      {
        if (dir.exists(env_path)) {
          unlink(env_path, recursive = TRUE, force = TRUE)

          if (!dir.exists(env_path)) {
            log_message(
              "Environment directory {.file {env_path}} removed successfully",
              message_type = "success",
              verbose = verbose
            )
            TRUE
          } else {
            log_message(
              "Failed to remove environment directory: {.file {env_path}}",
              message_type = "error"
            )
            FALSE
          }
        } else {
          log_message(
            "Environment directory does not exist: {.file {env_path}}",
            message_type = "warning",
            verbose = verbose
          )
          TRUE
        }
      },
      error = function(e) {
        log_message(
          "Direct removal failed: {.val {e$message}}",
          message_type = "error"
        )
        FALSE
      }
    )
  } else if (!is.null(env_path) && dir.exists(env_path)) {
    log_message(
      "Removing leftover environment directory: {.file {env_path}}",
      verbose = verbose
    )
    result <- tryCatch(
      {
        unlink(env_path, recursive = TRUE, force = TRUE)
        !dir.exists(env_path)
      },
      error = function(e) {
        log_message(
          "Direct cleanup failed: {.val {e$message}}",
          message_type = "warning",
          verbose = verbose
        )
        FALSE
      }
    )
  }

  if (result) {
    clear_env_cache()
    Sys.unsetenv("RETICULATE_PYTHON")
    log_message(
      "{.file {envname}} environment removed successfully",
      message_type = "success",
      verbose = verbose
    )
  } else {
    log_message(
      "Failed to remove {.file {envname}} environment",
      message_type = "error"
    )
  }

  return(invisible(result))
}

#' @title List conda-compatible Python environments
#'
#' @md
#' @inheritParams PrepareEnv
#' @return A data frame of conda-compatible Python environments.
#' @export
ListEnv <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = TRUE)) {
    return(invisible(NULL))
  }

  envs_dirs <- get_conda_envs_dirs(conda = conda)
  conda_envs <- tryCatch(
    reticulate::conda_list(conda = conda),
    error = function(...) data.frame(name = character(), python = character())
  )
  conda_envs <- filter_conda_env_table(
    conda_envs = conda_envs,
    envs_dirs = envs_dirs,
    root_prefix = conda_root_prefix(conda = conda)
  )
  for (envs_dir in envs_dirs) {
    if (is.null(envs_dir) || !dir.exists(envs_dir)) {
      next
    }
    existing_names <- conda_envs$name

    env_dirs <- list.dirs(envs_dir, full.names = FALSE, recursive = FALSE)

    for (env_dir in env_dirs) {
      if (!env_dir %in% existing_names) {
        env_path <- file.path(envs_dir, env_dir)
        python_path <- file.path(
          env_path,
          conda_env_python_suffix()
        )

        if (
          valid_conda_env_paths(env_path) &&
            file.exists(python_path)
        ) {
          new_row <- data.frame(
            name = env_dir,
            python = normalizePath(python_path, mustWork = FALSE),
            stringsAsFactors = FALSE
          )
          conda_envs <- rbind(conda_envs, new_row)
        }
      }
    }
  }

  if (nrow(conda_envs) > 0) {
    conda_envs <- conda_envs[order(conda_envs$name), , drop = FALSE]
  }

  return(conda_envs)
}

accept_conda_tos <- function(conda = "auto") {
  conda <- resolve_conda(conda)
  if (!ensure_conda(conda, error_if_missing = FALSE)) {
    return(invisible(FALSE))
  }

  manager <- conda_manager_type(conda)
  if (!identical(manager, "conda")) {
    log_message(
      "Skipping conda Terms of Service acceptance for {.pkg {manager}}."
    )
    return(invisible(FALSE))
  }

  system2t <- get_namespace_fun("reticulate", "system2t")
  conda_version <- get_conda_version(conda)
  supports_tos <- conda_supports_tos(conda)

  if (!supports_tos) {
    version_text <- if (is.null(conda_version)) {
      "unknown"
    } else {
      as.character(conda_version)
    }
    log_message(
      "Skipping conda ToS acceptance because {.pkg conda} does not support the {.val tos} command (version: {.val {version_text}}).",
      message_type = "warning"
    )
    return(invisible(FALSE))
  }

  channels <- c(
    "https://repo.anaconda.com/pkgs/main",
    "https://repo.anaconda.com/pkgs/r",
    "https://repo.anaconda.com/pkgs/msys2"
  )

  log_message(
    "Accepting conda Terms of Service for required channels..."
  )

  for (channel in channels) {
    tryCatch(
      {
        args <- c("tos", "accept", "--override-channels", "--channel", channel)
        status <- system2t(conda, args, stdout = FALSE, stderr = FALSE)
        if (status == 0L) {
          log_message(
            "Accepted ToS for channel: {.val {channel}}",
            message_type = "success"
          )
        } else {
          log_message(
            "Failed to accept ToS for channel: {.val {channel}} [error code {.val {status}}]",
            message_type = "warning"
          )
        }
      },
      error = function(e) {
        log_message(
          "Failed to accept ToS for channel: {.val {channel}} ({.val {e$message}})",
          message_type = "warning"
        )
      }
    )
  }

  return(invisible(TRUE))
}
