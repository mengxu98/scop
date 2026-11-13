grn_runtime_profiles <- function() {
  list(
    native = list(kind = "r", environment = NULL, packages = character()),
    r = list(kind = "r", environment = NULL, packages = character()),
    genie3 = list(kind = "r", environment = NULL, packages = "GENIE3"),
    ppcor = list(kind = "r", environment = NULL, packages = "ppcor"),
    pidc = list(kind = "r", environment = NULL, packages = c("pbapply", "purrr", "thisutils")),
    scorpion = list(kind = "r", environment = NULL, packages = "SCORPION"),
    leap = list(kind = "r", environment = NULL, packages = "LEAP"),
    scode = list(kind = "r", environment = NULL, packages = "MASS"),
    sincerities = list(kind = "r", environment = NULL, packages = c("glmnet", "ppcor", "thisutils")),
    multicsn = list(kind = "r", environment = NULL, packages = c("multiCSN", "inferCSN", "Seurat", "SeuratObject", "Signac", "Matrix", "GenomicRanges", "GenomeInfoDb", "motifmatchr")),
    multicsn_rna = list(kind = "r", environment = NULL, packages = c("multiCSN", "inferCSN", "Matrix")),
    pando = list(kind = "r", environment = NULL, packages = c("Pando", "Signac", "Seurat", "SeuratObject", "GenomicRanges", "motifmatchr", "glmnet", "dplyr")),
    grnboost2 = list(kind = "python", environment = "grnboost2_env", packages = c("arboreto", "numpy", "pandas", "dask", "distributed")),
    scenic = list(kind = "python", environment = "scenic_env", packages = c("pyscenic", "arboreto", "ctxcore")),
    gniplr = list(kind = "python", environment = "gniplr_env", packages = c("numpy", "pandas", "scipy", "scikit-learn", "statsmodels")),
    regdiffusion = list(kind = "python", environment = "regdiffusion_env", packages = c("regdiffusion", "torch", "numpy", "pandas")),
    scsgl = list(kind = "python", environment = "scsgl_env", packages = c("numpy", "pandas", "scipy", "numba")),
    linger = list(kind = "python", environment = "linger_env", packages = c("LingerGRN", "torch", "numpy", "pandas")),
    scenicplus = list(kind = "python", environment = "scenicplus_env", packages = c("scenicplus", "pycisTopic", "pycistarget", "snakemake", "pybiomart", "gseapy"))
  )
}

grn_pando_runtime_packages <- function(config) {
  if (isTRUE(config$cores > 1L)) c("foreach", "doParallel") else character()
}

grn_resolve_runtime <- function(config, profile) {
  spec <- if (is.list(profile)) profile else grn_runtime_profiles()[[profile]]
  if (is.null(spec) || !spec$kind %in% c("r", "python")) {
    stop("Unknown GRN runtime profile: ", paste(profile, collapse = ", "), call. = FALSE)
  }
  genome <- config$resources$genome
  if (identical(config$method, "pando") && is.null(genome)) genome <- "BSgenome.Hsapiens.UCSC.hg38"
  if (identical(spec$kind, "r") && is.character(genome) && length(genome) == 1L && startsWith(genome, "BSgenome.")) {
    if (!grepl("^BSgenome[.][A-Za-z0-9.]+$", genome)) stop("Invalid BSgenome package name.", call. = FALSE)
    spec$packages <- unique(c(spec$packages, genome))
  }
  if (identical(config$method, "pando") && identical(spec$kind, "r")) {
    spec$packages <- unique(c(spec$packages, grn_pando_runtime_packages(config)))
  }
  runtime <- config$runtime %||% list()
  if (inherits(runtime, "GRNRuntime")) {
    if (identical(attr(runtime, "profile"), grn_hash(spec))) {
      if (!file.exists(runtime$executable)) stop("The resolved GRN executable no longer exists: ", runtime$executable, call. = FALSE)
      return(runtime)
    }
    runtime <- runtime[intersect(names(runtime), c("python", "rscript", "r_library", "timeout"))]
  }
  invalid <- setdiff(names(runtime), c("python", "rscript", "r_library", "timeout"))
  if (length(invalid)) stop("Unknown GRN runtime option(s): ", paste(invalid, collapse = ", "), call. = FALSE)
  timeout <- runtime$timeout %||% Inf
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) {
    stop("GRN runtime timeout must be a positive number of seconds.", call. = FALSE)
  }
  executable <- if (identical(spec$kind, "python")) {
    runtime$python %||% tryCatch(
      conda_python(envname = config$envname %||% spec$environment, conda = config$conda %||% "auto"),
      error = function(e) NULL
    )
  } else {
    runtime$rscript %||% file.path(R.home("bin"), "Rscript")
  }
  if (is.null(executable) || length(executable) != 1L || !nzchar(executable) || !file.exists(executable)) {
    stop("The GRN runtime was not found. Supply method_args$runtime$", if (spec$kind == "python") "python" else "rscript", " or prepare the ", spec$environment %||% "R", " environment with PrepareEnv().", call. = FALSE)
  }
  executable <- file.path(normalizePath(dirname(executable), winslash = "/", mustWork = TRUE), basename(executable))
  r_library <- runtime$r_library %||% .libPaths()
  if (!is.character(r_library) || anyNA(r_library) || any(!dir.exists(r_library))) {
    stop("Every runtime r_library entry must name an existing library directory.", call. = FALSE)
  }
  r_library <- normalizePath(r_library, winslash = "/", mustWork = TRUE)
  packages <- spec$packages %||% character()
  package_names <- sub("[<>=!~].*$", "", packages)
  exact_versions <- ifelse(grepl("==", packages, fixed = TRUE), sub("^.*==", "", packages), "")
  if (any(!grepl("^[A-Za-z0-9_.-]+$", package_names)) || any(!grepl("^[A-Za-z0-9.+_-]*$", exact_versions)) ||
      any(grepl("[<>=!~]", packages) & !grepl("^[A-Za-z0-9_.-]+==[A-Za-z0-9.+_-]+$", packages))) {
    stop("GRN runtime package requirements must be names or exact name==version requirements.", call. = FALSE)
  }
  probe_config <- config
  probe_config$runtime <- list(timeout = min(timeout, 60))
  probe_config$work_dir <- file.path(config$work_dir %||% tempdir(), "runtime-check")
  if (spec$kind == "python") {
    quoted <- paste(sprintf("('%s','%s')", package_names, exact_versions), collapse = ",")
    code <- paste0("import importlib.metadata as m,sys\nprint('RUNTIME\\t'+sys.version.split()[0])\nmissing=[]\nfor p,expected in [", quoted, "]:\n try:\n  v=m.version(p)\n  print('PACKAGE\\t'+p+'\\t'+v)\n  if expected and v != expected: missing.append(p+'=='+expected+' (found '+v+')')\n except m.PackageNotFoundError: missing.append(p)\nif missing: raise RuntimeError('Missing or incompatible GRN packages: '+', '.join(missing))\n")
    observed <- grn_run_process(executable, c("-c", code), probe_config)
  } else {
    code <- paste0(".libPaths(c(", paste(encodeString(r_library, quote = '"'), collapse = ","), "), include.site=FALSE);cat('RUNTIME\t', as.character(getRversion()), '\n', sep=''); pkgs <- c(", paste(encodeString(package_names, quote = '"'), collapse = ","), "); missing <- pkgs[!vapply(pkgs,function(p)length(find.package(p,quiet=TRUE))==1L,logical(1))]; if(length(missing)) stop('Missing GRN packages: ',paste(missing,collapse=', '));for(p in pkgs) cat('PACKAGE\t', p, '\t',as.character(utils::packageVersion(p)), '\n',sep='')")
    observed <- grn_run_process(executable, c("--vanilla", "-e", code), probe_config)
  }
  lines <- readLines(observed$stdout, warn = FALSE)
  version <- sub("^RUNTIME\t", "", lines[startsWith(lines, "RUNTIME\t")])
  records <- strsplit(lines[startsWith(lines, "PACKAGE\t")], "\t", fixed = TRUE)
  versions <- if (length(records)) stats::setNames(vapply(records, `[[`, "", 3L), vapply(records, `[[`, "", 2L)) else character()
  mismatch <- package_names[nzchar(exact_versions) & versions[package_names] != exact_versions]
  if (length(mismatch)) stop("Incompatible GRN package version(s): ", paste(mismatch, collapse = ", "), call. = FALSE)
  resolved <- list(kind = spec$kind, executable = executable, python = if (spec$kind == "python") executable else NULL,
       rscript = if (spec$kind == "r") executable else NULL, r_library = r_library, timeout = timeout,
       version = version, package_versions = versions, environment = spec$environment)
  structure(resolved, class = c("GRNRuntime", "list"), profile = grn_hash(spec))
}

grn_run_process <- function(command, args, config, env = character()) {
  if (!requireNamespace("processx", quietly = TRUE)) stop("Install processx before running GRN subprocesses.", call. = FALSE)
  if (!is.character(command) || length(command) != 1L || !nzchar(command)) stop("GRN command must be one executable path.", call. = FALSE)
  if (!is.character(args) || anyNA(args)) stop("GRN process arguments must be character values.", call. = FALSE)
  if (length(env) && (is.null(names(env)) || any(!nzchar(names(env))))) stop("GRN process environment variables must be named.", call. = FALSE)
  executable_directory <- if (grepl("[/\\\\]", command)) dirname(command) else dirname(Sys.which(command))
  if (length(executable_directory) == 1L && nzchar(executable_directory) && dir.exists(executable_directory)) {
    path <- if ("PATH" %in% names(env)) env[["PATH"]] else Sys.getenv("PATH")
    paths <- unique(c(normalizePath(executable_directory, winslash = "/", mustWork = TRUE), strsplit(path, .Platform$path.sep, fixed = TRUE)[[1L]]))
    env <- c(env[names(env) != "PATH"], PATH = paste(paths[nzchar(paths)], collapse = .Platform$path.sep))
  }
  directory <- config$work_dir %||% tempdir()
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  token <- basename(tempfile("process-"))
  stdout <- file.path(directory, paste0(token, ".stdout.log"))
  stderr <- file.path(directory, paste0(token, ".stderr.log"))
  receipt <- file.path(directory, paste0(token, ".rds"))
  started <- Sys.time()
  timeout <- config$runtime$timeout %||% Inf
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout <= 0) stop("GRN process timeout must be positive seconds.", call. = FALSE)
  result <- list(command = command, args = args, stdout = stdout, stderr = stderr,
                 started_at = started, status = "running", exit_status = NA_integer_)
  process <- processx::process$new(command, args, stdout = stdout, stderr = stderr,
                                   env = if (length(env)) c("current", env) else NULL,
                                   wd = directory, cleanup = TRUE, cleanup_tree = FALSE, windows_hide_window = TRUE)
  on.exit({
    process$kill_tree()
    result$finished_at <- Sys.time()
    saveRDS(result, receipt)
  }, add = TRUE)
  tryCatch({
    repeat {
      if (!process$is_alive()) break
      if (as.numeric(difftime(Sys.time(), started, units = "secs")) > timeout) {
        result$status <- "timed_out"
        process$kill_tree()
        stop("GRN subprocess timed out. Logs: ", stdout, "; ", stderr, call. = FALSE)
      }
      process$wait(timeout = 100)
    }
    result$exit_status <- process$get_exit_status()
    result$status <- if (identical(result$exit_status, 0L)) "completed" else "failed"
    if (!identical(result$exit_status, 0L)) {
      detail <- paste(utils::tail(readLines(stderr, warn = FALSE), 10L), collapse = "\n")
      stop("GRN subprocess failed with exit status ", result$exit_status, ". Logs: ", stdout, "; ", stderr, if (nzchar(detail)) paste0("\n", detail), call. = FALSE)
    }
  }, interrupt = function(e) {
    result$status <<- "cancelled"
    process$kill_tree()
    stop(structure(list(message = paste0("GRN subprocess cancelled. Logs: ", stdout, "; ", stderr), call = NULL), class = c("interrupt", "condition")))
  }, error = function(e) {
    if (identical(result$status, "running")) result$status <<- "failed"
    stop(e)
  })
  result$finished_at <- Sys.time()
  result$receipt <- receipt
  result
}
