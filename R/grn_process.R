grn_python <- function(config, module) {
  envname <- config$envname %||% env_module_environment(module)
  python <- config$runtime$python %||% tryCatch(
    conda_python(envname = envname, conda = config$conda %||% "auto"),
    error = function(e) NULL
  )
  if (!is.character(python) || length(python) != 1L || is.na(python) ||
    !file.exists(python) || dir.exists(python)) {
    stop("Python runtime was not found. Run PrepareEnv(modules = \"", module,
      "\") or specify method_args$runtime$python.",
      call. = FALSE
    )
  }
  file.path(normalizePath(dirname(python), winslash = "/", mustWork = TRUE), basename(python))
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
  result <- list(
    command = command, args = args, stdout = stdout, stderr = stderr,
    started_at = started, status = "running", exit_status = NA_integer_
  )
  process <- processx::process$new(command, args,
    stdout = stdout, stderr = stderr,
    env = if (length(env)) c("current", env) else NULL,
    wd = directory, cleanup = TRUE, cleanup_tree = FALSE, windows_hide_window = TRUE
  )
  on.exit(
    {
      process$kill_tree()
      result$finished_at <- Sys.time()
      saveRDS(result, receipt)
    },
    add = TRUE
  )
  tryCatch(
    {
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
    },
    interrupt = function(e) {
      result$status <<- "cancelled"
      process$kill_tree()
      stop(structure(list(message = paste0("GRN subprocess cancelled. Logs: ", stdout, "; ", stderr), call = NULL), class = c("interrupt", "condition")))
    },
    error = function(e) {
      if (identical(result$status, "running")) result$status <<- "failed"
      stop(e)
    }
  )
  result$finished_at <- Sys.time()
  result$receipt <- receipt
  result
}
