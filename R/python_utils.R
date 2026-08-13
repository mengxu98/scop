scop_python_import <- function(module, convert = TRUE) {
  python_dir <- system.file("python", package = "thisutils", mustWork = FALSE)
  if (!nzchar(python_dir)) {
    log_message(
      "thisutils ({.code >= 0.4.8}) does not provide a Python {.file log_message} module",
      message_type = "error"
    )
  }
  if (isFALSE(reticulate::py_available(initialize = FALSE))) {
    configured_python <- normalize_python_runtime_path(
      Sys.getenv("RETICULATE_PYTHON", unset = "")
    )
    if (is.null(configured_python) || !file.exists(configured_python)) {
      log_message(
        "Python is not configured. Run {.code PrepareEnv()} before importing scop Python helpers",
        message_type = "error"
      )
    }
    initialization_error <- ""
    initialized <- tryCatch(
      {
        reticulate::use_python(configured_python, required = TRUE)
        isTRUE(reticulate::py_available(initialize = TRUE))
      },
      error = function(e) {
        initialization_error <<- conditionMessage(e)
        FALSE
      }
    )
    if (!initialized) {
      detail <- if (nzchar(initialization_error)) {
        paste0(": ", initialization_error)
      } else {
        ""
      }
      log_message(
        "Unable to initialize configured Python runtime {.file {configured_python}}{detail}",
        message_type = "error"
      )
    }
    assert_python_runtime_switchable(configured_python)
  }
  if (isFALSE(reticulate::py_eval(
    "'log_message' in __import__('sys').modules",
    convert = TRUE
  ))) {
    reticulate::import_from_path(
      "log_message",
      path = python_dir,
      convert = FALSE
    )
  }
  reticulate::import_from_path(
    module,
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = convert
  )
}

runner_script_path <- function(script, backend) {
  candidates <- c(
    system.file("python", script, package = "scop", mustWork = FALSE),
    file.path("inst", "python", script)
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    log_message(
      "Bundled {.pkg {backend}} Python runner was not found",
      message_type = "error"
    )
  }
  normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}

runner_write_json <- function(x, path) {
  check_r("jsonlite", verbose = FALSE)
  to_json <- get_namespace_fun("jsonlite", "toJSON")
  content <- as.character(to_json(
    x,
    auto_unbox = TRUE,
    null = "null",
    digits = NA,
    pretty = TRUE
  ))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(content, con = temporary, useBytes = TRUE)
  renamed <- suppressWarnings(file.rename(temporary, path))
  if (!isTRUE(renamed)) {
    log_message(
      "Unable to atomically write JSON file {.file {path}}",
      message_type = "error"
    )
  }
  invisible(path)
}

runner_read_json <- function(path) {
  check_r("jsonlite", verbose = FALSE)
  get_namespace_fun("jsonlite", "fromJSON")(path, simplifyVector = FALSE)
}

runner_system2 <- function(command, args, env, stdout, stderr) {
  if (length(env) && (is.null(names(env)) || any(!nzchar(names(env))))) {
    log_message(
      "Subprocess environment variables must be named",
      message_type = "error"
    )
  }
  log_message_dir <- system.file("python", package = "thisutils", mustWork = FALSE)
  if (nzchar(log_message_dir)) {
    python_path <- unique(c(
      log_message_dir,
      strsplit(Sys.getenv("PYTHONPATH", unset = ""), .Platform$path.sep, fixed = TRUE)[[1]]
    ))
    python_path <- python_path[nzchar(python_path)]
    env <- c(env, PYTHONPATH = paste(python_path, collapse = .Platform$path.sep))
  }
  if (!length(env)) {
    return(system2(
      command = command,
      args = args,
      stdout = stdout,
      stderr = stderr
    ))
  }
  old_env <- Sys.getenv(names(env), unset = NA_character_)
  names(old_env) <- names(env)
  if (length(env)) {
    do.call(Sys.setenv, as.list(env))
  }
  on.exit(
    {
      restore <- !is.na(old_env)
      if (any(restore)) {
        do.call(Sys.setenv, as.list(old_env[restore]))
      }
      if (any(!restore)) {
        Sys.unsetenv(names(old_env)[!restore])
      }
    },
    add = TRUE
  )
  system2(
    command = command,
    args = args,
    stdout = stdout,
    stderr = stderr
  )
}

runner_request_id <- function() {
  paste0(
    Sys.getpid(),
    "-",
    basename(tempfile("request_"))
  )
}

runner_acquire_lock <- function(path, backend) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  token <- runner_request_id()
  connection <- suppressWarnings(tryCatch(
    file(path, open = "wx", encoding = "UTF-8"),
    error = function(...) NULL
  ))
  if (is.null(connection)) {
    log_message(
      "Another {.pkg {backend}} run is using {.arg result_dir}. If no run is active, remove the stale lock file {.file {path}}.",
      message_type = "error"
    )
  }
  wrote <- tryCatch(
    {
      writeLines(token, connection, useBytes = TRUE)
      TRUE
    },
    error = function(...) FALSE
  )
  close(connection)
  if (!wrote) {
    unlink(path, force = TRUE)
    log_message(
      "Unable to create the {.pkg {backend}} result lock",
      message_type = "error"
    )
  }
  list(path = path, token = token)
}

runner_release_lock <- function(lock) {
  if (!file.exists(lock$path)) {
    return(invisible(FALSE))
  }
  observed <- tryCatch(
    readLines(lock$path, n = 1L, warn = FALSE),
    error = function(...) character()
  )
  if (!identical(observed, lock$token)) {
    return(invisible(FALSE))
  }
  invisible(unlink(lock$path, force = TRUE) == 0L)
}

runner_read_csv <- function(path, label, backend) {
  read_output <- function(...) {
    tryCatch(
      utils::read.csv(...),
      error = function(e) {
        log_message(
          "Unable to read {.pkg {backend}} {.val {label}} output: {conditionMessage(e)}",
          message_type = "error"
        )
      }
    )
  }
  header <- read_output(
    path,
    nrows = 0L,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(header) < 2L) {
    log_message(
      "{.pkg {backend}} returned invalid {.val {label}} output",
      message_type = "error"
    )
  }
  column_classes <- rep(NA_character_, ncol(header))
  column_classes[[1L]] <- "character"
  value <- read_output(
    path,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = column_classes,
    na.strings = character()
  )
  row_names <- value[[1L]]
  value <- value[-1L]
  value[] <- lapply(
    value,
    utils::type.convert,
    na.strings = "NA",
    as.is = TRUE
  )
  if (
    anyNA(row_names) ||
      any(!nzchar(row_names)) ||
      anyDuplicated(row_names)
  ) {
    log_message(
      "{.pkg {backend}} {.val {label}} must have unique non-empty row names",
      message_type = "error"
    )
  }
  rownames(value) <- row_names
  value
}

runner_read_numeric_csv <- function(path, label, backend) {
  value <- as.matrix(runner_read_csv(path, label, backend))
  suppressWarnings(storage.mode(value) <- "double")
  if (
    nrow(value) == 0L ||
      ncol(value) == 0L ||
      is.null(rownames(value)) ||
      is.null(colnames(value)) ||
      anyDuplicated(rownames(value)) ||
      anyDuplicated(colnames(value)) ||
      any(!is.finite(value))
  ) {
    log_message(
      "{.pkg {backend}} returned invalid {.val {label}} output",
      message_type = "error"
    )
  }
  value
}

runner_error <- function(
  status,
  stdout_path,
  stderr_path,
  backend,
  max_lines = 20L
) {
  stderr <- runner_tail_lines(stderr_path, max_lines = max_lines)
  stdout <- runner_tail_lines(stdout_path, max_lines = max_lines)
  details <- c(
    if (length(stderr)) {
      c("Python stderr:", stderr)
    },
    if (length(stdout)) {
      c("Python stdout:", stdout)
    }
  )
  if (!length(details)) {
    details <- "<no Python output captured>"
  }
  log_message(
    "{.pkg {backend}} Python runner failed with status {.val {status}}:\n{.code {paste(details, collapse = '\n')}}",
    message_type = "error"
  )
}

runner_tail_lines <- function(path, max_lines = 20L, chunk_size = 1000L) {
  if (!file.exists(path)) {
    return(character())
  }
  connection <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  output <- character()
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (!length(lines)) {
      break
    }
    lines <- lines[nzchar(trimws(lines))]
    output <- utils::tail(c(output, lines), max_lines)
  }
  output
}
