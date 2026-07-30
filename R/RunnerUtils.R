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
  writeLines(
    as.character(to_json(
      x,
      auto_unbox = TRUE,
      null = "null",
      digits = NA,
      pretty = TRUE
    )),
    con = path,
    useBytes = TRUE
  )
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
  old_env <- Sys.getenv(names(env), unset = NA_character_)
  names(old_env) <- names(env)
  if (length(env)) {
    do.call(Sys.setenv, as.list(env))
  }
  on.exit({
    restore <- !is.na(old_env)
    if (any(restore)) {
      do.call(Sys.setenv, as.list(old_env[restore]))
    }
    if (any(!restore)) {
      Sys.unsetenv(names(old_env)[!restore])
    }
  }, add = TRUE)
  system2(
    command = command,
    args = args,
    stdout = stdout,
    stderr = stderr
  )
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
  read_output <- function(path) {
    if (!file.exists(path)) {
      return(character())
    }
    lines <- readLines(path, warn = FALSE)
    lines[nzchar(trimws(lines))]
  }
  stderr <- read_output(stderr_path)
  stdout <- read_output(stdout_path)
  details <- c(
    if (length(stderr)) {
      c("Python stderr:", utils::tail(stderr, max_lines))
    },
    if (length(stdout)) {
      c("Python stdout:", utils::tail(stdout, max_lines))
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
