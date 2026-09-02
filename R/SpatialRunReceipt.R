# Internal spatial run receipt renderer --------------------------------------

spatial_run_receipt_lines <- function(x, name) {
  if (is.null(x)) {
    return(character())
  }
  if (!is.character(x)) {
    stop(sprintf("%s must be NULL or a character vector", name), call. = FALSE)
  }
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(character())
  }
  as.character(x)
}

spatial_run_receipt_quote <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("%s must be one non-empty character string", name), call. = FALSE)
  }
  deparse1(x, width.cutoff = 500L)
}

spatial_run_receipt <- function(
  done,
  scope = NULL,
  saved = NULL,
  plot = NULL,
  inspect = NULL,
  status = c("completed", "partial"),
  replaced = FALSE,
  verbose = TRUE,
  .envir = parent.frame()
) {
  verbose <- thisutils::get_verbose(verbose)
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  status <- match.arg(status)
  if (!is.character(done) || length(done) != 1L || is.na(done) || !nzchar(done)) {
    stop("done must be one non-empty character string", call. = FALSE)
  }
  scope <- spatial_run_receipt_lines(scope, "scope")
  saved <- spatial_run_receipt_lines(saved, "saved")
  plot <- spatial_run_receipt_lines(plot, "plot")
  inspect <- spatial_run_receipt_lines(inspect, "inspect")
  if (!is.logical(replaced) || length(replaced) != 1L || is.na(replaced)) {
    stop("replaced must be one non-missing logical value", call. = FALSE)
  }
  if (!is.environment(.envir)) {
    stop(".envir must be an environment", call. = FALSE)
  }

  message_type <- if (identical(status, "completed")) "success" else "running"
  log_message(
    done,
    message_type = message_type,
    verbose = TRUE,
    .envir = .envir
  )

  if (length(scope) > 0L) {
    log_message(
      paste(c("{.field Scope}", scope), collapse = " "),
      level = 2,
      timestamp = FALSE,
      verbose = TRUE,
      .envir = .envir
    )
  }

  if (length(saved) > 0L) {
    saved_label <- if (isTRUE(replaced)) "{.field Replaced}" else "{.field Saved}"
    log_message(
      paste(c(saved_label, saved), collapse = " "),
      level = 2,
      timestamp = FALSE,
      verbose = TRUE,
      .envir = .envir
    )
  }

  next_lines <- if (length(plot) > 0L) plot else inspect
  next_label <- if (length(plot) > 0L) "{.field Plot}" else "{.field Inspect}"
  if (length(next_lines) > 0L) {
    next_call <- next_lines[[1L]]
    log_message(
      paste(next_label, "{.code {next_call}}"),
      level = 2,
      timestamp = FALSE,
      verbose = TRUE,
      .envir = environment()
    )
  }

  invisible(NULL)
}
