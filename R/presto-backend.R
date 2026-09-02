# Runtime-optional Presto backend helpers ------------------------------------

.presto_repository <- "immunogenomics/presto"
.presto_package <- "presto"

presto_check_r <- function(install = FALSE) {
  if (!is.logical(install) || length(install) != 1L || is.na(install)) {
    log_message("{.arg install} must be TRUE or FALSE", message_type = "error")
  }
  package_request <- if (isTRUE(install)) {
    .presto_repository
  } else {
    .presto_package
  }
  status <- check_r(
    package_request,
    install = install,
    verbose = FALSE
  )
  status <- unlist(status, recursive = TRUE, use.names = TRUE)
  if (.presto_package %in% names(status)) {
    return(isTRUE(status[[.presto_package]]))
  }
  length(status) == 1L && isTRUE(status[[1L]])
}

presto_get_fun <- function(
  fun = "wilcoxauc",
  install = FALSE,
  error_on_missing = TRUE
) {
  if (!is.character(fun) || length(fun) != 1L || is.na(fun) || !nzchar(fun)) {
    log_message("{.arg fun} must be one non-empty string", message_type = "error")
  }
  if (!is.logical(error_on_missing) ||
    length(error_on_missing) != 1L ||
    is.na(error_on_missing)) {
    log_message(
      "{.arg error_on_missing} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  if (!isTRUE(presto_check_r(install = install))) {
    if (!isTRUE(error_on_missing)) {
      return(NULL)
    }
    log_message(
      "The optional {.pkg presto} backend is unavailable",
      message_type = "error"
    )
  }
  resolved <- tryCatch(
    get_namespace_fun(.presto_package, fun),
    error = identity
  )
  if (inherits(resolved, "error")) {
    if (!isTRUE(error_on_missing)) {
      return(NULL)
    }
    log_message(
      "The optional {.pkg presto} backend does not provide {.fn {fun}}: {conditionMessage(resolved)}",
      message_type = "error"
    )
  }
  resolved
}
