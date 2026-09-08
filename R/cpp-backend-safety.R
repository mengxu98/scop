# Map an R `cores` / `n_threads` argument onto C++ OpenMP `n_threads`.
# `NULL` or a non-positive value becomes 0, which means "use the process
# OpenMP default" (`omp_get_max_threads()`, i.e. OMP_NUM_THREADS).
scop_n_threads <- function(cores = NULL) {
  if (is.null(cores) || length(cores) < 1L) {
    return(0L)
  }
  n <- suppressWarnings(as.integer(cores[[1L]]))
  if (length(n) != 1L || is.na(n) || n < 0L) {
    return(0L)
  }
  n
}

# OpenMP team for a C++ kernel that is the only parallel region in the call
# (one R worker). Public `cores = 1` means "one R worker", not "one OpenMP
# thread": NULL or cores <= 1 uses the process OpenMP default. cores >= 2
# caps the team. Nested R workers must pass the inner budget through
# scop_n_threads() instead, so cores = 1 stays serial inside each worker.
scop_inner_n_threads <- function(cores = NULL) {
  if (is.null(cores) || length(cores) < 1L) {
    return(0L)
  }
  n <- suppressWarnings(as.integer(cores[[1L]]))
  if (length(n) != 1L || is.na(n) || n <= 1L) {
    return(0L)
  }
  n
}

cpp_dense_gib <- function(n_rows, n_cols, copies = 1) {
  values <- c(n_rows, n_cols, copies)
  if (
    length(values) != 3L ||
      any(!is.finite(values)) ||
      any(values < 0)
  ) {
    log_message("Dense-memory dimensions and copies must be finite non-negative values.", message_type = "error")
  }
  as.numeric(n_rows) * as.numeric(n_cols) * 8 * as.numeric(copies) / 1024^3
}

assert_cpp_dense_budget <- function(
  n_rows,
  n_cols,
  copies,
  max_dense_gib,
  context
) {
  if (
    length(max_dense_gib) != 1L ||
      is.na(max_dense_gib) ||
      max_dense_gib <= 0
  ) {
    log_message("max_dense_gib must be one positive number or Inf.", message_type = "error")
  }
  estimated_gib <- cpp_dense_gib(n_rows, n_cols, copies)
  if (is.finite(max_dense_gib) && estimated_gib > max_dense_gib) {
    log_message(
      sprintf(
        "%s would require at least %.2f GiB for dense working matrices, exceeding max_dense_gib = %.2f. Use a smaller input, choose the reference backend, or explicitly increase max_dense_gib.",
        context,
        estimated_gib,
        max_dense_gib
      ),
      message_type = "error"
    )
  }
  estimated_gib
}

assert_cpp_approximation_opt_in <- function(
  allow_approximate,
  context,
  reference_backend = "python"
) {
  if (!isTRUE(allow_approximate)) {
    log_message(
      sprintf(
        "%s is an approximate implementation. Set allow_approximate = TRUE to opt in, or use backend = \"%s\" for the reference workflow.",
        context,
        reference_backend
      ),
      message_type = "error"
    )
  }
  invisible(TRUE)
}

reject_unsupported_cpp_arguments <- function(arguments, context) {
  arguments <- unique(as.character(arguments))
  arguments <- arguments[nzchar(arguments)]
  if (length(arguments) > 0L) {
    log_message(
      sprintf(
        "%s does not support: %s. Use the reference backend for these arguments.",
        context,
        paste(arguments, collapse = ", ")
      ),
      message_type = "error"
    )
  }
  invisible(TRUE)
}
