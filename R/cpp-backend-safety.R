cpp_dense_gib <- function(n_rows, n_cols, copies = 1) {
  values <- c(n_rows, n_cols, copies)
  if (
    length(values) != 3L ||
      any(!is.finite(values)) ||
      any(values < 0)
  ) {
    stop("Dense-memory dimensions and copies must be finite non-negative values.", call. = FALSE)
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
    stop("max_dense_gib must be one positive number or Inf.", call. = FALSE)
  }
  estimated_gib <- cpp_dense_gib(n_rows, n_cols, copies)
  if (is.finite(max_dense_gib) && estimated_gib > max_dense_gib) {
    stop(
      sprintf(
        "%s would require at least %.2f GiB for dense working matrices, exceeding max_dense_gib = %.2f. Use a smaller input, choose the reference backend, or explicitly increase max_dense_gib.",
        context,
        estimated_gib,
        max_dense_gib
      ),
      call. = FALSE
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
    stop(
      sprintf(
        "%s is an approximate implementation. Set allow_approximate = TRUE to opt in, or use backend = \"%s\" for the reference workflow.",
        context,
        reference_backend
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

reject_unsupported_cpp_arguments <- function(arguments, context) {
  arguments <- unique(as.character(arguments))
  arguments <- arguments[nzchar(arguments)]
  if (length(arguments) > 0L) {
    stop(
      sprintf(
        "%s does not support: %s. Use the reference backend for these arguments.",
        context,
        paste(arguments, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
