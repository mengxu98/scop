run_pidc <- function(expression, regulators, targets, arguments, config) {
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    stop("Install the pbapply package before running pidc.", call. = FALSE)
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Install purrr before running PIDC.", call. = FALSE)
  }
  nonconstant <- apply(expression, 1L, stats::sd) > 0
  expression <- expression[nonconstant, , drop = FALSE]
  regulators <- intersect(regulators, rownames(expression))
  targets <- intersect(targets, rownames(expression))
  if (length(regulators) < 2L || !length(targets)) {
    stop("PIDC needs at least two nonconstant regulators and one target.", call. = FALSE)
  }
  upstream <- new.env(parent = baseenv())
  upstream$is <- methods::is
  sys.source(config$external_sources[["entry"]], envir = upstream)
  log_scale <- arguments$log_scale %||% FALSE
  if (!is.logical(log_scale) || length(log_scale) != 1L || is.na(log_scale)) {
    stop("PIDC log_scale must be TRUE or FALSE.", call. = FALSE)
  }
  fit <- upstream$PIDC(
    expression,
    regulators = regulators,
    targets = targets,
    logScale = log_scale,
    ncores = config$cores,
    diag = "zero",
    verbose = config$verbose
  )
  list(
    edges = grn_rna_matrix_edges(fit),
    parameters = list(log_scale = log_scale),
    implementation = "multiCSN-manuscript::PIDC-R",
    package_versions = c(
      pbapply = as.character(utils::packageVersion("pbapply")),
      purrr = as.character(utils::packageVersion("purrr"))
    )
  )
}
