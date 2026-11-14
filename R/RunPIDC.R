#' @rdname RunGRN-methods
#' @export
RunPIDC <- function(object, ...) {
  RunGRN(object, methods = "pidc", ...)
}

infer_pidc <- function(expression, regulators, targets, arguments, config) {
  check_r("PIDC", strict = TRUE, install = FALSE, verbose = FALSE)
  pidc <- get_namespace_fun("PIDC", "PIDC")
  nonconstant <- apply(expression, 1L, stats::sd) > 0
  expression <- expression[nonconstant, , drop = FALSE]
  regulators <- intersect(regulators, rownames(expression))
  targets <- intersect(targets, rownames(expression))
  if (length(regulators) < 2L || !length(targets)) {
    stop("PIDC needs at least two nonconstant regulators and one target.", call. = FALSE)
  }
  log_scale <- arguments$log_scale %||% FALSE
  if (!is.logical(log_scale) || length(log_scale) != 1L || is.na(log_scale)) {
    stop("PIDC log_scale must be TRUE or FALSE.", call. = FALSE)
  }
  fit <- pidc(
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
    implementation = "PIDC::PIDC",
    package_versions = c(PIDC = as.character(utils::packageVersion("PIDC")))
  )
}
