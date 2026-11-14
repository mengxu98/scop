#' @rdname RunGRN-methods
#' @export
RunPPCOR <- function(object, ...) {
  RunGRN(object, methods = "ppcor", ...)
}

infer_ppcor <- function(expression, arguments) {
  check_r("ppcor", strict = TRUE, verbose = FALSE)
  method <- arguments$method %||% "spearman"
  if (!method %in% c("spearman", "pearson", "kendall")) {
    stop("Invalid ppcor method.", call. = FALSE)
  }
  expression <- expression[apply(expression, 1L, stats::sd) > 0, , drop = FALSE]
  if (nrow(expression) < 2L) {
    stop("PPCOR needs at least two nonconstant genes.", call. = FALSE)
  }
  fit <- get_namespace_fun("ppcor", "pcor")(t(expression), method = method)
  dimnames(fit$estimate) <- list(rownames(expression), rownames(expression))
  edges <- grn_rna_matrix_edges(fit$estimate)
  index <- cbind(
    match(edges$regulator, rownames(fit$estimate)),
    match(edges$target, colnames(fit$estimate))
  )
  edges$p_value <- as.numeric(fit$p.value[index])
  edges$statistic <- as.numeric(fit$statistic[index])
  list(
    edges = edges,
    extra = list(sample_size = fit$n, variables = fit$gp, method = fit$method),
    parameters = list(method = method),
    implementation = "ppcor::pcor",
    package_versions = c(ppcor = as.character(utils::packageVersion("ppcor")))
  )
}
