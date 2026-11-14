infer_ppcor <- function(expression, arguments) {
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop("Install the ppcor package before running ppcor.", call. = FALSE)
  }
  method <- arguments$method %||% "spearman"
  if (!method %in% c("spearman", "pearson", "kendall")) {
    stop("Invalid ppcor method.", call. = FALSE)
  }
  fit <- ppcor::pcor(t(expression), method = method)
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
