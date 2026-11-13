run_leap <- function(expression, pseudotime, max_lag_prop = 1 / 3, MAC_cutoff = 0) {
  if (!requireNamespace("LEAP", quietly = TRUE)) {
    stop("Install LEAP to run this method.", call. = FALSE)
  }
  if (length(max_lag_prop) != 1L || !is.finite(max_lag_prop) ||
    max_lag_prop <= 0 || max_lag_prop >= 1 ||
    length(MAC_cutoff) != 1L || !is.finite(MAC_cutoff) || MAC_cutoff < 0) {
    stop("LEAP requires 0 < max_lag_prop < 1 and a nonnegative MAC_cutoff.", call. = FALSE)
  }

  branches <- lapply(seq_len(ncol(pseudotime)), function(branch) {
    cells <- which(is.finite(pseudotime[, branch]))
    cells[order(pseudotime[cells, branch], colnames(expression)[cells], method = "radix")]
  })
  branches <- unique(branches[lengths(branches) >= 3L])
  if (!length(branches)) {
    stop("LEAP needs a trajectory with at least three finite-time cells.", call. = FALSE)
  }

  fits <- lapply(branches, function(cells) {
    LEAP::MAC_counter(
      expression[, cells, drop = FALSE],
      max_lag_prop = max_lag_prop,
      MAC_cutoff = MAC_cutoff,
      lag_matrix = FALSE,
      symmetric = FALSE
    )
  })
  ordered_cells <- lapply(branches, function(cells) colnames(expression)[cells])
  tables <- lapply(fits, function(fit) {
    fit <- as.matrix(fit)
    if (!nrow(fit)) {
      return(data.frame(
        regulator = character(), target = character(), weight = numeric(),
        lag = numeric(), signed_weight = numeric()
      ))
    }
    data.frame(
      regulator = rownames(expression)[fit[, "Row gene index"]],
      target = rownames(expression)[fit[, "Column gene index"]],
      weight = as.numeric(fit[, "Correlation"]),
      lag = as.numeric(fit[, "Lag"]),
      signed_weight = as.numeric(fit[, "Correlation"])
    )
  })
  edges <- do.call(rbind, tables)
  edges <- edges[is.finite(edges$weight), , drop = FALSE]
  edges <- edges[order(-abs(edges$weight), method = "radix"), , drop = FALSE]
  edges <- edges[!duplicated(edges[c("regulator", "target")]), , drop = FALSE]

  list(
    branches = fits,
    ordered_cells = ordered_cells,
    implementation = "LEAP::MAC_counter",
    package_versions = c(LEAP = as.character(utils::packageVersion("LEAP"))),
    edges = edges,
    parameters = list(max_lag_prop = max_lag_prop, MAC_cutoff = MAC_cutoff, branches = ordered_cells),
    genes = rownames(expression),
    cells = colnames(expression)
  )
}
