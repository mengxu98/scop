RunBulk_de_edgeR <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = FALSE,
  logfc.threshold = 0,
  p.adjust.method = "BH",
  verbose = TRUE
) {
  pkg_check <- .bulk_check_r_packages("edgeR", "method 'edgeR'")
  if (!is.null(pkg_check)) {
    return(pkg_check)
  }

  pair <- tryCatch(
    .bulk_resolve_condition_pair(
      condition = condition,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    ),
    error = function(e) e
  )
  if (inherits(pair, "error")) {
    return(.bulk_method_failed(reason = pair$message))
  }

  out <- tryCatch(
    RunDEtest_edgeR(
      count_matrix = count_matrix,
      condition = condition,
      condition1 = pair$condition1,
      condition2 = pair$condition2,
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    ),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(.bulk_method_failed(reason = out$message))
  }
  if (is.null(out) || nrow(out) == 0) {
    return(.bulk_method_success(
      results = data.frame(),
      details = list(condition_pair = pair),
      parameters = list(
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        p.adjust.method = p.adjust.method
      )
    ))
  }

  .bulk_method_success(
    results = out,
    details = list(condition_pair = pair),
    parameters = list(
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    )
  )
}
