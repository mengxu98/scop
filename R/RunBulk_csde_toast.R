RunBulk_csde_toast <- function(
  count_matrix,
  condition,
  proportions,
  condition1 = NULL,
  condition2 = NULL,
  backend = c("limma_interaction", "native"),
  p.adjust.method = "BH",
  min_prop_variance = 1e-8,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  if (identical(backend, "native")) {
    return(.bulk_method_failed(
      reason = paste(
        "The native TOAST backend is not implemented yet.",
        "Use backend = 'limma_interaction' for SCOP internal CSDE fitting."
      ),
      parameters = list(backend = backend),
      results = .bulk_empty_csde_result()
    ))
  }
  pkg_check <- .bulk_check_r_packages("limma", "limma_interaction backend")
  if (!is.null(pkg_check)) {
    return(pkg_check)
  }
  lmFit <- get_namespace_fun("limma", "lmFit")
  eBayes <- get_namespace_fun("limma", "eBayes")
  topTable <- get_namespace_fun("limma", "topTable")
  log_message(
    "{.val csde_TOAST} currently uses a limma interaction backend.",
    verbose = verbose
  )

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

  cond <- as.character(condition)
  names(cond) <- names(condition) %||% colnames(count_matrix)
  sample_keep <- intersect(
    colnames(count_matrix)[cond[colnames(count_matrix)] %in% c(pair$condition1, pair$condition2)],
    rownames(proportions)
  )
  if (length(sample_keep) < 4) {
    return(.bulk_method_failed(
      reason = "TOAST requires at least 4 matched samples between counts and proportions.",
      details = list(condition_pair = pair)
    ))
  }

  count_use <- as.matrix(count_matrix[, sample_keep, drop = FALSE])
  cond_use <- factor(cond[sample_keep], levels = c(pair$condition1, pair$condition2))
  if (any(table(cond_use) < 2)) {
    return(.bulk_method_failed(
      reason = "TOAST requires at least 2 samples per condition.",
      details = list(condition_pair = pair)
    ))
  }
  prop_use <- as.matrix(proportions[sample_keep, , drop = FALSE])
  expr_use <- log1p(count_use)

  csde_results <- list()
  for (ct in colnames(prop_use)) {
    prop_vec <- as.numeric(prop_use[, ct])
    if (stats::var(prop_vec, na.rm = TRUE) <= min_prop_variance) {
      next
    }

    design <- stats::model.matrix(
      ~ prop + cond,
      data = data.frame(
        prop = prop_vec,
        cond = cond_use
      )
    )
    interaction_col <- prop_vec * as.numeric(cond_use == pair$condition2)
    design <- cbind(design, prop_cond = interaction_col)
    fit <- tryCatch(
      lmFit(expr_use, design),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      next
    }
    fit <- eBayes(fit)
    if (!"prop_cond" %in% colnames(fit$coefficients)) {
      next
    }
    tt <- topTable(
      fit = fit,
      coef = "prop_cond",
      number = Inf,
      sort.by = "none"
    )
    if (nrow(tt) == 0) {
      next
    }
    out_df <- data.frame(
      gene = rownames(tt),
      cell_type = ct,
      group1 = pair$condition2,
      group2 = pair$condition1,
      effect = tt$logFC,
      p_val = tt$P.Value,
      p_val_adj = stats::p.adjust(tt$P.Value, method = p.adjust.method),
      method = "TOAST",
      stringsAsFactors = FALSE
    )
    csde_results[[ct]] <- out_df
  }

  if (length(csde_results) == 0) {
    return(.bulk_method_failed(
      reason = "No valid CSDE result was produced for method 'TOAST'.",
      details = list(condition_pair = pair),
      results = .bulk_empty_csde_result()
    ))
  }

  res <- do.call(rbind, csde_results)
  rownames(res) <- NULL
  .bulk_method_success(
    results = res,
    details = list(
      condition_pair = pair,
      engine = "limma_interaction",
      package_backend = FALSE,
      note = paste(
        "This runner currently uses a limma interaction model rather than",
        "TOAST's native csDE API."
      )
    ),
    parameters = list(
      backend = backend,
      p.adjust.method = p.adjust.method,
      min_prop_variance = min_prop_variance
    )
  )
}
