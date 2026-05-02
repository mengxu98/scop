RunBulk_de_dream <- function(
  count_matrix,
  condition,
  sample_data = NULL,
  condition1 = NULL,
  condition2 = NULL,
  dream_formula = NULL,
  only.pos = FALSE,
  logfc.threshold = 0,
  p.adjust.method = "BH",
  verbose = TRUE
) {
  pkg_check <- .bulk_check_r_packages(
    c("variancePartition", "limma", "edgeR"),
    "method 'dream'"
  )
  if (!is.null(pkg_check)) {
    return(pkg_check)
  }
  DGEList <- get_namespace_fun("edgeR", "DGEList")
  filterByExpr <- get_namespace_fun("edgeR", "filterByExpr")
  calcNormFactors <- get_namespace_fun("edgeR", "calcNormFactors")
  voomWithDreamWeights <- get_namespace_fun("variancePartition", "voomWithDreamWeights")
  dream <- get_namespace_fun("variancePartition", "dream")
  eBayes <- get_namespace_fun("limma", "eBayes")
  topTable <- get_namespace_fun("limma", "topTable")

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
  keep <- cond %in% c(pair$condition1, pair$condition2)
  count_use <- as.matrix(count_matrix[, keep, drop = FALSE])
  cond_use <- factor(cond[keep], levels = c(pair$condition1, pair$condition2))
  if (ncol(count_use) < 2 || any(table(cond_use) < 2)) {
    return(.bulk_method_failed(
      reason = "dream requires at least 2 samples per condition.",
      details = list(condition_pair = pair)
    ))
  }

  sample_df <- if (is.null(sample_data)) {
    data.frame(row.names = colnames(count_use))
  } else {
    if (is.null(rownames(sample_data)) || any(rownames(sample_data) == "")) {
      return(.bulk_method_failed(
        reason = "{.arg sample_data} must contain row names matching sample IDs."
      ))
    }
    if (!all(colnames(count_use) %in% rownames(sample_data))) {
      return(.bulk_method_failed(
        reason = "{.arg sample_data} row names must cover all bulk sample IDs."
      ))
    }
    sample_data[colnames(count_use), , drop = FALSE]
  }
  sample_df$condition <- cond_use
  sample_df$sample <- rownames(sample_df)

  if (is.null(dream_formula)) {
    dream_formula <- stats::as.formula("~ condition")
  } else if (is.character(dream_formula) && length(dream_formula) == 1) {
    dream_formula <- stats::as.formula(dream_formula)
  }
  if (!inherits(dream_formula, "formula")) {
    return(.bulk_method_failed(
      reason = "{.arg dream_formula} must be a formula or a one-length character."
    ))
  }
  if (!"condition" %in% all.vars(dream_formula)) {
    return(.bulk_method_failed(
      reason = "{.arg dream_formula} must include the {.val condition} term."
    ))
  }

  out <- tryCatch(
    {
      dge <- DGEList(counts = count_use)
      keep_features <- filterByExpr(dge, group = sample_df$condition)
      if (!any(keep_features)) {
        data.frame()
      } else {
        dge <- dge[keep_features, , keep.lib.sizes = FALSE]
        dge <- calcNormFactors(dge)

        vobj <- voomWithDreamWeights(
          dge,
          formula = dream_formula,
          data = sample_df
        )
        fit <- dream(
          exprObj = vobj,
          formula = dream_formula,
          data = sample_df
        )
        fit <- eBayes(fit)
        coef_names <- colnames(fit$coefficients)
        coef_use <- grep(
          pattern = paste0("^condition", make.names(pair$condition2), "$"),
          x = coef_names,
          value = TRUE
        )
        if (length(coef_use) == 0) {
          coef_use <- grep("^condition", coef_names, value = TRUE)
        }
        if (length(coef_use) == 0) {
          stop("Cannot find condition coefficient in dream model.")
        }
        tt <- topTable(
          fit = fit,
          coef = coef_use[[1]],
          number = Inf,
          sort.by = "none"
        )
        if (nrow(tt) == 0) {
          data.frame()
        } else {
          detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
          pct.1 <- round(rowMeans(detect[, sample_df$condition == pair$condition1, drop = FALSE]), 3)
          pct.2 <- round(rowMeans(detect[, sample_df$condition == pair$condition2, drop = FALSE]), 3)
          out_df <- data.frame(
            p_val = tt$P.Value,
            avg_log2FC = tt$logFC,
            ave_expr = tt$AveExpr,
            pct.1 = pct.1,
            pct.2 = pct.2,
            row.names = rownames(tt)
          )
          out_df$p_val_adj <- stats::p.adjust(out_df$p_val, method = p.adjust.method)
          if (isTRUE(only.pos)) {
            out_df <- out_df[out_df$avg_log2FC >= logfc.threshold, , drop = FALSE]
          } else {
            out_df <- out_df[abs(out_df$avg_log2FC) >= logfc.threshold, , drop = FALSE]
          }
          out_df
        }
      }
    },
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(.bulk_method_failed(reason = out$message))
  }

  .bulk_method_success(
    results = out,
    details = list(condition_pair = pair),
    parameters = list(
      dream_formula = deparse(dream_formula),
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    )
  )
}
