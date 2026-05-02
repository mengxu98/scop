RunBulk_de_deseq2 <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = FALSE,
  logfc.threshold = 0,
  p.adjust.method = "BH",
  verbose = TRUE
) {
  pkg_check <- .bulk_check_r_packages("DESeq2", "method 'DESeq2'")
  if (!is.null(pkg_check)) {
    return(pkg_check)
  }
  DESeqDataSetFromMatrix <- get_namespace_fun("DESeq2", "DESeqDataSetFromMatrix")
  DESeq <- get_namespace_fun("DESeq2", "DESeq")
  DESeq_results <- get_namespace_fun("DESeq2", "results")

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
      reason = "DESeq2 requires at least 2 samples per condition.",
      details = list(condition_pair = pair)
    ))
  }

  out <- tryCatch(
    {
      col_data <- data.frame(
        condition = cond_use,
        row.names = colnames(count_use),
        stringsAsFactors = FALSE
      )
      dds <- DESeqDataSetFromMatrix(
        countData = round(count_use),
        colData = col_data,
        design = ~condition
      )
      dds <- DESeq(dds, quiet = TRUE)
      res <- DESeq_results(
        dds,
        contrast = c("condition", pair$condition2, pair$condition1)
      )
      res_df <- as.data.frame(res)
      if (nrow(res_df) == 0) {
        data.frame()
      } else {
        detect <- count_use[rownames(res_df), , drop = FALSE] > 0
        pct.1 <- round(rowMeans(detect[, cond_use == pair$condition1, drop = FALSE]), 3)
        pct.2 <- round(rowMeans(detect[, cond_use == pair$condition2, drop = FALSE]), 3)

        out_df <- data.frame(
          p_val = res_df$pvalue,
          avg_log2FC = res_df$log2FoldChange,
          ave_expr = log10(res_df$baseMean + 1),
          pct.1 = pct.1,
          pct.2 = pct.2,
          p_val_adj = res_df$padj,
          row.names = rownames(res_df)
        )
        out_df <- out_df[!is.na(out_df$p_val), , drop = FALSE]
        out_df$p_val_adj <- stats::p.adjust(out_df$p_val, method = p.adjust.method)
        if (isTRUE(only.pos)) {
          out_df <- out_df[out_df$avg_log2FC >= logfc.threshold, , drop = FALSE]
        } else {
          out_df <- out_df[abs(out_df$avg_log2FC) >= logfc.threshold, , drop = FALSE]
        }
        out_df
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
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    )
  )
}
