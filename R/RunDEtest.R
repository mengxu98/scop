# Internal sample-level and bulk DE method implementations for RunDEtest().
NULL

RunLimmaVoom <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = FALSE,
  logfc.threshold = 0,
  p.adjust.method = "BH",
  verbose = TRUE
) {
  check_r(c("limma", "edgeR"), verbose = FALSE)

  pair <- tryCatch(
    resolve_condition_pair(
      condition = condition,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    ),
    error = function(e) e
  )
  if (inherits(pair, "error")) {
    return(list(
      status = "failed",
      reason = pair$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }

  out <- tryCatch(
    RunDEtest_limma(
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
    return(list(
      status = "failed",
      reason = out$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  if (is.null(out) || nrow(out) == 0) {
    return(list(
      status = "success",
      reason = NULL,
      results = data.frame(),
      details = list(condition_pair = pair),
      parameters = list(
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        p.adjust.method = p.adjust.method
      )
    ))
  }

  list(
    status = "success",
    reason = NULL,
    results = out,
    details = list(condition_pair = pair),
    parameters = list(
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    )
  )
}


RunEdgeR <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = FALSE,
  logfc.threshold = 0,
  p.adjust.method = "BH",
  verbose = TRUE
) {
  check_r("edgeR", verbose = FALSE)

  pair <- tryCatch(
    resolve_condition_pair(
      condition = condition,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    ),
    error = function(e) e
  )
  if (inherits(pair, "error")) {
    return(list(
      status = "failed",
      reason = pair$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
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
    return(list(
      status = "failed",
      reason = out$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  if (is.null(out) || nrow(out) == 0) {
    return(list(
      status = "success",
      reason = NULL,
      results = data.frame(),
      details = list(condition_pair = pair),
      parameters = list(
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        p.adjust.method = p.adjust.method
      )
    ))
  }

  list(
    status = "success",
    reason = NULL,
    results = out,
    details = list(condition_pair = pair),
    parameters = list(
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    )
  )
}


RunDESeq2 <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = FALSE,
  logfc.threshold = 0,
  p.adjust.method = "BH",
  verbose = TRUE
) {
  check_r("DESeq2", verbose = FALSE)

  pair <- tryCatch(
    resolve_condition_pair(
      condition = condition,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    ),
    error = function(e) e
  )
  if (inherits(pair, "error")) {
    return(list(
      status = "failed",
      reason = pair$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }

  cond <- as.character(condition)
  keep <- cond %in% c(pair$condition1, pair$condition2)
  count_use <- as.matrix(count_matrix[, keep, drop = FALSE])
  cond_use <- factor(cond[keep], levels = c(pair$condition1, pair$condition2))
  if (ncol(count_use) < 2 || any(table(cond_use) < 2)) {
    return(list(
      status = "failed",
      reason = "DESeq2 requires at least 2 samples per condition.",
      details = list(condition_pair = pair),
      results = data.frame(),
      parameters = list()
    ))
  }

  out <- tryCatch(
    RunDEtest_DESeq2(
      count_matrix = count_use,
      condition = cond_use,
      condition1 = pair$condition1,
      condition2 = pair$condition2,
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    ),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(list(
      status = "failed",
      reason = out$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  if (is.null(out) || nrow(out) == 0) {
    return(list(
      status = "success",
      reason = NULL,
      results = data.frame(),
      details = list(condition_pair = pair),
      parameters = list(
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        p.adjust.method = p.adjust.method
      )
    ))
  }

  list(
    status = "success",
    reason = NULL,
    results = out,
    details = list(condition_pair = pair),
    parameters = list(
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    )
  )
}


RunDream <- function(
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
  check_r(c("variancePartition", "limma", "edgeR"), verbose = FALSE)

  pair <- tryCatch(
    resolve_condition_pair(
      condition = condition,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    ),
    error = function(e) e
  )
  if (inherits(pair, "error")) {
    return(list(
      status = "failed",
      reason = pair$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }

  cond <- as.character(condition)
  keep <- cond %in% c(pair$condition1, pair$condition2)
  count_use <- as.matrix(count_matrix[, keep, drop = FALSE])
  cond_use <- factor(cond[keep], levels = c(pair$condition1, pair$condition2))
  if (ncol(count_use) < 2 || any(table(cond_use) < 2)) {
    return(list(
      status = "failed",
      reason = "dream requires at least 2 samples per condition.",
      details = list(condition_pair = pair),
      results = data.frame(),
      parameters = list()
    ))
  }

  if (!is.null(sample_data)) {
    if (is.null(rownames(sample_data)) || any(rownames(sample_data) == "")) {
      return(list(
        status = "failed",
        reason = "{.arg sample_data} must contain row names matching sample IDs.",
        results = data.frame(),
        details = list(),
        parameters = list()
      ))
    }
    if (!all(colnames(count_use) %in% rownames(sample_data))) {
      return(list(
        status = "failed",
        reason = "{.arg sample_data} row names must cover all bulk sample IDs.",
        results = data.frame(),
        details = list(),
        parameters = list()
      ))
    }
  }

  if (is.null(dream_formula)) {
    dream_formula <- stats::as.formula("~ condition")
  } else if (is.character(dream_formula) && length(dream_formula) == 1) {
    dream_formula <- stats::as.formula(dream_formula)
  }
  if (!inherits(dream_formula, "formula")) {
    return(list(
      status = "failed",
      reason = "{.arg dream_formula} must be a formula or a one-length character.",
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  if (!"condition" %in% all.vars(dream_formula)) {
    return(list(
      status = "failed",
      reason = "{.arg dream_formula} must include the {.val condition} term.",
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }

  out <- tryCatch(
    RunDEtest_dream(
      count_matrix = count_use,
      condition = cond_use,
      condition1 = pair$condition1,
      condition2 = pair$condition2,
      sample_data = sample_data,
      dream_formula = dream_formula,
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      p.adjust.method = p.adjust.method
    ),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(list(
      status = "failed",
      reason = out$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  if (is.null(out) || nrow(out) == 0) {
    return(list(
      status = "success",
      reason = NULL,
      results = data.frame(),
      details = list(condition_pair = pair),
      parameters = list(
        dream_formula = deparse(dream_formula),
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        p.adjust.method = p.adjust.method
      )
    ))
  }

  list(
    status = "success",
    reason = NULL,
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


dispatch_de <- function(
  ctx,
  method_name,
  condition1 = NULL,
  condition2 = NULL,
  de_markers_type = "single",
  de_args = list(),
  verbose = TRUE
) {
  method_map <- list(
    de_limma_voom = RunLimmaVoom,
    de_edgeR_qlf = RunEdgeR,
    de_DESeq2 = RunDESeq2,
    de_dream = RunDream
  )
  method_fun <- method_map[[method_name]]
  if (is.null(method_fun)) {
    log_message(
      paste0("Unsupported DE method: ", method_name),
      message_type = "error"
    )
  }

  per_group <- list()
  results_list <- list()
  fail_messages <- character(0)

  for (grp in names(ctx$counts_by_group)) {
    condition_vec <- ctx$condition_by_group[[grp]]
    comparison_specs <- tryCatch(
      prepare_de_comparisons(
        condition = condition_vec,
        condition1 = condition1,
        condition2 = condition2,
        markers_type = de_markers_type
      ),
      error = function(e) e
    )
    if (inherits(comparison_specs, "error")) {
      fail_messages <- c(
        fail_messages,
        paste0(grp, ": ", comparison_specs$message)
      )
      next
    }

    group_bundles <- list()
    for (spec in comparison_specs) {
      args <- utils::modifyList(
        list(
          count_matrix = ctx$counts_by_group[[grp]],
          condition = spec$condition,
          sample_data = ctx$sample_meta_by_group[[grp]],
          condition1 = spec$condition1,
          condition2 = spec$condition2,
          verbose = verbose
        ),
        de_args
      )
      bundle <- tryCatch(
        invoke_fun(
          method_fun,
          args[names(args) %in% names(formals(method_fun))]
        ),
        error = function(e) {
          list(
            status = "failed",
            reason = e$message,
            results = data.frame(),
            details = list(),
            parameters = list()
          )
        }
      )
      if (!is.list(bundle) || is.null(bundle$status)) {
        bundle <- list(
          status = "failed",
          reason = "Method did not return a valid bundle.",
          results = data.frame(),
          details = list(),
          parameters = list()
        )
      }
      group_bundles[[spec$label]] <- bundle

      if (identical(bundle$status, "success")) {
        df <- as.data.frame(bundle$results)
        if (nrow(df) > 0) {
          if (!"gene" %in% colnames(df)) {
            df$gene <- rownames(df)
          }
          df$group1 <- grp
          df$group2 <- spec$label
          df$comparison <- compose_comparison_label(grp, spec$label)
          df$method <- method_name
          results_list[[paste0(
            grp,
            "::",
            spec$label
          )]] <- coerce_de_schema(df)
        }
      } else {
        fail_messages <- c(
          fail_messages,
          paste0(grp, "[", spec$label, "]: ", bundle$reason %||% "failed")
        )
      }
    }

    per_group[[grp]] <- list(
      markers_type = de_markers_type,
      comparisons = group_bundles
    )
  }

  results <- if (length(results_list) == 0) {
    data.frame(
      gene = character(),
      group1 = character(),
      group2 = character(),
      avg_log2FC = numeric(),
      p_val = numeric(),
      p_val_adj = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, results_list)
    rownames(out) <- NULL
    out
  }

  if (nrow(results) > 0) {
    list(
      status = "success",
      reason = NULL,
      results = results,
      details = list(groups = per_group),
      parameters = list(
        method = method_name,
        de_markers_type = de_markers_type
      )
    )
  } else {
    list(
      status = "failed",
      reason = paste(fail_messages, collapse = "; "),
      details = list(groups = per_group),
      parameters = list(
        method = method_name,
        de_markers_type = de_markers_type
      ),
      results = data.frame(
        gene = character(),
        group1 = character(),
        group2 = character(),
        avg_log2FC = numeric(),
        p_val = numeric(),
        p_val_adj = numeric(),
        method = character(),
        stringsAsFactors = FALSE
      )
    )
  }
}


build_context <- function(
  mode = c("pure_bulk"),
  bulk_se,
  condition.by = NULL,
  group.by = NULL,
  bulk_assay = "counts"
) {
  mode <- match.arg(mode)
  if (!methods::is(bulk_se, "SummarizedExperiment")) {
    log_message(
      "{.arg bulk_se} must be a {.cls SummarizedExperiment}",
      message_type = "error"
    )
  }
  assay_names <- SummarizedExperiment::assayNames(bulk_se)
  if (!bulk_assay %in% assay_names) {
    log_message(
      "{.arg bulk_assay} must be one of {.val {assay_names}}",
      message_type = "error"
    )
  }

  counts <- SummarizedExperiment::assay(bulk_se, bulk_assay)
  if (!inherits(counts, c("matrix", "Matrix"))) {
    counts <- as.matrix(counts)
  }
  sample_meta <- as.data.frame(SummarizedExperiment::colData(bulk_se))
  if (is.null(rownames(sample_meta)) || any(!nzchar(rownames(sample_meta)))) {
    rownames(sample_meta) <- colnames(counts)
  }

  if (is.null(group.by)) {
    group_factor <- factor(rep("all", ncol(counts)), levels = "all")
  } else {
    if (!group.by %in% colnames(sample_meta)) {
      log_message(
        "{.arg group.by} must be present in {.arg colData(bulk_se)}",
        message_type = "error"
      )
    }
    group_factor <- factor(sample_meta[[group.by]])
  }

  condition <- NULL
  if (!is.null(condition.by)) {
    if (!condition.by %in% colnames(sample_meta)) {
      log_message(
        "{.arg condition.by} must be present in {.arg colData(bulk_se)}",
        message_type = "error"
      )
    }
    condition <- factor(sample_meta[[condition.by]])
  }

  group_levels <- levels(group_factor)
  counts_by_group <- stats::setNames(vector("list", length(group_levels)), group_levels)
  condition_by_group <- stats::setNames(vector("list", length(group_levels)), group_levels)
  sample_meta_by_group <- stats::setNames(vector("list", length(group_levels)), group_levels)

  for (grp in group_levels) {
    keep <- !is.na(group_factor) & group_factor == grp
    counts_by_group[[grp]] <- counts[, keep, drop = FALSE]
    sample_meta_by_group[[grp]] <- sample_meta[keep, , drop = FALSE]
    condition_by_group[[grp]] <- if (is.null(condition)) NULL else condition[keep]
  }

  list(
    mode = mode,
    counts_global = counts,
    sample_meta = sample_meta,
    condition_global = condition,
    group_global = group_factor,
    counts_by_group = counts_by_group,
    condition_by_group = condition_by_group,
    sample_meta_by_group = sample_meta_by_group
  )
}


resolve_condition_pair <- function(
  condition,
  condition1 = NULL,
  condition2 = NULL,
  strict_two_levels = FALSE
) {
  condition <- factor(condition)
  levels_use <- levels(droplevels(condition))
  if (length(levels_use) < 2L) {
    stop("At least two condition levels are required.", call. = FALSE)
  }
  if (isTRUE(strict_two_levels) && length(levels_use) != 2L && (is.null(condition1) || is.null(condition2))) {
    stop("Exactly two condition levels are required unless condition1 and condition2 are provided.", call. = FALSE)
  }
  condition1 <- condition1 %||% levels_use[[1]]
  condition2 <- condition2 %||% levels_use[[2]]
  if (!all(c(condition1, condition2) %in% levels_use)) {
    stop("condition1 and condition2 must be present in condition.", call. = FALSE)
  }
  if (identical(condition1, condition2)) {
    stop("condition1 and condition2 must be different.", call. = FALSE)
  }
  list(condition1 = condition1, condition2 = condition2)
}


compose_comparison_label <- function(group, comparison) {
  if (identical(group, "all")) {
    return(as.character(comparison))
  }
  paste(group, comparison, sep = "_")
}


prepare_de_comparisons <- function(
  condition,
  condition1 = NULL,
  condition2 = NULL,
  markers_type = "single"
) {
  markers_type <- normalize_de_markers_type(markers_type)
  cond <- as.character(condition)
  cond_use <- cond[!is.na(cond)]
  levels_use <- unique(cond_use)

  if (identical(markers_type, "single")) {
    pair <- resolve_condition_pair(
      condition = cond,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    )
    return(list(list(
      condition = cond,
      condition1 = pair$condition1,
      condition2 = pair$condition2,
      label = paste0(pair$condition2, "_vs_", pair$condition1)
    )))
  }

  if (length(levels_use) < 2) {
    log_message(
      "At least two condition levels are required for DE comparison.",
      message_type = "error"
    )
  }

  if (identical(markers_type, "all")) {
    return(lapply(levels_use, function(lv) {
      cond_binary <- ifelse(!is.na(cond) & cond == lv, lv, "others")
      list(
        condition = cond_binary,
        condition1 = "others",
        condition2 = lv,
        label = paste0(lv, "_vs_others")
      )
    }))
  }

  pairs <- utils::combn(levels_use, 2, simplify = FALSE)
  lapply(pairs, function(pr) {
    list(
      condition = cond,
      condition1 = pr[[1]],
      condition2 = pr[[2]],
      label = paste0(pr[[2]], "_vs_", pr[[1]])
    )
  })
}


normalize_de_markers_type <- function(markers_type) {
  if (!is.character(markers_type) || length(markers_type) != 1) {
    log_message(
      "{.arg de_markers_type} must be a single string.",
      message_type = "error"
    )
  }
  key <- gsub("[^a-z]", "", tolower(markers_type))
  map <- c(
    single = "single",
    all = "all",
    onevsrest = "all",
    onevsall = "all",
    paired = "paired",
    pairwise = "paired"
  )
  out <- unname(map[[key]])
  if (is.null(out) || is.na(out)) {
    log_message(
      "Unsupported {.arg de_markers_type}. Use one of {.val single}, {.val all}, or {.val paired}.",
      message_type = "error"
    )
  }
  out
}


coerce_de_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      gene = character(),
      group1 = character(),
      group2 = character(),
      avg_log2FC = numeric(),
      p_val = numeric(),
      p_val_adj = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    ))
  }
  if (!"gene" %in% colnames(df)) {
    df$gene <- rownames(df)
  }
  required <- c(
    "group1",
    "group2",
    "avg_log2FC",
    "p_val",
    "p_val_adj",
    "method"
  )
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    for (nm in missing) {
      if (nm %in% c("avg_log2FC", "p_val", "p_val_adj")) {
        df[[nm]] <- NA_real_
      } else {
        df[[nm]] <- NA_character_
      }
    }
  }
  keep <- c(
    "gene",
    "group1",
    "group2",
    "avg_log2FC",
    "p_val",
    "p_val_adj",
    "method"
  )
  optional <- intersect(
    c("pct.1", "pct.2", "ave_expr", "comparison"),
    colnames(df)
  )
  out <- df[, c(keep, optional), drop = FALSE]
  rownames(out) <- NULL
  out
}

filter_de_results <- function(de_results, DE_threshold) {
  de_results <- detest_res(de_results)
  if (nrow(de_results) == 0) {
    return(de_results)
  }

  if ("pct.1" %in% colnames(de_results) && "pct.2" %in% colnames(de_results)) {
    de_results[["diff_pct"]] <- de_results[["pct.1"]] - de_results[["pct.2"]]
  }
  if ("p_val_adj" %in% colnames(de_results)) {
    de_results[["minus_log10padj"]] <- -log10(de_results[["p_val_adj"]])
    de_results[["-log10padj"]] <- de_results[["minus_log10padj"]]
  }

  keep <- tryCatch(
    with(de_results, eval(rlang::parse_expr(DE_threshold))),
    error = function(e) {
      log_message(
        paste0(
          "Failed to evaluate {.arg DE_threshold}: ",
          conditionMessage(e)
        ),
        message_type = "error"
      )
    }
  )
  keep <- keep %in% TRUE
  de_results[keep, , drop = FALSE]
}

prepare_de_for_pathway <- function(de_results, require_score = FALSE) {
  de_results <- detest_res(de_results)
  if (nrow(de_results) == 0) {
    return(de_results[0, , drop = FALSE])
  }
  if (!"gene" %in% colnames(de_results)) {
    log_message(
      "Cannot find {.field gene} in DEtest results.",
      message_type = "error"
    )
  }
  if (isTRUE(require_score) && !"avg_log2FC" %in% colnames(de_results)) {
    log_message(
      "Cannot find {.field avg_log2FC} in DEtest results for GSEA.",
      message_type = "error"
    )
  }

  comparison <- rep(NA_character_, nrow(de_results))
  if ("comparison" %in% colnames(de_results)) {
    comparison <- as.character(de_results[["comparison"]])
  }
  missing_comparison <- is.na(comparison) | !nzchar(comparison)
  if (any(missing_comparison) && "group1" %in% colnames(de_results)) {
    comparison[missing_comparison] <- as.character(de_results[["group1"]])[
      missing_comparison
    ]
  }
  missing_comparison <- is.na(comparison) | !nzchar(comparison)
  if (any(missing_comparison) && "cluster" %in% colnames(de_results)) {
    comparison[missing_comparison] <- as.character(de_results[["cluster"]])[
      missing_comparison
    ]
  }
  missing_comparison <- is.na(comparison) | !nzchar(comparison)
  if (any(missing_comparison)) {
    comparison[missing_comparison] <- "All"
  }

  de_results[["comparison"]] <- comparison
  de_results <- de_results[
    !is.na(de_results[["gene"]]) & nzchar(as.character(de_results[["gene"]])),
    ,
    drop = FALSE
  ]
  if (isTRUE(require_score)) {
    de_results <- de_results[
      !is.na(de_results[["avg_log2FC"]]),
      ,
      drop = FALSE
    ]
  }
  de_results
}

#' Calculate expression fold change
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return A data frame of fold-change statistics.
#' @export
FoldChange <- function(object, ...) {
  UseMethod("FoldChange")
}

#' @export
FoldChange.default <- function(
  object,
  cells.1,
  cells.2,
  mean.fxn,
  fc.name,
  features = NULL,
  ...
) {
  features <- features %||% rownames(x = object)
  thresh.min <- 0
  pct.1 <- round(
    x = Matrix::rowSums(
      x = object[features, cells.1, drop = FALSE] > thresh.min,
      na.rm = TRUE
    ) /
      length(cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = Matrix::rowSums(
      x = object[features, cells.2, drop = FALSE] > thresh.min,
      na.rm = TRUE
    ) /
      length(cells.2),
    digits = 3
  )
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}


PerformDE <- function(
  object,
  cells.1,
  cells.2,
  features,
  test.use,
  verbose,
  min.cells.feature,
  latent.vars,
  min.expression = 0,
  ...
) {
  if (
    !(test.use %in% c("negbinom", "poisson", "MAST", "LR")) &&
      !is.null(latent.vars)
  ) {
    log_message(
      "'latent.vars' is only used for the following tests: ",
      paste(c("negbinom", "poisson", "MAST", "LR"), collapse = ", "),
      message_type = "warning",
      verbose = verbose
    )
  }
  data.use <- object[
    features,
    c(cells.1, cells.2),
    drop = FALSE
  ]
  dense_data_use <- function() {
    data.use_dense <- as_matrix(data.use)
    data.use_dense[data.use_dense <= min.expression] <- NA
    data.use_dense
  }

  de.results <- switch(
    EXPR = test.use,
    "wilcox" = WilcoxDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.expression = min.expression,
      verbose = verbose,
      ...
    ),
    "bimod" = get_namespace_fun(
      "Seurat",
      "DiffExpTest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "roc" = get_namespace_fun(
      "Seurat",
      "MarkerTest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "t" = get_namespace_fun(
      "Seurat",
      "DiffTTest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "negbinom" = get_namespace_fun(
      "Seurat",
      "GLMDETest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "poisson" = get_namespace_fun(
      "Seurat",
      "GLMDETest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "MAST" = get_namespace_fun(
      "Seurat",
      "MASTDETest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = get_namespace_fun(
      "Seurat",
      "DESeq2DETest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = get_namespace_fun(
      "Seurat",
      "LRDETest"
    )(
      data.use = dense_data_use(),
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    log_message(
      "Unknown test: {.pkg {test.use}}",
      message_type = "error"
    )
  )
  de.results
}

WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  min.expression = 0,
  verbose = TRUE,
  ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  if (
    inherits(data.use, "sparseMatrix") &&
      isTRUE(min.expression >= 0)
  ) {
    p_val <- run_sparse_wilcox(
      x = data.use,
      n_group1 = length(cells.1),
      min.expression = min.expression
    )
    return(data.frame(
      p_val = p_val,
      row.names = names(p_val)
    ))
  }

  data.use <- as_matrix(data.use)
  data.use[data.use <= min.expression] <- NA
  check_r("limma", verbose = FALSE)
  p_val <- parallelize_fun(
    seq_len(nrow(data.use)),
    fun = function(x) {
      keep <- colnames(data.use)[!is.na(data.use[x, ])]
      j <- seq_len(length.out = length(intersect(cells.1, keep)))
      statistics <- data.use[x, keep]
      min(
        2 *
          min(
            limma::rankSumTestWithCorrelation(
              index = j,
              statistics = statistics
            )
          ),
        1
      )
    }
  ) |>
    purrr::list_c()
  data.frame(
    p_val,
    row.names = rownames(data.use)
  )
}

RunDEtestFindMarkers <- function(
  srt,
  assay,
  layer,
  cells.1,
  cells.2,
  features,
  test.use,
  logfc.threshold,
  base,
  min.pct,
  min.diff.pct,
  max.cells.per.ident,
  min.cells.feature,
  min.cells.group,
  latent.vars,
  only.pos,
  norm.method,
  pseudocount.use,
  mean.fxn,
  verbose,
  random.seed = 1,
  ...
) {
  extra_args <- list(...)
  use_sparse_wilcox <- identical(test.use, "wilcox") &&
    is.null(latent.vars) &&
    is.null(mean.fxn) &&
    layer %in% c("data", "counts") &&
    identical(norm.method, "LogNormalize") &&
    length(extra_args) == 0 &&
    !requireNamespace(paste0("pre", "sto"), quietly = TRUE)
  if (isTRUE(use_sparse_wilcox)) {
    return(FindMarkers(
      object = srt,
      assay = assay,
      layer = layer,
      cells.1 = cells.1,
      cells.2 = cells.2,
      features = features,
      test.use = test.use,
      logfc.threshold = logfc.threshold,
      base = base,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      max.cells.per.ident = max.cells.per.ident,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      only.pos = only.pos,
      norm.method = norm.method,
      pseudocount.use = pseudocount.use,
      random.seed = random.seed,
      verbose = verbose
    ))
  }

  do.call(
    Seurat::FindMarkers,
    c(
      list(
        object = Seurat::GetAssay(srt, assay),
        layer = layer,
        cells.1 = cells.1,
        cells.2 = cells.2,
        features = features,
        test.use = test.use,
        logfc.threshold = logfc.threshold,
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = verbose,
        random.seed = random.seed
      ),
      extra_args
    )
  )
}

RunDEtestSparseWilcoxMarkers <- function(
  srt,
  assay,
  layer,
  cells.1,
  cells.2,
  features,
  logfc.threshold,
  base,
  min.pct,
  min.diff.pct,
  max.cells.per.ident,
  min.cells.group,
  only.pos,
  pseudocount.use,
  random.seed,
  verbose
) {
  data.use <- GetAssayData5(
    object = srt,
    assay = assay,
    layer = layer
  )
  features <- features %||% rownames(data.use)
  missing_features <- setdiff(features, rownames(data.use))
  if (length(missing_features) > 0) {
    log_message(
      "{.arg features} contains features not present in the selected assay/layer",
      message_type = "error"
    )
  }

  ValidateCellGroups <- get_namespace_fun("Seurat", "ValidateCellGroups")
  ValidateCellGroups(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )

  fc.results <- RunDEtestSparseFoldChange(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    layer = layer,
    pseudocount.use = pseudocount.use,
    base = base
  )

  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(alpha.min) <- rownames(fc.results)
  features <- names(alpha.min)[alpha.min >= min.pct]
  if (length(features) == 0) {
    return(fc.results[features, , drop = FALSE])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(alpha.min)[
    alpha.min >= min.pct & alpha.diff >= min.diff.pct
  ]
  if (length(features) == 0) {
    return(fc.results[features, , drop = FALSE])
  }

  fc_col <- colnames(fc.results)[1]
  total.diff <- fc.results[[fc_col]]
  names(total.diff) <- rownames(fc.results)
  features.diff <- if (isTRUE(only.pos)) {
    names(total.diff)[total.diff >= logfc.threshold]
  } else {
    names(total.diff)[abs(total.diff) >= logfc.threshold]
  }
  features <- intersect(features, features.diff)
  if (length(features) == 0) {
    return(fc.results[features, , drop = FALSE])
  }

  if (max.cells.per.ident < Inf) {
    if (!is.null(random.seed)) {
      set.seed(random.seed)
    }
    if (length(cells.1) > max.cells.per.ident) {
      cells.1 <- sample(cells.1, size = max.cells.per.ident)
    }
    if (length(cells.2) > max.cells.per.ident) {
      cells.2 <- sample(cells.2, size = max.cells.per.ident)
    }
  }

  data.de <- data.use[features, c(cells.1, cells.2), drop = FALSE]
  p_val <- run_sparse_wilcox_all_cells(
    x = data.de,
    n_group1 = length(cells.1)
  )
  de.results <- data.frame(
    p_val = p_val,
    row.names = names(p_val)
  )
  de.results <- cbind(
    de.results,
    fc.results[rownames(de.results), , drop = FALSE]
  )
  de.results <- de.results[order(de.results$p_val, -de.results[[fc_col]]), ]
  de.results$p_val_adj <- stats::p.adjust(
    p = de.results$p_val,
    method = "bonferroni",
    n = nrow(data.use)
  )
  de.results
}

RunDEtestSparseFoldChange <- function(
  object,
  cells.1,
  cells.2,
  features,
  layer,
  pseudocount.use,
  base
) {
  base.text <- ifelse(base == exp(1), "", base)
  fc.name <- ifelse(
    layer == "scale.data",
    "avg_diff",
    paste0("avg_log", base.text, "FC")
  )
  group1 <- RunDEtestSparseGroupStats(
    object = object,
    features = features,
    cells = cells.1,
    layer = layer,
    pseudocount.use = pseudocount.use,
    base = base
  )
  group2 <- RunDEtestSparseGroupStats(
    object = object,
    features = features,
    cells = cells.2,
    layer = layer,
    pseudocount.use = pseudocount.use,
    base = base
  )
  fc.results <- as.data.frame(
    cbind(
      group1[["mean"]] - group2[["mean"]],
      group1[["pct"]],
      group2[["pct"]]
    )
  )
  rownames(fc.results) <- features
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  fc.results
}

RunDEtestSparseGroupStats <- function(
  object,
  features,
  cells,
  layer,
  pseudocount.use,
  base
) {
  mat <- object[features, cells, drop = FALSE]
  mat <- methods::as(mat, "dgCMatrix")
  detected <- Matrix::rowSums(mat > 0)
  pct <- round(detected / length(cells), digits = 3)

  if (identical(layer, "data")) {
    mat_sum <- mat
    mat_sum@x <- expm1(mat_sum@x)
    sum_vals <- Matrix::rowSums(mat_sum)
    mean_vals <- log((sum_vals + pseudocount.use) / length(cells), base = base)
  } else {
    sum_vals <- Matrix::rowSums(mat)
    mean_vals <- log((sum_vals + pseudocount.use) / length(cells), base = base)
  }
  names(mean_vals) <- features
  names(pct) <- features
  list(mean = mean_vals, pct = pct)
}

aggregate_counts_by_group <- function(counts, groups) {
  groups <- as.character(groups)
  groups <- groups[!is.na(groups)]
  group_levels <- unique(groups)
  agg <- lapply(group_levels, function(group) {
    Matrix::rowSums(counts[, groups == group, drop = FALSE])
  })
  agg <- do.call(cbind, agg)
  if (is.null(dim(agg))) {
    agg <- matrix(agg, ncol = 1)
  }
  rownames(agg) <- rownames(counts)
  colnames(agg) <- group_levels
  agg
}

RunDEtest_limma <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = TRUE,
  logfc.threshold = 0,
  p.adjust.method = "bonferroni"
) {
  check_r("limma", verbose = FALSE)
  check_r("edgeR", verbose = FALSE)

  condition_all <- as.character(condition)
  if (is.null(condition1) || is.null(condition2)) {
    condition_levels <- unique(condition_all)
    if (length(condition_levels) < 2) {
      return(NULL)
    }
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
  }

  keep <- condition_all %in% c(condition1, condition2)
  count_matrix <- count_matrix[, keep, drop = FALSE]
  condition <- factor(condition_all[keep], levels = c(condition1, condition2))
  if (ncol(count_matrix) < 2 || any(table(condition) < 2)) {
    return(NULL)
  }

  dge <- edgeR::DGEList(counts = count_matrix)
  keep_features <- edgeR::filterByExpr(dge, group = condition)
  if (!any(keep_features)) {
    return(data.frame())
  }
  dge <- dge[keep_features, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~condition)
  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(
    fit,
    coef = ncol(design),
    number = Inf,
    sort.by = "none"
  )
  if (is.null(tt) || nrow(tt) == 0) {
    return(data.frame())
  }

  detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
  pct.1 <- round(rowMeans(detect[, condition == condition1, drop = FALSE]), 3)
  pct.2 <- round(rowMeans(detect[, condition == condition2, drop = FALSE]), 3)
  out <- data.frame(
    p_val = tt$P.Value,
    avg_log2FC = tt$logFC,
    ave_expr = tt$AveExpr,
    pct.1 = pct.1,
    pct.2 = pct.2,
    row.names = rownames(tt)
  )
  out$p_val_adj <- stats::p.adjust(out$p_val, method = p.adjust.method)
  if (isTRUE(only.pos)) {
    out <- out[out$avg_log2FC >= logfc.threshold, , drop = FALSE]
  } else {
    out <- out[abs(out$avg_log2FC) >= logfc.threshold, , drop = FALSE]
  }
  out
}

RunDEtest_edgeR <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = TRUE,
  logfc.threshold = 0,
  p.adjust.method = "bonferroni"
) {
  check_r("edgeR", verbose = FALSE)

  condition_all <- as.character(condition)
  if (is.null(condition1) || is.null(condition2)) {
    condition_levels <- unique(condition_all)
    if (length(condition_levels) < 2) {
      return(NULL)
    }
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
  }

  keep <- condition_all %in% c(condition1, condition2)
  count_matrix <- count_matrix[, keep, drop = FALSE]
  condition <- factor(condition_all[keep], levels = c(condition1, condition2))
  if (ncol(count_matrix) < 2 || any(table(condition) < 2)) {
    return(NULL)
  }

  dge <- edgeR::DGEList(counts = count_matrix)
  keep_features <- edgeR::filterByExpr(dge, group = condition)
  if (!any(keep_features)) {
    return(data.frame())
  }
  dge <- dge[keep_features, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~condition)
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  if (is.null(tt) || nrow(tt) == 0) {
    return(data.frame())
  }

  detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
  pct.1 <- round(rowMeans(detect[, condition == condition1, drop = FALSE]), 3)
  pct.2 <- round(rowMeans(detect[, condition == condition2, drop = FALSE]), 3)
  out <- data.frame(
    p_val = tt$PValue,
    avg_log2FC = tt$logFC,
    ave_expr = tt$logCPM,
    pct.1 = pct.1,
    pct.2 = pct.2,
    row.names = rownames(tt)
  )
  out$p_val_adj <- stats::p.adjust(out$p_val, method = p.adjust.method)
  if (isTRUE(only.pos)) {
    out <- out[out$avg_log2FC >= logfc.threshold, , drop = FALSE]
  } else {
    out <- out[abs(out$avg_log2FC) >= logfc.threshold, , drop = FALSE]
  }
  out
}

RunDEtest_pseudobulk <- function(
  srt,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  features = NULL,
  feature_type = "gene",
  markers_type = "all",
  test.use = "limma",
  only.pos = TRUE,
  fc.threshold = 1.5,
  base = 2,
  sample_col = NULL,
  condition_col = NULL,
  p.adjust.method = "bonferroni",
  layer = "counts",
  assay = NULL,
  verbose = TRUE,
  cores = 1,
  ...
) {
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(srt)
  }
  if (!markers_type %in% c("all")) {
    log_message(
      "Sample-level differential testing currently supports only {.val all} markers.",
      message_type = "error"
    )
  }
  if (!test.use %in% c("limma", "edgeR")) {
    log_message(
      "Sample-level differential testing currently supports only {.val limma} and {.val edgeR}.",
      message_type = "error"
    )
  }
  if (is.null(sample_col) || !sample_col %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg sample_col} must be a metadata column in {.cls Seurat} for sample-level differential testing",
      message_type = "error"
    )
  }
  if (is.null(condition_col) || !condition_col %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg condition_col} must be a metadata column in {.cls Seurat} for sample-level differential testing",
      message_type = "error"
    )
  }
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be a metadata column in {.cls Seurat} for sample-level differential testing",
      message_type = "error"
    )
  }

  counts <- GetAssayData5(srt, layer = layer, assay = assay)
  features <- features %||% rownames(counts)
  features <- intersect(features, rownames(counts))
  if (length(features) == 0) {
    log_message(
      "No valid features available for sample-level differential testing",
      message_type = "error"
    )
  }
  counts <- counts[features, , drop = FALSE]

  condition_levels_all <- unique(as.character(stats::na.omit(srt[[
    condition_col,
    drop = TRUE
  ]])))
  condition1 <- group1 %||% NULL
  condition2 <- group2 %||% NULL
  if (is.null(condition1) || is.null(condition2)) {
    if (length(condition_levels_all) < 2) {
      log_message(
        "Sample-level differential testing requires at least two condition levels",
        message_type = "error"
      )
    }
    condition1 <- condition1 %||% condition_levels_all[[1]]
    condition2 <- condition2 %||% condition_levels_all[[2]]
  }
  if (identical(condition1, condition2)) {
    log_message(
      "{.arg group1} and {.arg group2} must refer to two different condition labels in sample-level differential testing",
      message_type = "error"
    )
  }

  target_groups <- if (is.null(group.by)) {
    "All"
  } else {
    unique(as.character(stats::na.omit(srt[[group.by, drop = TRUE]])))
  }

  log_message(
    "Start sample-level differential testing",
    verbose = verbose
  )

  run_group_fun <- if (isTRUE(cores == 1)) {
    function(X, FUN) lapply(X, FUN)
  } else {
    function(X, FUN) parallelize_fun(X, FUN, cores = cores, verbose = verbose)
  }

  all_markers <- run_group_fun(
    target_groups,
    function(current_group) {
      cells_use <- colnames(srt)
      if (!is.null(group.by)) {
        cells_use <- colnames(srt)[
          srt[[group.by, drop = TRUE]] %in% current_group
        ]
      }
      meta_use <- srt@meta.data[
        cells_use,
        c(sample_col, condition_col),
        drop = FALSE
      ]
      meta_use <- meta_use[stats::complete.cases(meta_use), , drop = FALSE]
      cells_use <- rownames(meta_use)
      if (length(cells_use) == 0) {
        return(NULL)
      }

      sample_condition <- unique(meta_use[,
        c(sample_col, condition_col),
        drop = FALSE
      ])
      sample_tab <- table(meta_use[[sample_col]], meta_use[[condition_col]])
      if (any(rowSums(sample_tab > 0) > 1)) {
        log_message(
          "Each sample must map to a single condition in sample-level differential testing",
          message_type = "error"
        )
      }

      counts_use <- counts[, cells_use, drop = FALSE]
      count_matrix <- aggregate_counts_by_group(
        counts = counts_use,
        groups = meta_use[[sample_col]]
      )
      condition_map <- stats::setNames(
        as.character(sample_condition[[condition_col]]),
        sample_condition[[sample_col]]
      )
      condition_use <- condition_map[colnames(count_matrix)]
      markers <- switch(test.use,
        limma = RunDEtest_limma(
          count_matrix = count_matrix,
          condition = condition_use,
          condition1 = condition1,
          condition2 = condition2,
          only.pos = only.pos,
          logfc.threshold = log(fc.threshold, base = base),
          p.adjust.method = p.adjust.method
        ),
        edgeR = RunDEtest_edgeR(
          count_matrix = count_matrix,
          condition = condition_use,
          condition1 = condition1,
          condition2 = condition2,
          only.pos = only.pos,
          logfc.threshold = log(fc.threshold, base = base),
          p.adjust.method = p.adjust.method
        ),
        DESeq2 = RunDEtest_DESeq2(
          count_matrix = count_matrix,
          condition = condition_use,
          condition1 = condition1,
          condition2 = condition2,
          only.pos = only.pos,
          logfc.threshold = log(fc.threshold, base = base),
          p.adjust.method = p.adjust.method
        ),
        dream = RunDEtest_dream(
          count_matrix = count_matrix,
          condition = condition_use,
          condition1 = condition1,
          condition2 = condition2,
          only.pos = only.pos,
          logfc.threshold = log(fc.threshold, base = base),
          p.adjust.method = p.adjust.method
        )
      )
      if (is.null(markers) || nrow(markers) == 0) {
        return(NULL)
      }
      markers[, "gene"] <- rownames(markers)
      markers[, "feature"] <- rownames(markers)
      markers[, "feature_type"] <- feature_type
      markers[, "group1"] <- current_group
      markers[, "group2"] <- paste0(condition2, "_vs_", condition1)
      markers[, "condition1"] <- condition1
      markers[, "condition2"] <- condition2
      markers[, "sample_number1"] <- sum(
        condition_use == condition1,
        na.rm = TRUE
      )
      markers[, "sample_number2"] <- sum(
        condition_use == condition2,
        na.rm = TRUE
      )
      markers
    }
  )
  all_markers <- do.call(rbind.data.frame, all_markers)
  tool_name <- if (is.null(group.by)) {
    "DEtest_pseudobulk"
  } else {
    paste0("DEtest_", group.by)
  }
  if (is.null(srt@tools[[tool_name]])) {
    srt@tools[[tool_name]] <- list()
  }
  srt@tools[[tool_name]][["feature_type"]] <- feature_type
  srt@tools[[tool_name]][["sample_col"]] <- sample_col
  srt@tools[[tool_name]][["condition_col"]] <- condition_col
  srt@tools[[tool_name]][["condition1"]] <- condition1
  srt@tools[[tool_name]][["condition2"]] <- condition2
  srt@tools[[tool_name]][["test.use"]] <- test.use
  srt@tools[[tool_name]][["assay"]] <- assay
  srt@tools[[tool_name]][["layer"]] <- layer
  srt@tools[[tool_name]][["group.by"]] <- group.by
  srt@tools[[tool_name]][["groups_tested"]] <- target_groups

  if (is.null(all_markers) || nrow(all_markers) == 0) {
    srt@tools[[tool_name]][[paste0("AllMarkers_", test.use)]] <- data.frame()
  } else {
    rownames(all_markers) <- NULL
    all_markers[, "group1"] <- factor(
      all_markers[, "group1"],
      levels = target_groups
    )
    all_markers[, "test_group_number"] <- as.integer(
      table(all_markers[["gene"]])[all_markers[, "gene"]]
    )
    all_markers_matrix <- as.data.frame.matrix(
      table(all_markers[, c("gene", "group1")])
    )
    all_markers[, "test_group"] <- apply(all_markers_matrix, 1, function(x) {
      paste0(colnames(all_markers_matrix)[x > 0], collapse = ";")
    })[all_markers[, "gene"]]
    srt@tools[[tool_name]][[paste0("AllMarkers_", test.use)]] <- all_markers
  }

  log_message(
    "Sample-level differential testing completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Differential gene test
#'
#' @description
#' Perform differential expression testing on a `Seurat` object or a
#' `SummarizedExperiment` bulk object.
#' Users have the flexibility to specify custom cell groups, marker types, and
#' various options for DE analysis.
#'
#' @md
#' @inheritParams thisutils::parallelize_fun
#' @inheritParams Seurat::FindMarkers
#' @inheritParams standard_scop
#' @inheritParams FeatureDimPlot
#' @param object A `Seurat` object or a `SummarizedExperiment` object.
#' @param srt Compatibility alias for `object`.
#' @param group.by A grouping variable in the dataset to define the groups or conditions for the differential test.
#' If not provided, the function uses the "active.ident" variable in the Seurat object.
#' @param group1 A vector of cell IDs or a character vector specifying the cells that belong to the first group.
#' If both group.by and group1 are provided, group1 takes precedence.
#' For sample-level methods (`"edgeR"`, `"limma"`, `"DESeq2"`, and `"dream"`),
#' this parameter is interpreted as the first condition label.
#' @param group2 A vector of cell IDs or a character vector specifying the cells that belong to the second group.
#' This parameter is only used when group.by or group1 is provided.
#' For sample-level methods (`"edgeR"`, `"limma"`, `"DESeq2"`, and `"dream"`),
#' this parameter is interpreted as the second condition label.
#' @param ident.1,ident.2 Seurat-style aliases for `group1` and `group2`.
#' @param cells1 A vector of cell IDs specifying the cells that belong to group1. If provided, group1 is ignored.
#' @param cells2 A vector of cell IDs specifying the cells that belong to group2.
#' This parameter is only used when cells1 is provided.
#' @param cells.1,cells.2 Seurat-style aliases for `cells1` and `cells2`.
#' @param features A vector of feature names specifying the features to consider for the differential test.
#' If not provided, all features in the dataset are considered.
#' @param feature_type Feature type used for differential testing.
#' Default is `"gene"`.
#' @param markers_type A character value specifying the type of markers to find.
#' Possible values are "all", "paired", "conserved", and "disturbed".
#' Sample-level methods (`"edgeR"`, `"limma"`, `"DESeq2"`, and `"dream"`)
#' currently support only `"all"`.
#' @param grouping.var A character value specifying the grouping variable for finding conserved or disturbed markers.
#' This parameter is only used when markers_type is "conserved" or "disturbed".
#' @param fc.threshold A numeric value used to filter genes for testing based on their average fold change between/among the two groups.
#' Default is `1.5`.
#' @param logfc.threshold Seurat-style log fold-change threshold. When provided,
#' it is converted to `fc.threshold = base^logfc.threshold`.
#' @param meta.method A character value specifying the method to use for combining p-values in the conserved markers test.
#' Possible values are "maximump", "minimump", "wilkinsonp", "meanp", "sump", and "votep".
#' @param norm.method Normalization method for fold change calculation when layer is 'data'.
#' Default is `"LogNormalize"`.
#' @param sample_col Metadata column storing biological sample IDs.
#' Required when `test.use` is a sample-level pseudobulk method on `Seurat`.
#' @param condition_col Metadata column storing condition labels.
#' Required when `test.use` is a sample-level pseudobulk method on `Seurat`, and
#' required for `SummarizedExperiment` input.
#' @param bulk_assay Assay name used as the bulk counts matrix for
#' `SummarizedExperiment` input.
#' @param p.adjust.method A character value specifying the method to use for adjusting p-values.
#' Default is `"bonferroni"`.
#' @param test.use Differential testing method.
#' `"edgeR"`, `"limma"`, `"DESeq2"`, and `"dream"` run sample-level
#' pseudobulk differential testing on `Seurat` input and bulk DE on
#' `SummarizedExperiment` input.
#' @param ... Additional arguments to pass to the [Seurat::FindMarkers] function.
#'
#' @export
#'
#' @seealso
#' [VolcanoPlot], [RunEnrichment], [RunGSEA], [GroupHeatmap]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   only.pos = FALSE
#' )
#' AllMarkers <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_SubCellType$AllMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht1 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = AllMarkers$gene,
#'   feature_split = AllMarkers$group1,
#'   group.by = "SubCellType"
#' )
#' ht1$plot
#'
#' TopMarkers <- AllMarkers |>
#'   dplyr::group_by(gene) |>
#'   dplyr::top_n(1, avg_log2FC) |>
#'   dplyr::group_by(group1) |>
#'   dplyr::top_n(3, avg_log2FC)
#' ht2 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = TopMarkers$gene,
#'   feature_split = TopMarkers$group1,
#'   group.by = "SubCellType",
#'   show_row_names = TRUE
#' )
#' ht2$plot
#'
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   markers_type = "paired",
#'   cores = 2
#' )
#' PairedMarkers <- dplyr::filter(
#'   pancreas_sub@tools$DEtest_SubCellType$PairedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht3 <- GroupHeatmap(
#'   pancreas_sub,
#'   features = PairedMarkers$gene,
#'   feature_split = PairedMarkers$group1,
#'   group.by = "SubCellType"
#' )
#' ht3$plot
#'
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Uncorrected"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("celltype", "tech")
#' )
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group.by = "celltype",
#'   grouping.var = "tech",
#'   markers_type = "conserved",
#'   cores = 2
#' )
#' ConservedMarkers1 <- dplyr::filter(
#'   panc8_sub@tools$DEtest_celltype$ConservedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht4 <- GroupHeatmap(
#'   panc8_sub,
#'   layer = "data",
#'   features = ConservedMarkers1$gene,
#'   feature_split = ConservedMarkers1$group1,
#'   group.by = "tech",
#'   split.by = "celltype",
#'   within_groups = TRUE
#' )
#' ht4$plot
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group.by = "tech",
#'   grouping.var = "celltype",
#'   markers_type = "conserved",
#'   cores = 2
#' )
#' ConservedMarkers2 <- dplyr::filter(
#'   panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1
#' )
#' ht4 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   features = ConservedMarkers2$gene,
#'   feature_split = ConservedMarkers2$group1,
#'   group.by = "tech",
#'   split.by = "celltype"
#' )
#' ht4$plot
#'
#' panc8_sub <- RunDEtest(
#'   srt = panc8_sub,
#'   group.by = "celltype",
#'   grouping.var = "tech",
#'   markers_type = "disturbed",
#'   cores = 2
#' )
#' DisturbedMarkers <- dplyr::filter(
#'   panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
#'   p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
#' )
#' ht5 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   features = DisturbedMarkers$gene,
#'   feature_split = DisturbedMarkers$group1,
#'   group.by = "celltype",
#'   split.by = "tech"
#' )
#' ht5$plot
#'
#' gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
#' DisturbedMarkers_specific <- DisturbedMarkers[
#'   DisturbedMarkers$gene %in% gene_specific,
#' ]
#' ht6 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   features = DisturbedMarkers_specific$gene,
#'   feature_split = DisturbedMarkers_specific$group1,
#'   group.by = "celltype",
#'   split.by = "tech"
#' )
#' ht6$plot
#'
#' ht7 <- GroupHeatmap(
#'   srt = panc8_sub,
#'   layer = "data",
#'   aggregate_fun = function(x) mean(expm1(x)) + 1,
#'   features = DisturbedMarkers_specific$gene,
#'   feature_split = DisturbedMarkers_specific$group1,
#'   group.by = "celltype",
#'   grouping.var = "tech",
#'   numerator = "smartseq2"
#' )
#' ht7$plot
#'
#' cell_index <- ave(
#'   seq_along(pancreas_sub$CellType),
#'   pancreas_sub$CellType,
#'   FUN = seq_along
#' )
#' pancreas_sub[["sample"]] <- paste0(
#'   "S",
#'   (cell_index - 1) %% 4 + 1
#' )
#' pancreas_sub[["condition"]] <- ifelse(
#'   pancreas_sub$sample %in% c("S1", "S2"),
#'   "ctrl",
#'   "case"
#' )
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   test.use = "limma",
#'   fc.threshold = 1,
#'   layer = "counts",
#'   only.pos = FALSE
#' )
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   test.use = "limma",
#'   group_use = "Ductal",
#'   plot_type = "volcano",
#'   x_metric = "avg_log2FC",
#'   y_metric = "p_val"
#' )
#'
#' pancreas_sub <- RunDEtest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   test.use = "edgeR",
#'   layer = "counts",
#'   fc.threshold = 1,
#'   only.pos = FALSE
#' )
#' DEtestPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   test.use = "edgeR",
#'   group_use = "Ductal",
#'   plot_type = "volcano",
#'   x_metric = "avg_log2FC",
#'   y_metric = "p_val",
#'   DE_threshold = "abs(avg_log2FC) > log2(1.5) & p_val < 0.05"
#' )
#'
#' data(islet_bulk)
#' bulk_out <- RunDEtest(
#'   islet_bulk,
#'   condition_col = "condition",
#'   group1 = "control",
#'   group2 = "bfa",
#'   test.use = "edgeR",
#'   only.pos = FALSE,
#'   fc.threshold = 1
#' )
#' DEtestPlot(
#'   bulk_out,
#'   test.use = "edgeR",
#'   plot_type = "volcano",
#'   x_metric = "avg_log2FC",
#'   y_metric = "p_val"
#' )
RunDEtest <- function(
  object = NULL,
  ...,
  srt = NULL
) {
  if (is.null(object)) {
    object <- srt
    if (methods::is(object, "SummarizedExperiment")) {
      return(RunDEtest.SummarizedExperiment(object, ...))
    }
    if (methods::is(object, "Seurat")) {
      return(RunDEtest.Seurat(object, ...))
    }
  }
  if (methods::is(object, "SummarizedExperiment")) {
    return(RunDEtest.SummarizedExperiment(object, ...))
  }
  UseMethod(generic = "RunDEtest", object = object)
}

#' @rdname RunDEtest
#' @export
RunDEtest.Seurat <- function(
  object,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  ident.1 = NULL,
  ident.2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  feature_type = c("gene", "peak", "cCRE"),
  markers_type = c(
    "all",
    "paired",
    "conserved",
    "disturbed"
  ),
  grouping.var = NULL,
  meta.method = c(
    "maximump",
    "minimump",
    "wilkinsonp",
    "meanp",
    "sump",
    "votep"
  ),
  test.use = "wilcox",
  only.pos = TRUE,
  fc.threshold = 1.5,
  logfc.threshold = NULL,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  norm.method = "LogNormalize",
  sample_col = NULL,
  condition_col = NULL,
  bulk_assay = "counts",
  p.adjust.method = "bonferroni",
  layer = "data",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
) {
  srt <- object
  set.seed(seed)
  feature_type <- match.arg(feature_type)
  markers_type <- match.arg(markers_type)
  meta.method <- match.arg(meta.method)
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(srt)
  }
  if (!is.null(ident.1)) {
    if (
      !is.null(group1) &&
        !identical(as.character(group1), as.character(ident.1))
    ) {
      log_message(
        "Both {.arg group1} and Seurat-style {.arg ident.1} were provided; using {.arg ident.1}.",
        message_type = "warning",
        verbose = verbose
      )
    }
    group1 <- ident.1
  }
  if (!is.null(ident.2)) {
    if (
      !is.null(group2) &&
        !identical(as.character(group2), as.character(ident.2))
    ) {
      log_message(
        "Both {.arg group2} and Seurat-style {.arg ident.2} were provided; using {.arg ident.2}.",
        message_type = "warning",
        verbose = verbose
      )
    }
    group2 <- ident.2
  }
  if (!is.null(cells.1)) {
    if (
      !is.null(cells1) &&
        !identical(as.character(cells1), as.character(cells.1))
    ) {
      log_message(
        "Both {.arg cells1} and Seurat-style {.arg cells.1} were provided; using {.arg cells.1}.",
        message_type = "warning",
        verbose = verbose
      )
    }
    cells1 <- cells.1
  }
  if (!is.null(cells.2)) {
    if (
      !is.null(cells2) &&
        !identical(as.character(cells2), as.character(cells.2))
    ) {
      log_message(
        "Both {.arg cells2} and Seurat-style {.arg cells.2} were provided; using {.arg cells.2}.",
        message_type = "warning",
        verbose = verbose
      )
    }
    cells2 <- cells.2
  }
  if (!is.null(logfc.threshold)) {
    if (!missing(fc.threshold)) {
      log_message(
        "Both {.arg fc.threshold} and Seurat-style {.arg logfc.threshold} were provided; using {.arg logfc.threshold}.",
        message_type = "warning",
        verbose = verbose
      )
    }
    fc.threshold <- base^logfc.threshold
  }

  sample_level_methods <- c("edgeR", "limma", "DESeq2", "dream")
  is_sample_level <- test.use %in% sample_level_methods
  if (is_sample_level) {
    if (markers_type != "all") {
      log_message(
        "Sample-level differential testing currently supports only {.val all} markers.",
        message_type = "error"
      )
    }
    if (!is.null(cells1) || !is.null(cells2)) {
      log_message(
        "{.arg cells1} and {.arg cells2} are not supported for sample-level differential testing. Use {.arg group1} and {.arg group2} to specify condition labels.",
        message_type = "error"
      )
    }
    if (!is.null(grouping.var)) {
      log_message(
        "{.arg grouping.var} is not supported for sample-level differential testing.",
        message_type = "error"
      )
    }
    if (!identical(layer, "counts")) {
      log_message(
        "Sample-level differential testing uses the {.arg counts} layer. Reset {.arg layer = 'counts'}.",
        message_type = "warning",
        verbose = verbose
      )
      layer <- "counts"
    }
    return(RunDEtest_pseudobulk(
      srt = srt,
      group.by = group.by,
      group1 = group1,
      group2 = group2,
      features = features,
      feature_type = feature_type,
      markers_type = markers_type,
      test.use = test.use,
      only.pos = only.pos,
      fc.threshold = fc.threshold,
      base = base,
      sample_col = sample_col,
      condition_col = condition_col,
      p.adjust.method = p.adjust.method,
      layer = layer,
      assay = assay,
      verbose = verbose,
      cores = cores,
      ...
    ))
  }

  if (!is.null(sample_col) || !is.null(condition_col)) {
    log_message(
      "{.arg sample_col} and {.arg condition_col} are only used by {.val edgeR} and {.val limma}. Ignoring them for cell-level differential testing.",
      message_type = "warning",
      verbose = verbose
    )
    sample_col <- NULL
    condition_col <- NULL
  }
  if (markers_type %in% c("conserved", "disturbed")) {
    if (is.null(grouping.var)) {
      log_message(
        "'grouping.var' must be provided when finding conserved or disturbed markers",
        message_type = "error"
      )
    }
  }

  skip_presto_check <- identical(test.use, "wilcox") &&
    markers_type %in% c("all", "paired") &&
    is.null(latent.vars) &&
    is.null(mean.fxn) &&
    layer %in% c("data", "counts") &&
    identical(norm.method, "LogNormalize") &&
    length(list(...)) == 0
  if (!isTRUE(skip_presto_check)) {
    check_r("immunogenomics/presto", verbose = FALSE)
  }

  status <- CheckDataType(srt, layer = layer, assay = assay)
  if (layer == "counts" && status != "raw_counts") {
    log_message(
      "Data in the {.arg counts} layer is not raw counts",
      message_type = "error"
    )
  }
  if (layer == "data" && status != "log_normalized_counts") {
    if (status == "raw_counts") {
      log_message(
        "Data in the {.arg data} layer is raw counts. Perform {.fun NormalizeData}({.val LogNormalize})",
        message_type = "warning",
        verbose = verbose
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
    if (status == "raw_normalized_counts") {
      log_message(
        "Data in the {.arg data} layer is raw_normalized_counts. Perform {.fun NormalizeData}({.val LogNormalize})",
        message_type = "warning",
        verbose = verbose
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
    if (status == "unknown") {
      log_message(
        "Data in the 'data' layer is unknown. Please check the data type",
        message_type = "warning",
        verbose = verbose
      )
    }
  }

  log_message(
    "Start differential expression test",
    verbose = verbose
  )
  if (fc.threshold < 1) {
    log_message(
      "{.arg fc.threshold} must be greater than or equal to 1",
      message_type = "error"
    )
  }

  if (!is.null(cells1) || !is.null(group1)) {
    if (is.null(cells1)) {
      if (is.null(group.by)) {
        log_message(
          "{.arg group.by} must be provided when {.arg group1} specified",
          message_type = "error"
        )
      }
      cells1 <- colnames(srt)[srt[[group.by, drop = TRUE]] %in% group1]
    }
    if (is.null(cells2) && !is.null(group2)) {
      cells2 <- colnames(srt)[srt[[group.by, drop = TRUE]] %in% group2]
    }
    if (!all(cells1 %in% colnames(srt))) {
      log_message(
        "{.arg cells1} has some cells not in {.cls Seurat}",
        message_type = "error"
      )
    }
    if (is.null(cells2)) {
      cells2 <- colnames(srt)[!colnames(srt) %in% cells1]
      group2 <- "others"
    }
    if (!all(cells2 %in% colnames(srt))) {
      log_message(
        "{.arg cells2} has some cells not in {.cls Seurat}",
        message_type = "error"
      )
    }
    if (length(cells1) < 3 || length(cells2) < 3) {
      log_message(
        "Cell groups must have more than 3 cells",
        message_type = "error"
      )
    }

    log_message(
      "Find ",
      markers_type,
      " markers(",
      test.use,
      ") for custom cell groups...",
      verbose = verbose
    )

    if (markers_type == "all") {
      markers <- RunDEtestFindMarkers(
        srt = srt,
        assay = assay,
        layer = layer,
        cells.1 = cells1,
        cells.2 = cells2,
        features = features,
        test.use = test.use,
        logfc.threshold = log(fc.threshold, base = base),
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = FALSE,
        ...
      )

      if (!is.null(markers) && nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        group1_str <- if (is.null(group1)) {
          "group1"
        } else if (length(group1) > 1) {
          paste(group1, collapse = ";")
        } else {
          as.character(group1)
        }
        group2_str <- if (is.null(group2)) {
          "group2"
        } else if (length(group2) > 1) {
          paste(group2, collapse = ";")
        } else {
          as.character(group2)
        }
        markers[, "group1"] <- group1_str
        markers[, "group2"] <- group2_str
        rownames(markers) <- NULL
        markers[, "group1"] <- factor(
          markers[, "group1"],
          levels = unique(markers[, "group1"])
        )
        if ("p_val" %in% colnames(markers)) {
          markers[, "p_val_adj"] <- stats::p.adjust(
            markers[, "p_val"],
            method = p.adjust.method
          )
        }
        markers[, "test_group_number"] <- 1L
        markers[, "test_group"] <- as.character(markers[, "group1"])
        srt@tools[["DEtest_custom"]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
        srt@tools[["DEtest_custom"]][["cells2"]] <- cells2
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
      }
    }

    if (markers_type == "conserved") {
      markers <- FindConservedMarkers2(
        object = srt,
        assay = assay,
        layer = layer,
        cells.1 = cells1,
        cells.2 = cells2,
        features = features,
        grouping.var = grouping.var,
        test.use = test.use,
        logfc.threshold = log(fc.threshold, base = base),
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        meta.method = meta.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = FALSE,
        ...
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        group1_str <- if (is.null(group1)) {
          "group1"
        } else if (length(group1) > 1) {
          paste(group1, collapse = ";")
        } else {
          as.character(group1)
        }
        group2_str <- if (is.null(group2)) {
          "group2"
        } else if (length(group2) > 1) {
          paste(group2, collapse = ";")
        } else {
          as.character(group2)
        }
        markers[, "group1"] <- group1_str
        markers[, "group2"] <- group2_str
        rownames(markers) <- NULL
        markers[, "group1"] <- factor(
          markers[, "group1"],
          levels = unique(markers[, "group1"])
        )
        if ("p_val" %in% colnames(markers)) {
          markers[, "p_val_adj"] <- stats::p.adjust(
            markers[, "p_val"],
            method = p.adjust.method
          )
        }
        markers[, "test_group_number"] <- as.integer(
          table(markers[["gene"]])[markers[, "gene"]]
        )
        markers_matrix <- as.data.frame.matrix(
          table(markers[, c("gene", "group1")])
        )
        markers[, "test_group"] <- apply(markers_matrix, 1, function(x) {
          paste0(colnames(markers_matrix)[x > 0], collapse = ";")
        })[markers[, "gene"]]
        srt@tools[["DEtest_custom"]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
        srt@tools[["DEtest_custom"]][["cells2"]] <- cells2
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
      }
    }
    if (markers_type == "disturbed") {
      srt_tmp <- srt
      srt_tmp[[grouping.var, drop = TRUE]][setdiff(
        colnames(srt_tmp),
        cells1
      )] <- NA
      srt_tmp <- RunDEtest(
        object = srt_tmp,
        assay = assay,
        layer = layer,
        group.by = grouping.var,
        markers_type = "all",
        features = features,
        test.use = test.use,
        fc.threshold = fc.threshold,
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        p.adjust.method = p.adjust.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        cores = cores,
        seed = seed,
        verbose = FALSE,
        ...
      )
      markers <- srt_tmp@tools[[paste0("DEtest_", grouping.var)]][[paste0(
        "AllMarkers_",
        test.use
      )]]
      if (!is.null(markers) && nrow(markers) > 0) {
        colnames(markers) <- gsub("group", "var", colnames(markers))
        markers[["group1"]] <- group1 %||% "group1"
        srt@tools[["DEtest_custom"]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
      } else {
        log_message(
          "No markers found",
          message_type = "warning",
          verbose = verbose
        )
      }
    }
  } else {
    if (is.null(group.by)) {
      cell_group <- Seurat::Idents(srt)
      group.by <- "active.ident"
    } else {
      cell_group <- srt[[group.by, drop = TRUE]]
    }
    if (!is.factor(cell_group)) {
      cell_group <- factor(cell_group, levels = unique(cell_group))
    }
    cell_group <- lapply(levels(cell_group), function(x) {
      cell <- cell_group[cell_group == x]
      out <- sample(
        cell,
        size = min(max.cells.per.ident, length(cell)),
        replace = FALSE
      )
      return(out)
    })
    cell_group <- stats::setNames(
      unlist(
        lapply(cell_group, function(x) x),
        use.names = FALSE
      ),
      unlist(lapply(cell_group, names))
    )

    args1 <- list(
      srt = srt,
      assay = assay,
      layer = layer,
      features = features,
      test.use = test.use,
      logfc.threshold = log(fc.threshold, base = base),
      base = base,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      max.cells.per.ident = max.cells.per.ident,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      latent.vars = latent.vars,
      only.pos = only.pos,
      norm.method = norm.method,
      pseudocount.use = pseudocount.use,
      mean.fxn = mean.fxn,
      verbose = FALSE,
      ...
    )

    log_message(
      "Find ",
      markers_type,
      " markers(",
      test.use,
      ") among ",
      nlevels(cell_group),
      " groups..."
    )

    if (markers_type == "all") {
      AllMarkers <- parallelize_fun(
        levels(cell_group),
        function(group) {
          cells.1 <- names(cell_group)[which(cell_group == group)]
          cells.2 <- names(cell_group)[which(cell_group != group)]
          if (length(cells.1) < 3 || length(cells.2) < 3) {
            return(NULL)
          } else {
            args1[["cells.1"]] <- cells.1
            args1[["cells.2"]] <- cells.2
            markers <- do.call(RunDEtestFindMarkers, args1)
            if (!is.null(markers) && nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(group)
              markers[, "group2"] <- "others"
              if ("p_val" %in% colnames(markers)) {
                markers[, "p_val_adj"] <- stats::p.adjust(
                  markers[, "p_val"],
                  method = p.adjust.method
                )
              }
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )
      AllMarkers <- do.call(rbind.data.frame, AllMarkers)
      if (!is.null(AllMarkers) && nrow(AllMarkers) > 0) {
        rownames(AllMarkers) <- NULL
        AllMarkers[, "group1"] <- factor(
          AllMarkers[, "group1"],
          levels = levels(cell_group)
        )
        AllMarkers[, "test_group_number"] <- as.integer(
          table(AllMarkers[["gene"]])[AllMarkers[, "gene"]]
        )
        AllMarkersMatrix <- as.data.frame.matrix(
          table(AllMarkers[, c("gene", "group1")])
        )
        AllMarkers[, "test_group"] <- apply(
          AllMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(AllMarkersMatrix)[x > 0], collapse = ";")
          }
        )[AllMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- AllMarkers
      } else {
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "AllMarkers_",
          test.use
        )]] <- data.frame()
      }
    }

    if (markers_type == "paired") {
      pair <- expand.grid(x = levels(cell_group), y = levels(cell_group))
      pair <- pair[pair[, 1] != pair[, 2], , drop = FALSE]
      cell_index <- split(names(cell_group), cell_group)
      PairedMarkers <- parallelize_fun(
        seq_len(nrow(pair)),
        function(i) {
          cells.1 <- cell_index[[as.character(pair[i, 1])]]
          cells.2 <- cell_index[[as.character(pair[i, 2])]]
          if (length(cells.1) < 3 || length(cells.2) < 3) {
            return(NULL)
          } else {
            args1[["cells.1"]] <- cells.1
            args1[["cells.2"]] <- cells.2
            markers <- do.call(RunDEtestFindMarkers, args1)
            if (!is.null(markers) && nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(pair[i, 1])
              markers[, "group2"] <- as.character(pair[i, 2])
              if ("p_val" %in% colnames(markers)) {
                markers[, "p_val_adj"] <- stats::p.adjust(
                  markers[, "p_val"],
                  method = p.adjust.method
                )
              }
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )
      PairedMarkers <- do.call(rbind.data.frame, PairedMarkers)
      if (!is.null(PairedMarkers) && nrow(PairedMarkers) > 0) {
        rownames(PairedMarkers) <- NULL
        PairedMarkers[, "group1"] <- factor(
          PairedMarkers[, "group1"],
          levels = levels(cell_group)
        )
        PairedMarkers[, "test_group_number"] <- as.integer(
          table(PairedMarkers[["gene"]])[PairedMarkers[, "gene"]]
        )
        PairedMarkersMatrix <- as.data.frame.matrix(
          table(
            PairedMarkers[, c("gene", "group1")]
          )
        )
        PairedMarkers[, "test_group"] <- apply(
          PairedMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(PairedMarkersMatrix)[x > 0], collapse = ";")
          }
        )[PairedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkers_",
          test.use
        )]] <- PairedMarkers
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkersMatrix_",
          test.use
        )]] <- PairedMarkersMatrix
      } else {
        log_message(
          "No markers found",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkers_",
          test.use
        )]] <- data.frame()
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "PairedMarkersMatrix_",
          test.use
        )]] <- NULL
      }
    }

    if (markers_type == "conserved") {
      ConservedMarkers <- parallelize_fun(
        levels(cell_group),
        function(group) {
          cells.1 <- names(cell_group)[which(cell_group == group)]
          cells.2 <- names(cell_group)[which(cell_group != group)]
          if (length(cells.1) < 3 || length(cells.2) < 3) {
            return(NULL)
          } else {
            args1[["cells.1"]] <- cells.1
            args1[["cells.2"]] <- cells.2
            args1[["object"]] <- srt
            args1[["assay"]] <- assay
            args1[["grouping.var"]] <- grouping.var
            args1[["meta.method"]] <- meta.method
            markers <- do.call(FindConservedMarkers2, args1)
            if (!is.null(markers) && nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(group)
              markers[, "group2"] <- "others"
              if ("p_val" %in% colnames(markers)) {
                markers[, "p_val_adj"] <- stats::p.adjust(
                  markers[, "p_val"],
                  method = p.adjust.method
                )
              }
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )
      ConservedMarkers <- do.call(
        rbind.data.frame,
        lapply(
          ConservedMarkers,
          function(x) {
            x[, c(
              "avg_log2FC",
              "pct.1",
              "pct.2",
              "max_pval",
              "p_val",
              "p_val_adj",
              "gene",
              "group1",
              "group2"
            )]
          }
        )
      )
      if (!is.null(ConservedMarkers) && nrow(ConservedMarkers) > 0) {
        rownames(ConservedMarkers) <- NULL
        ConservedMarkers[, "group1"] <- factor(
          ConservedMarkers[, "group1"],
          levels = levels(cell_group)
        )
        ConservedMarkers[
          ,
          "test_group_number"
        ] <- as.integer(table(ConservedMarkers[["gene"]])[ConservedMarkers[
          ,
          "gene"
        ]])
        ConservedMarkersMatrix <- as.data.frame.matrix(table(ConservedMarkers[, c(
          "gene",
          "group1"
        )]))
        ConservedMarkers[, "test_group"] <- apply(
          ConservedMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(ConservedMarkersMatrix)[x > 0], collapse = ";")
          }
        )[ConservedMarkers[, "gene"]]
        ConservedMarkers <- ConservedMarkers[, c(
          "avg_log2FC",
          "pct.1",
          "pct.2",
          "max_pval",
          "p_val",
          "p_val_adj",
          "gene",
          "group1",
          "group2",
          "test_group_number",
          "test_group"
        )]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- ConservedMarkers
      } else {
        log_message(
          "No markers found",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "ConservedMarkers_",
          test.use
        )]] <- data.frame()
      }
    }

    if (markers_type == "disturbed") {
      DisturbedMarkers <- parallelize_fun(
        levels(cell_group),
        function(group) {
          cells.1 <- names(cell_group)[which(cell_group == group)]
          srt_tmp <- srt
          srt_tmp[[grouping.var, drop = TRUE]][setdiff(
            colnames(srt_tmp),
            cells.1
          )] <- NA
          if (
            length(stats::na.omit(unique(srt_tmp[[
              grouping.var,
              drop = TRUE
            ]]))) <
              2
          ) {
            return(NULL)
          } else {
            srt_tmp <- RunDEtest(
              object = srt_tmp,
              assay = assay,
              layer = layer,
              group.by = grouping.var,
              markers_type = "all",
              features = features,
              test.use = test.use,
              fc.threshold = fc.threshold,
              base = base,
              min.pct = min.pct,
              min.diff.pct = min.diff.pct,
              max.cells.per.ident = max.cells.per.ident,
              min.cells.feature = min.cells.feature,
              min.cells.group = min.cells.group,
              latent.vars = latent.vars,
              only.pos = only.pos,
              norm.method = norm.method,
              p.adjust.method = p.adjust.method,
              pseudocount.use = pseudocount.use,
              mean.fxn = mean.fxn,
              cores = cores,
              seed = seed,
              verbose = FALSE,
              ...
            )
            markers <- srt_tmp@tools[[paste0("DEtest_", grouping.var)]][[paste0(
              "AllMarkers_",
              test.use
            )]]
            if (!is.null(markers) && nrow(markers) > 0) {
              colnames(markers) <- gsub("group", "var", colnames(markers))
              markers[["group1"]] <- as.character(group)
              return(markers)
            } else {
              return(NULL)
            }
          }
        },
        cores = cores,
        verbose = verbose
      )

      DisturbedMarkers <- do.call(rbind.data.frame, DisturbedMarkers)
      if (!is.null(DisturbedMarkers) && nrow(DisturbedMarkers) > 0) {
        rownames(DisturbedMarkers) <- NULL
        DisturbedMarkers[, "group1"] <- factor(
          DisturbedMarkers[, "group1"],
          levels = levels(cell_group)
        )
        DisturbedMarkers[
          ,
          "test_group_number"
        ] <- as.integer(table(unique(DisturbedMarkers[, c("gene", "group1")])[[
          "gene"
        ]])[DisturbedMarkers[, "gene"]])
        DisturbedMarkersMatrix <- as.data.frame.matrix(table(DisturbedMarkers[, c(
          "gene",
          "group1"
        )]))
        DisturbedMarkers[, "test_group"] <- apply(
          DisturbedMarkersMatrix,
          1,
          function(x) {
            paste0(colnames(DisturbedMarkersMatrix)[x > 0], collapse = ";")
          }
        )[DisturbedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- DisturbedMarkers
      } else {
        log_message(
          "No markers found.",
          message_type = "warning",
          verbose = verbose
        )
        srt@tools[[paste0("DEtest_", group.by)]][[paste0(
          "DisturbedMarkers_",
          test.use
        )]] <- data.frame()
      }
    }
  }

  log_message(
    "Differential expression test completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

#' @rdname RunDEtest
#' @export
RunDEtest.SummarizedExperiment <- function(
  object,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  features = NULL,
  feature_type = c("gene", "peak", "cCRE"),
  markers_type = c(
    "all",
    "paired",
    "conserved",
    "disturbed"
  ),
  grouping.var = NULL,
  meta.method = c(
    "maximump",
    "minimump",
    "wilkinsonp",
    "meanp",
    "sump",
    "votep"
  ),
  test.use = "edgeR",
  only.pos = TRUE,
  fc.threshold = 1.5,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  norm.method = "LogNormalize",
  sample_col = NULL,
  condition_col = NULL,
  bulk_assay = "counts",
  p.adjust.method = "bonferroni",
  layer = "counts",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
) {
  set.seed(seed)
  feature_type <- match.arg(feature_type)
  markers_type <- match.arg(markers_type)
  meta.method <- match.arg(meta.method)
  bulk_se <- object

  if (!identical(feature_type, "gene")) {
    log_message(
      "{.arg feature_type} currently supports only {.val gene} for {.cls SummarizedExperiment} input.",
      message_type = "error"
    )
  }
  if (markers_type %in% c("conserved", "disturbed")) {
    log_message(
      "{.arg markers_type} values {.val conserved} and {.val disturbed} are not supported for {.cls SummarizedExperiment} input.",
      message_type = "error"
    )
  }
  if (!is.null(cells1) || !is.null(cells2)) {
    log_message(
      "{.arg cells1} and {.arg cells2} are not supported for {.cls SummarizedExperiment} input.",
      message_type = "error"
    )
  }
  if (!is.null(sample_col)) {
    log_message(
      "{.arg sample_col} is ignored for {.cls SummarizedExperiment} input because bulk sample IDs are taken from column names.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (!identical(layer, "counts")) {
    log_message(
      "{.arg layer} is ignored for {.cls SummarizedExperiment} input. Use {.arg bulk_assay} instead.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (!is.null(assay)) {
    log_message(
      "{.arg assay} is ignored for {.cls SummarizedExperiment} input. Use {.arg bulk_assay} instead.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (!is.null(grouping.var)) {
    log_message(
      "{.arg grouping.var} is not supported for {.cls SummarizedExperiment} input.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (
    !is.null(latent.vars) ||
      !isTRUE(identical(norm.method, "LogNormalize")) ||
      !isTRUE(identical(pseudocount.use, 1)) ||
      !isTRUE(identical(min.pct, 0.1)) ||
      !isTRUE(identical(min.diff.pct, -Inf)) ||
      !isTRUE(identical(max.cells.per.ident, Inf)) ||
      !isTRUE(identical(min.cells.group, 3)) ||
      !isTRUE(identical(meta.method, "maximump")) ||
      !is.null(mean.fxn)
  ) {
    log_message(
      "Several cell-level specific DE arguments are ignored for {.cls SummarizedExperiment} input.",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (is.null(condition_col) || !nzchar(condition_col)) {
    log_message(
      "{.arg condition_col} is required for {.cls SummarizedExperiment} input.",
      message_type = "error"
    )
  }
  if (!test.use %in% detest_sample_methods()) {
    log_message(
      paste0(
        "{.arg test.use} must be one of ",
        paste(detest_sample_methods(), collapse = ", "),
        " for {.cls SummarizedExperiment} input."
      ),
      message_type = "error"
    )
  }
  if (fc.threshold < 1) {
    log_message(
      "{.arg fc.threshold} must be greater than or equal to 1",
      message_type = "error"
    )
  }

  ctx <- build_context(
    mode = "pure_bulk",
    bulk_se = bulk_se,
    condition.by = condition_col,
    group.by = group.by,
    bulk_assay = bulk_assay
  )
  if (!is.null(features)) {
    features <- intersect(as.character(features), rownames(ctx$counts_global))
    if (length(features) == 0) {
      log_message(
        "No requested {.arg features} were found in the bulk assay.",
        message_type = "error"
      )
    }
    ctx$counts_global <- ctx$counts_global[features, , drop = FALSE]
    ctx$counts_by_group <- lapply(
      ctx$counts_by_group,
      function(x) x[features, , drop = FALSE]
    )
  }

  de_markers_type <- if (identical(markers_type, "paired")) {
    "paired"
  } else {
    unique_condition <- unique(as.character(ctx$condition_global))
    if (!is.null(group1) || !is.null(group2) || length(unique_condition) <= 2) {
      "single"
    } else {
      "all"
    }
  }
  method_name <- detest_method_key(test.use)
  bundle <- dispatch_de(
    ctx = ctx,
    method_name = method_name,
    condition1 = group1,
    condition2 = group2,
    de_markers_type = de_markers_type,
    de_args = utils::modifyList(
      list(
        only.pos = only.pos,
        logfc.threshold = log(fc.threshold, base = base),
        p.adjust.method = p.adjust.method
      ),
      list(...)
    ),
    verbose = verbose
  )
  bundle$results <- coerce_de_schema(bundle$results)
  store <- list(
    input = list(
      condition_col = condition_col,
      group.by = group.by,
      bulk_assay = bulk_assay,
      features = features
    ),
    active_method = method_name,
    test.use = test.use,
    methods = stats::setNames(list(bundle), method_name),
    results = bundle$results,
    parameters = utils::modifyList(
      bundle$parameters %||% list(),
      list(
        test.use = test.use,
        condition1 = group1,
        condition2 = group2,
        markers_type = de_markers_type
      )
    ),
    status = list(
      method = method_name,
      status = bundle$status %||% "failed",
      reason = bundle$reason %||% NULL
    )
  )
  store_meta(bulk_se, "DEtest", store)
}

RunDEtest_DESeq2 <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  only.pos = TRUE,
  logfc.threshold = 0,
  p.adjust.method = "bonferroni"
) {
  check_r("DESeq2", verbose = FALSE)
  DESeqDataSetFromMatrix <- get_namespace_fun(
    "DESeq2",
    "DESeqDataSetFromMatrix"
  )
  DESeq <- get_namespace_fun("DESeq2", "DESeq")
  DESeq_results <- get_namespace_fun("DESeq2", "results")

  condition_all <- as.character(condition)
  if (is.null(condition1) || is.null(condition2)) {
    condition_levels <- unique(condition_all)
    if (length(condition_levels) < 2) {
      return(NULL)
    }
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
  }

  keep <- condition_all %in% c(condition1, condition2)
  count_use <- as.matrix(count_matrix[, keep, drop = FALSE])
  cond_use <- factor(condition_all[keep], levels = c(condition1, condition2))
  if (ncol(count_use) < 2 || any(table(cond_use) < 2)) {
    return(NULL)
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
        contrast = c("condition", condition2, condition1)
      )
      res_df <- as.data.frame(res)
      if (nrow(res_df) == 0) {
        return(data.frame())
      }
      detect <- count_use[rownames(res_df), , drop = FALSE] > 0
      pct.1 <- round(
        rowMeans(detect[, cond_use == condition1, drop = FALSE]),
        3
      )
      pct.2 <- round(
        rowMeans(detect[, cond_use == condition2, drop = FALSE]),
        3
      )
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
      out_df$p_val_adj <- stats::p.adjust(
        out_df$p_val,
        method = p.adjust.method
      )
      if (isTRUE(only.pos)) {
        out_df <- out_df[out_df$avg_log2FC >= logfc.threshold, , drop = FALSE]
      } else {
        out_df <- out_df[
          abs(out_df$avg_log2FC) >= logfc.threshold, ,
          drop = FALSE
        ]
      }
      out_df
    },
    error = function(e) NULL
  )
  out
}

RunDEtest_dream <- function(
  count_matrix,
  condition,
  condition1 = NULL,
  condition2 = NULL,
  sample_data = NULL,
  dream_formula = NULL,
  only.pos = TRUE,
  logfc.threshold = 0,
  p.adjust.method = "bonferroni"
) {
  check_r(c("variancePartition", "limma", "edgeR"), verbose = FALSE)
  DGEList <- get_namespace_fun("edgeR", "DGEList")
  filterByExpr <- get_namespace_fun("edgeR", "filterByExpr")
  calcNormFactors <- get_namespace_fun("edgeR", "calcNormFactors")
  voomWithDreamWeights <- get_namespace_fun(
    "variancePartition",
    "voomWithDreamWeights"
  )
  dream <- get_namespace_fun("variancePartition", "dream")
  eBayes <- get_namespace_fun("limma", "eBayes")
  topTable <- get_namespace_fun("limma", "topTable")

  condition_all <- as.character(condition)
  if (is.null(condition1) || is.null(condition2)) {
    condition_levels <- unique(condition_all)
    if (length(condition_levels) < 2) {
      return(NULL)
    }
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
  }

  keep <- condition_all %in% c(condition1, condition2)
  count_use <- as.matrix(count_matrix[, keep, drop = FALSE])
  cond_use <- factor(condition_all[keep], levels = c(condition1, condition2))
  if (ncol(count_use) < 2 || any(table(cond_use) < 2)) {
    return(NULL)
  }

  sample_df <- if (is.null(sample_data)) {
    data.frame(row.names = colnames(count_use))
  } else {
    if (!all(colnames(count_use) %in% rownames(sample_data))) {
      return(NULL)
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
  if (
    !inherits(dream_formula, "formula") ||
      !"condition" %in% all.vars(dream_formula)
  ) {
    return(NULL)
  }

  out <- tryCatch(
    {
      dge <- DGEList(counts = count_use)
      keep_features <- filterByExpr(dge, group = sample_df$condition)
      if (!any(keep_features)) {
        return(data.frame())
      }
      dge <- dge[keep_features, , keep.lib.sizes = FALSE]
      dge <- calcNormFactors(dge)
      vobj <- voomWithDreamWeights(
        dge,
        formula = dream_formula,
        data = sample_df
      )
      fit <- dream(exprObj = vobj, formula = dream_formula, data = sample_df)
      fit <- eBayes(fit)
      coef_names <- colnames(fit$coefficients)
      coef_use <- grep(
        pattern = paste0("^condition", make.names(condition2), "$"),
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
        return(data.frame())
      }
      detect <- dge$counts[rownames(tt), , drop = FALSE] > 0
      pct.1 <- round(
        rowMeans(detect[, sample_df$condition == condition1, drop = FALSE]),
        3
      )
      pct.2 <- round(
        rowMeans(detect[, sample_df$condition == condition2, drop = FALSE]),
        3
      )
      out_df <- data.frame(
        p_val = tt$P.Value,
        avg_log2FC = tt$logFC,
        ave_expr = tt$AveExpr,
        pct.1 = pct.1,
        pct.2 = pct.2,
        row.names = rownames(tt)
      )
      out_df$p_val_adj <- stats::p.adjust(
        out_df$p_val,
        method = p.adjust.method
      )
      if (isTRUE(only.pos)) {
        out_df <- out_df[out_df$avg_log2FC >= logfc.threshold, , drop = FALSE]
      } else {
        out_df <- out_df[
          abs(out_df$avg_log2FC) >= logfc.threshold, ,
          drop = FALSE
        ]
      }
      out_df
    },
    error = function(e) NULL
  )
  out
}

detest_sample_methods <- function() {
  c("edgeR", "limma", "DESeq2", "dream")
}

detest_method_key <- function(test.use) {
  key <- tolower(as.character(test.use) %||% "")
  method_map <- c(
    limma = "de_limma_voom",
    edger = "de_edgeR_qlf",
    deseq2 = "de_DESeq2",
    dream = "de_dream"
  )
  out <- unname(method_map[[key]])
  if (is.null(out) || is.na(out)) {
    log_message(
      paste0(
        "Unsupported bulk DE {.arg test.use}: ",
        test.use,
        ". Use one of limma, edgeR, DESeq2, or dream."
      ),
      message_type = "error"
    )
  }
  out
}

detest_key_method <- function(method_name) {
  key <- as.character(method_name) %||% ""
  method_map <- c(
    de_limma_voom = "limma",
    de_edgeR_qlf = "edgeR",
    de_DESeq2 = "DESeq2",
    de_dream = "dream"
  )
  unname(method_map[[key]]) %||% key
}

store_meta <- function(object, name, value) {
  meta <- S4Vectors::metadata(object)
  meta[[name]] <- value
  S4Vectors::metadata(object) <- meta
  object
}

detest_res <- function(res) {
  if (is.null(res)) {
    return(coerce_de_schema(data.frame()))
  }
  out <- as.data.frame(res, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(out) && nrow(out) > 0 && !is.null(rownames(out))) {
    out$gene <- rownames(out)
  }
  coerce_de_schema(out)
}

resolve_detest_result <- function(
  object = NULL,
  group.by = NULL,
  test.use = "wilcox",
  res = NULL
) {
  if (!is.null(res)) {
    return(detest_res(res))
  }

  if (inherits(object, "Seurat")) {
    group.by <- group.by %||% "custom"
    layer <- paste0("DEtest_", group.by)
    if (
      !layer %in% names(object@tools) ||
        length(grep(pattern = "AllMarkers", names(object@tools[[layer]]))) == 0
    ) {
      log_message(
        "Cannot find the DEtest result for the group {.val {group.by}}. Perform {.fn RunDEtest} first",
        message_type = "error"
      )
    }
    index <- grep(
      pattern = paste0("AllMarkers_", test.use),
      names(object@tools[[layer]])
    )[1]
    if (is.na(index)) {
      log_message(
        "Cannot find the {.val AllMarkers_{test.use}} in the DEtest result",
        message_type = "error"
      )
    }
    return(detest_res(object@tools[[layer]][[index]]))
  }

  if (inherits(object, "SummarizedExperiment")) {
    bundle <- S4Vectors::metadata(object)[["DEtest"]]
    if (is.null(bundle)) {
      log_message(
        "Cannot find bulk DEtest results in {.fn S4Vectors::metadata(object)[['DEtest']]}",
        message_type = "error"
      )
    }
    if (!is.null(bundle$results) && nrow(bundle$results) > 0) {
      return(detest_res(bundle$results))
    }
    if (!is.null(bundle$methods) && length(bundle$methods) > 0) {
      methods <- bundle$methods
      method_use <- detest_method_key(test.use)
      if (!method_use %in% names(methods)) {
        method_use <- names(methods)[1]
      }
      return(detest_res(methods[[method_use]][["results"]]))
    }
    return(detest_res(data.frame()))
  }

  log_message(
    "{.arg object} must be a {.cls Seurat} or {.cls SummarizedExperiment} object when {.arg res} is NULL.",
    message_type = "error"
  )
}

resolve_deconvolution_result <- function(object = NULL, res = NULL) {
  if (!is.null(res)) {
    return(deconv_schema(as.data.frame(res, stringsAsFactors = FALSE)))
  }
  if (!inherits(object, "SummarizedExperiment")) {
    log_message(
      "{.arg object} must be a {.cls SummarizedExperiment} object when {.arg res} is NULL.",
      message_type = "error"
    )
  }
  bundle <- S4Vectors::metadata(object)[["Deconvolution"]]
  if (is.null(bundle)) {
    log_message(
      "Cannot find deconvolution results in {.fn S4Vectors::metadata(object)[['Deconvolution']]}",
      message_type = "error"
    )
  }
  deconv_schema(bundle$results %||% data.frame())
}
