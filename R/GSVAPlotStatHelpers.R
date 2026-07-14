gsva_plot_resolve_scores <- function(
  srt = NULL,
  res = NULL,
  assay_name = "GSVA",
  group.by = NULL
) {
  enrichment <- NULL
  is_group_level <- FALSE

  if (!is.null(res)) {
    if (is.list(res) && "scores" %in% names(res)) {
      scores <- res[["scores"]]
      enrichment <- res[["enrichment"]]
      is_group_level <- !is.null(res[["group.by"]])
    } else if (is.matrix(res) || inherits(res, "Matrix")) {
      scores <- res
    } else {
      log_message(
        "Invalid {.arg res} format",
        message_type = "error"
      )
    }
  } else {
    if (is.null(srt) || !inherits(srt, "Seurat")) {
      log_message(
        "Either {.arg srt} or {.arg res} must be provided",
        message_type = "error"
      )
    }
    if (assay_name %in% SeuratObject::Assays(srt)) {
      scores <- GetAssayData5(srt, assay = assay_name, layer = "data")
    } else {
      tool_name <- gsva_plot_find_tool_name(srt, assay_name, group.by)
      if (is.null(tool_name)) {
        log_message(
          "GSVA results not found. Please run RunGSVA first",
          message_type = "error"
        )
      }
      gsva_result <- srt@tools[[tool_name]]
      scores <- gsva_result[["scores"]]
      enrichment <- gsva_result[["enrichment"]]
      is_group_level <- !is.null(gsva_result[["group.by"]])
    }
  }

  if (!is.matrix(scores)) {
    scores <- as.matrix(scores)
  }
  if (is.null(rownames(scores)) || is.null(colnames(scores))) {
    log_message(
      "GSVA scores must have row and column names",
      message_type = "error"
    )
  }
  list(
    scores = scores,
    enrichment = enrichment,
    is_group_level = isTRUE(is_group_level)
  )
}

gsva_plot_find_tool_name <- function(srt, assay_name, group.by = NULL) {
  candidates <- character(0)
  if (!is.null(group.by)) {
    candidates <- c(
      paste("GSVA", group.by, "gsva", sep = "_"),
      paste0("GSVA_", assay_name, "_", group.by),
      grep(paste0("^GSVA_", group.by, "_"), names(srt@tools), value = TRUE)
    )
  }
  candidates <- c(
    paste("GSVA_cell", "gsva", sep = "_"),
    paste0("GSVA_", assay_name),
    grep("^GSVA_cell_", names(srt@tools), value = TRUE),
    grep("^GSVA_", names(srt@tools), value = TRUE)
  )
  candidates <- unique(candidates)
  hit <- candidates[candidates %in% names(srt@tools)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1]]
}

gsva_plot_match_features <- function(features, scores, enrichment = NULL) {
  features <- unique(as.character(unlist(features)))
  features <- features[!is.na(features) & nzchar(features)]
  direct <- intersect(features, rownames(scores))
  if (length(direct) > 0L) {
    return(direct)
  }
  if (is.null(enrichment) || nrow(enrichment) == 0L) {
    return(character(0))
  }
  enrichment <- as.data.frame(enrichment)
  id <- as.character(enrichment[["ID"]] %||% enrichment[["Description"]])
  desc <- as.character(enrichment[["Description"]] %||% enrichment[["ID"]])
  row_key <- unique(desc[id %in% features | desc %in% features])
  intersect(row_key, rownames(scores))
}

gsva_plot_term_info <- function(scores, enrichment = NULL) {
  info <- data.frame(
    ID = rownames(scores),
    Description = rownames(scores),
    Database = "GSVA",
    stringsAsFactors = FALSE
  )
  if (is.null(enrichment) || nrow(enrichment) == 0L) {
    return(info)
  }
  enrichment <- as.data.frame(enrichment)
  desc <- as.character(enrichment[["Description"]] %||% enrichment[["ID"]])
  id <- as.character(enrichment[["ID"]] %||% desc)
  db <- as.character(enrichment[["Database"]] %||% "GSVA")
  mapped <- data.frame(
    ID = id,
    Description = desc,
    Database = db,
    stringsAsFactors = FALSE
  )
  mapped <- mapped[!duplicated(mapped[["Description"]]), , drop = FALSE]
  idx <- match(info[["Description"]], mapped[["Description"]])
  hit <- !is.na(idx)
  info[["ID"]][hit] <- mapped[["ID"]][idx[hit]]
  info[["Database"]][hit] <- mapped[["Database"]][idx[hit]]
  info
}

gsva_plot_score_table <- function(
  scores,
  enrichment = NULL,
  srt = NULL,
  group.by = NULL,
  group_use = NULL,
  sort.by = "abs",
  topTerm = 20,
  score_cutoff = NULL
) {
  if (!is.null(srt) && !is.null(group.by) && all(colnames(scores) %in% colnames(srt))) {
    if (!group.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group.by} {.val {group.by}} not found in meta.data",
        message_type = "error"
      )
    }
    cells <- colnames(scores)
    group_raw <- srt@meta.data[cells, group.by, drop = TRUE]
    group <- as.character(group_raw)
    keep_cell <- !is.na(group) & nzchar(group)
    if (!is.null(group_use)) {
      keep_cell <- keep_cell & group %in% group_use
    }
    if (!any(keep_cell)) {
      log_message(
        "No cells remain after applying {.arg group.by} and {.arg group_use}",
        message_type = "error"
      )
    }
    scores_use <- scores[, keep_cell, drop = FALSE]
    group <- group[keep_cell]
    if (is.factor(group_raw)) {
      group_levels <- levels(group_raw)
      group_levels <- group_levels[group_levels %in% unique(group)]
    } else {
      group_levels <- unique(group)
    }
    rows <- lapply(group_levels, function(g) {
      mat <- scores_use[, group == g, drop = FALSE]
      data.frame(
        Description = rownames(scores_use),
        Groups = g,
        GSVA_Score = rowMeans(mat, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, rows)
  } else if (!is.null(enrichment) && nrow(enrichment) > 0L && "GSVA_Score" %in% colnames(enrichment)) {
    enrichment <- as.data.frame(enrichment)
    df <- data.frame(
      Description = as.character(enrichment[["Description"]] %||% enrichment[["ID"]]),
      Groups = as.character(enrichment[["Groups"]] %||% "GSVA"),
      GSVA_Score = as.numeric(enrichment[["GSVA_Score"]]),
      stringsAsFactors = FALSE
    )
    if (!is.null(group_use)) {
      df <- df[df[["Groups"]] %in% group_use, , drop = FALSE]
    }
  } else {
    idx <- max.col(abs(scores), ties.method = "first")
    df <- data.frame(
      Description = rownames(scores),
      Groups = colnames(scores)[idx],
      GSVA_Score = scores[cbind(seq_len(nrow(scores)), idx)],
      stringsAsFactors = FALSE
    )
  }

  info <- gsva_plot_term_info(scores, enrichment)
  df <- merge(df, info, by = "Description", all.x = TRUE, sort = FALSE)
  df <- df[, c("ID", "Description", "Database", "Groups", "GSVA_Score"), drop = FALSE]
  if (!is.null(score_cutoff)) {
    df <- df[abs(df[["GSVA_Score"]]) >= abs(score_cutoff), , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No gene sets pass the selected filters",
      message_type = "error"
    )
  }
  rank_metric <- if (identical(sort.by, "abs")) abs(df[["GSVA_Score"]]) else df[["GSVA_Score"]]
  df <- df[order(rank_metric, decreasing = TRUE, na.last = NA), , drop = FALSE]
  if (!is.null(topTerm) && is.finite(topTerm) && nrow(df) > topTerm) {
    df <- df[seq_len(topTerm), , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

gsva_plot_diff_table <- function(
  scores,
  enrichment = NULL,
  srt,
  group.by,
  sample.by,
  group_use = NULL,
  aggregate.fun = "mean",
  test.use = "wilcox",
  p.adjust.method = "BH",
  is_group_level = FALSE
) {
  if (is.null(srt) || !inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be provided when {.arg mode = 'diff'}",
      message_type = "error"
    )
  }
  if (is.null(group.by) || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  if (is.null(sample.by) || !sample.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg sample.by} must be one metadata column in {.arg srt} for sample-level GSVA statistics",
      message_type = "error"
    )
  }
  score_columns <- colnames(scores)
  sample_values <- as.character(srt@meta.data[[sample.by]])
  cell_level <- all(score_columns %in% colnames(srt))
  sample_level <- all(score_columns %in% unique(sample_values))
  if (isTRUE(is_group_level) || (!cell_level && !sample_level)) {
    log_message(
      paste0(
        "Differential GSVA statistics require cell-level or sample-level scores ",
        "that match {.arg srt} cells or {.arg sample.by} values. Results from {.fn RunGSVA} with ",
        "{.arg group.by} are already aggregated and cannot produce real p-values; ",
        "rerun {.fn RunGSVA} with {.arg group.by = NULL}."
      ),
      message_type = "error"
    )
  }

  if (cell_level) {
    meta <- srt@meta.data[score_columns, c(group.by, sample.by), drop = FALSE]
    group <- as.character(meta[[group.by]])
    sample <- as.character(meta[[sample.by]])
  } else {
    sample_group <- lapply(score_columns, function(s) {
      vals <- unique(as.character(srt@meta.data[sample_values == s, group.by, drop = TRUE]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      vals
    })
    if (any(lengths(sample_group) != 1L)) {
      log_message(
        "Each sample-level GSVA score column must map to exactly one {.arg group.by} value",
        message_type = "error"
      )
    }
    sample <- score_columns
    group <- vapply(sample_group, `[`, character(1), 1L)
  }
  keep <- !is.na(group) & nzchar(group) & !is.na(sample) & nzchar(sample)
  if (!is.null(group_use)) {
    group_use <- as.character(group_use)
    if (length(group_use) != 2L) {
      log_message(
        "{.arg group_use} must contain exactly two groups when {.arg mode = 'diff'}",
        message_type = "error"
      )
    }
  } else {
    if (is.factor(srt@meta.data[[group.by]])) {
      group_use <- levels(srt@meta.data[[group.by]])
      group_use <- group_use[group_use %in% unique(group[keep])]
    } else {
      group_use <- unique(group[keep])
    }
    group_use <- group_use[seq_len(min(2L, length(group_use)))]
  }
  if (length(group_use) != 2L) {
    log_message(
      "{.arg mode = 'diff'} requires exactly two non-missing groups",
      message_type = "error"
    )
  }
  keep <- keep & group %in% group_use
  if (!all(group_use %in% unique(group[keep]))) {
    log_message(
      "{.arg group_use} contains groups not present in {.arg group.by}",
      message_type = "error"
    )
  }
  group <- group[keep]
  sample <- sample[keep]
  score_columns <- score_columns[keep]
  if (length(score_columns) == 0L) {
    log_message(
      "No GSVA score columns remain after applying {.arg group_use}",
      message_type = "error"
    )
  }

  if (cell_level) {
    key <- paste(sample, group, sep = "\r")
    key_factor <- factor(key, levels = unique(key))
    x <- t(scores[, score_columns, drop = FALSE])
    if (identical(aggregate.fun, "mean")) {
      sums <- rowsum(x, key_factor, reorder = FALSE, na.rm = TRUE)
      counts <- as.numeric(tabulate(key_factor, nbins = nlevels(key_factor)))
      agg <- sweep(sums, 1, counts, "/")
    } else {
      split_index <- split(seq_along(key), key_factor)
      agg <- do.call(rbind, lapply(split_index, function(i) {
        apply(x[i, , drop = FALSE], 2, stats::median, na.rm = TRUE)
      }))
    }
    agg_meta <- data.frame(
      key = levels(key_factor),
      stringsAsFactors = FALSE
    )
    key_parts <- strsplit(agg_meta[["key"]], "\r", fixed = TRUE)
    agg_meta[["sample"]] <- vapply(key_parts, `[`, character(1), 1L)
    agg_meta[["group"]] <- vapply(key_parts, `[`, character(1), 2L)
  } else {
    agg <- t(scores[, score_columns, drop = FALSE])
    agg_meta <- data.frame(
      sample = sample,
      group = group,
      stringsAsFactors = FALSE
    )
  }

  n_by_group <- table(agg_meta[["group"]])
  if (any(n_by_group[group_use] < 2L)) {
    log_message(
      "Each group must contain at least two {.arg sample.by} levels for GSVA differential statistics",
      message_type = "error"
    )
  }

  g1 <- group_use[[1]]
  g2 <- group_use[[2]]
  idx1 <- agg_meta[["group"]] == g1
  idx2 <- agg_meta[["group"]] == g2
  pvalue <- vapply(seq_len(ncol(agg)), function(j) {
    x1 <- agg[idx1, j]
    x2 <- agg[idx2, j]
    x1 <- x1[is.finite(x1)]
    x2 <- x2[is.finite(x2)]
    if (length(x1) < 2L || length(x2) < 2L) {
      return(NA_real_)
    }
    if (identical(test.use, "wilcox")) {
      out <- try(stats::wilcox.test(x1, x2)$p.value, silent = TRUE)
    } else {
      out <- try(stats::t.test(x1, x2)$p.value, silent = TRUE)
    }
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
  }, numeric(1))
  mean1 <- colMeans(agg[idx1, , drop = FALSE], na.rm = TRUE)
  mean2 <- colMeans(agg[idx2, , drop = FALSE], na.rm = TRUE)

  info <- gsva_plot_term_info(scores, enrichment)
  out <- data.frame(
    ID = info[["ID"]],
    Description = info[["Description"]],
    Database = info[["Database"]],
    group1 = g1,
    group2 = g2,
    mean_group1 = as.numeric(mean1),
    mean_group2 = as.numeric(mean2),
    diff = as.numeric(mean1 - mean2),
    pvalue = pvalue,
    p.adjust = stats::p.adjust(pvalue, method = p.adjust.method),
    test.use = test.use,
    n_group1 = sum(idx1),
    n_group2 = sum(idx2),
    stringsAsFactors = FALSE
  )
  out <- out[order(out[["p.adjust"]], -abs(out[["diff"]]), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

gsva_plot_diff_bar <- function(
  df,
  topTerm = 20,
  palette = "simspec",
  palcolor = NULL,
  theme_use,
  theme_args
) {
  if (!is.null(topTerm) && is.finite(topTerm) && nrow(df) > topTerm) {
    df <- df[seq_len(topTerm), , drop = FALSE]
  }
  df[["term"]] <- factor(df[["Description"]], levels = rev(df[["Description"]]))
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["term"]], y = .data[["diff"]], fill = .data[["diff"]])
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::scale_fill_gradientn(
      colours = palette_colors(palette = palette, palcolor = palcolor),
      name = "Mean diff"
    ) +
    ggplot2::labs(x = NULL, y = paste0(unique(df[["group1"]]), " - ", unique(df[["group2"]]))) +
    do.call(theme_use, theme_args) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

gsva_plot_diff_volcano <- function(
  df,
  padjustCutoff = NULL,
  palette = "simspec",
  palcolor = NULL,
  theme_use,
  theme_args
) {
  cutoff <- padjustCutoff %||% 0.05
  df[["neglog10_padj"]] <- -log10(pmax(df[["p.adjust"]], .Machine$double.xmin))
  df[["significant"]] <- ifelse(df[["p.adjust"]] <= cutoff, "FDR", "NS")
  cols <- palette_colors(c("NS", "FDR"), palette = palette, palcolor = palcolor)
  if (is.null(names(cols)) || !all(c("NS", "FDR") %in% names(cols))) {
    cols <- unname(cols)[seq_len(min(2L, length(cols)))]
    if (length(cols) < 2L) {
      cols <- rep(cols, length.out = 2L)
    }
    names(cols) <- c("NS", "FDR")
  } else {
    cols <- cols[c("NS", "FDR")]
  }
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["diff"]], y = .data[["neglog10_padj"]], color = .data[["significant"]])
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = -log10(cutoff), color = "grey75", linewidth = 0.35, linetype = 2) +
    ggplot2::scale_color_manual(values = cols[c("NS", "FDR")], drop = FALSE) +
    ggplot2::labs(
      x = paste0(unique(df[["group1"]]), " - ", unique(df[["group2"]])),
      y = "-log10(adjusted p-value)",
      color = NULL
    ) +
    do.call(theme_use, theme_args)
}
