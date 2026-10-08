ccc_ov_palette_28 <- function() {
  c(
    "#1F577B", "#A56BA7", "#E0A7C8", "#7CBB5F", "#368650",
    "#01A0A7", "#75C8CC", "#E069A6", "#941456", "#FCBC10",
    "#EF7B77", "#279AD7", "#F0EEF0", "#EAEFC5", "#A499CC",
    "#5E4D9A", "#78C2ED", "#866017", "#9F987F", "#E0DFED",
    "#F0D7BC", "#D5B26C", "#D5DA48", "#B6B812", "#9DC3C3",
    "#A89C92", "#FEE00C", "#FEF2A1"
  )
}

ccc_ov_palette_56 <- function() {
  c(
    "#d60000", "#8c3bff", "#018700", "#00acc6", "#97ff00",
    "#ff7ed1", "#6b004f", "#ffa52f", "#00009c", "#857067",
    "#004942", "#4f2a00", "#00fdcf", "#bcb6ff", "#95b379",
    "#bf03b8", "#2466a1", "#280041", "#dbb3af", "#fdf490",
    "#4f445b", "#a37c00", "#ff7066", "#3f806e", "#82000c",
    "#a37bb3", "#344d00", "#9ae4ff", "#eb0077", "#2d000a",
    "#5d90ff", "#00c61f", "#5701aa", "#001d00", "#9a4600",
    "#959ea5", "#9a425b", "#001f31", "#c8c300", "#ffcfff",
    "#00bd9a", "#3615ff", "#2d2424", "#df57ff", "#bde6bf",
    "#7e4497", "#524f3b", "#d86600", "#647438", "#c17287",
    "#6e7489", "#809c03", "#bd8a64", "#623338", "#cacdda",
    "#6beb82"
  )
}

ccc_ov_palette <- function(n = NULL) {
  cols <- if (!is.null(n) && is.finite(n) && n > length(ccc_ov_palette_28())) {
    ccc_ov_palette_56()
  } else {
    ccc_ov_palette_28()
  }
  if (!is.null(n) && is.finite(n) && n > length(cols)) {
    cols <- rep(cols, length.out = n)
  }
  cols
}

ccc_is_ov_palette <- function(palette) {
  if (is.null(palette) || length(palette) == 0L) {
    return(FALSE)
  }
  if (!is.character(palette) || length(palette) != 1L || is.na(palette)) {
    return(FALSE)
  }
  tolower(palette) %in% c("ov", "omicverse", "palette_28", "palette28")
}

ccc_resolve_category_palette <- function(palette = NULL, palcolor = NULL) {
  if (ccc_is_ov_palette(palette) && is.null(palcolor)) {
    return(list(palette = "Chinese", palcolor = ccc_ov_palette()))
  }
  if (ccc_is_ov_palette(palette)) {
    palette <- "Chinese"
  }
  list(palette = palette %||% "Chinese", palcolor = palcolor)
}

ccc_resolve_value_palette <- function(value_palette, fallback = "Reds") {
  if (
    is.null(value_palette) ||
      length(value_palette) == 0L ||
      is.na(value_palette[1])
  ) {
    return(fallback)
  }
  value_palette <- as.character(value_palette[1])
  switch(tolower(value_palette),
    "rdylbu_r" = "RdYlBu",
    "rdbu_r" = "RdBu",
    "reds_r" = "Reds",
    "blues_r" = "Blues",
    "greens_r" = "Greens",
    "greys_r" = "Greys",
    value_palette
  )
}

patchwork_cc <- function(
  p,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  font.size = NULL
) {
  if (!inherits(p, "patchwork")) {
    return(finalize_cc_plot(
      p,
      title = title,
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size
    ))
  }
  theme_obj <- apply_plot_theme(
    theme_use = theme_use,
    theme_args = theme_args,
    allow_null = TRUE
  )
  if (!is.null(title) || !is.null(subtitle)) {
    p <- p + patchwork::plot_annotation(title = title, subtitle = subtitle)
  }
  if (!is.null(theme_obj)) {
    p <- p &
      theme_obj &
      ggplot2::theme(
        legend.position = legend.position,
        legend.direction = legend.direction
      )
  }
  p
}

plot_cc_list <- function(
  plot_list,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0L) {
    return(NULL)
  }
  if (length(plot_list) == 1L) {
    return(plot_list[[1]])
  }
  if (!isTRUE(combine)) {
    return(plot_list)
  }
  if (
    all(vapply(
      plot_list,
      function(x) inherits(x, c("gg", "ggplot", "patchwork")),
      logical(1)
    ))
  ) {
    return(patchwork::wrap_plots(
      plot_list,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow
    ))
  }
  plot_list
}

simplify_cc_plot_list <- function(plot_list) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 1L) {
    return(plot_list[[1]])
  }
  plot_list
}

ccc_alias_arg <- function(dots, names, current = NULL) {
  if (is.null(dots) || length(dots) == 0L) {
    return(current)
  }
  for (nm in names) {
    if (nm %in% names(dots) && !is.null(dots[[nm]])) {
      return(dots[[nm]])
    }
  }
  current
}

ccc_normalize_use_arg <- function(values) {
  if (is.null(values)) {
    return(NULL)
  }
  values <- as.character(values)
  values <- trimws(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    return(NULL)
  }
  unique(values)
}

ccc_is_provided <- function(value) {
  if (is.null(value)) {
    return(FALSE)
  }
  if (is.character(value)) {
    return(any(!is.na(value) & nzchar(trimws(value))))
  }
  if (is.list(value) || length(value) > 1L) {
    return(length(value) > 0L)
  }
  if (length(value) == 0L) {
    return(FALSE)
  }
  if (is.logical(value) && length(value) == 1L) {
    return(!is.na(value) && isTRUE(value))
  }
  !all(is.na(value))
}

ccc_assert_unsupported <- function(plot_type, ...) {
  args <- list(...)
  if (length(args) == 0L) {
    return(invisible(TRUE))
  }
  unsupported <- names(args)[vapply(args, ccc_is_provided, logical(1))]
  if (length(unsupported) > 0L) {
    formatted <- paste0("{.arg ", unsupported, "}", collapse = ", ")
    log_message(
      "{.val plot_type = {plot_type}} does not support: {formatted}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

ccc_assert_allowed_metrics <- function(
  plot_type,
  color.by = "score",
  allowed_color.by = c("score", "pvalue"),
  value = "sum",
  allowed_value = c("sum", "mean", "max", "count")
) {
  if (!color.by %in% allowed_color.by) {
    log_message(
      "{.val plot_type = {plot_type}} only supports {.arg color.by}: {.val {allowed_color.by}}",
      message_type = "error"
    )
  }
  if (!value %in% allowed_value) {
    log_message(
      "{.val plot_type = {plot_type}} only supports {.arg value}: {.val {allowed_value}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

ccc_mark_significance <- function(df, thresh = 0.05) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  if (!"score" %in% colnames(df)) {
    score_col <- c("prob", "means", "prioritization_score")[
      c("prob", "means", "prioritization_score") %in% colnames(df)
    ][1]
    df$score <- if (!is.null(score_col) && !is.na(score_col)) {
      suppressWarnings(as.numeric(df[[score_col]]))
    } else {
      NA_real_
    }
  }
  if (!"pvalue" %in% colnames(df)) {
    df$pvalue <- if ("pval" %in% colnames(df)) df$pval else NA_real_
  }
  score <- suppressWarnings(as.numeric(df$score))
  pvalue <- suppressWarnings(as.numeric(df$pvalue))
  thresh <- suppressWarnings(as.numeric(thresh)[1])
  if (!is.finite(thresh)) {
    thresh <- 0.05
  }
  has_pvalue <- is.finite(pvalue)
  df$significant <- ifelse(
    has_pvalue,
    pvalue < thresh,
    is.finite(score) & score > 0
  )
  df$significant[is.na(df$significant)] <- FALSE
  if ("significance_basis" %in% colnames(df)) {
    not_tested <- as.character(df$significance_basis) %in% c("not_tested", "none")
    df$significant[not_tested] <- NA
  }
  df$neglog10_pvalue <- ifelse(
    is.finite(pvalue) & pvalue > 0,
    -log10(pmax(pvalue, 1e-300)),
    NA_real_
  )
  df
}

filter_long_df <- function(
  df,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  signaling = NULL,
  pairLR.use = NULL
) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  if (!is.null(sender.use)) {
    df <- df[df$sender %in% sender.use, , drop = FALSE]
  }
  if (!is.null(receiver.use)) {
    df <- df[df$receiver %in% receiver.use, , drop = FALSE]
  }
  if (!is.null(ligand.use)) {
    df <- df[df$ligand %in% ligand.use, , drop = FALSE]
  }
  if (!is.null(receptor.use)) {
    df <- df[df$receptor %in% receptor.use, , drop = FALSE]
  }
  if (!is.null(interaction.use)) {
    keep <- ccc_match_identifiers(
      df = df,
      values = interaction.use,
      columns = c(
        "interaction_name",
        "interaction_label",
        "interaction_display",
        "interaction_name_2",
        "interacting_pair"
      ),
      include_lr = TRUE
    )
    df <- df[keep, , drop = FALSE]
  }
  if (!is.null(pairLR.use)) {
    keep <- ccc_match_identifiers(
      df = df,
      values = pairLR.use,
      columns = c(
        "pair_lr",
        "pairLR",
        "interacting_pair",
        "interaction_name_2",
        "interaction_name",
        "interaction_label",
        "interaction_display"
      ),
      include_lr = TRUE
    )
    df <- df[keep, , drop = FALSE]
  }
  if (!is.null(signaling) && "pathway_name" %in% colnames(df)) {
    signaling_key <- ccc_identifier_key(signaling)
    pathway_key <- ccc_identifier_key(df$pathway_name)
    class_key <- if ("classification" %in% colnames(df)) {
      ccc_identifier_key(df$classification)
    } else {
      pathway_key
    }
    df <- df[
      df$pathway_name %in% signaling |
        pathway_key %in% signaling_key |
        class_key %in% signaling_key, ,
      drop = FALSE
    ]
  }
  df
}

prepare_plot_df <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  out <- df
  out$sender <- as.character(out$sender)
  out$receiver <- as.character(out$receiver)
  out$interaction_name <- as.character(out$interaction_name)
  out$ligand <- as.character(out$ligand)
  out$receptor <- as.character(out$receptor)
  out$pair <- paste(out$sender, out$receiver, sep = " -> ")
  if (!"interaction_label" %in% colnames(out)) {
    out$interaction_label <- out$interaction_name
  }
  out$interaction_label <- as.character(out$interaction_label)
  miss_label <- is.na(out$interaction_label) | !nzchar(out$interaction_label)
  out$interaction_label[miss_label] <- paste(
    out$ligand[miss_label],
    out$receptor[miss_label],
    sep = " - "
  )
  out$interaction_label <- ccc_display_interaction(out$interaction_label)
  if (!"classification" %in% colnames(out)) {
    out$classification <- out$pathway_name
  }
  out$classification <- as.character(out$classification)
  out$classification[is.na(out$classification) | !nzchar(out$classification)] <- "Unclassified"
  out$pathway_name <- out$classification
  if (!"interaction_display" %in% colnames(out)) {
    out$interaction_display <- out$interaction_label
  }
  if (!"ligand_display" %in% colnames(out)) {
    out$ligand_display <- ccc_display_gene(out$ligand)
  }
  if (!"receptor_display" %in% colnames(out)) {
    out$receptor_display <- ccc_display_gene(out$receptor)
  }
  out$pvalue <- suppressWarnings(as.numeric(out$pvalue))
  out$specificity <- ifelse(
    is.finite(out$pvalue) & out$pvalue > 0,
    -log10(out$pvalue),
    NA_real_
  )
  out
}

group_summary <- function(
  df,
  group_cols,
  value_col,
  fun,
  out_col = value_col
) {
  if (is.null(df) || nrow(df) == 0L) {
    out <- as.data.frame(stats::setNames(
      replicate(length(group_cols), character(0), simplify = FALSE),
      group_cols
    ))
    out[[out_col]] <- numeric(0)
    return(out)
  }
  group_df <- df[, group_cols, drop = FALSE]
  for (nm in group_cols) {
    group_df[[nm]] <- as.character(group_df[[nm]])
    group_df[[nm]][is.na(group_df[[nm]])] <- ""
  }
  key <- do.call(paste, c(group_df, sep = "\r"))
  split_idx <- split(seq_len(nrow(df)), key, drop = TRUE)
  pieces <- lapply(split_idx, function(idx) {
    row <- group_df[idx[1], , drop = FALSE]
    row[[out_col]] <- fun(df[[value_col]][idx])
    row
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

pair_plot_df <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  if (!all(c("significant", "neglog10_pvalue") %in% colnames(df))) {
    df <- ccc_mark_significance(df)
  }
  if (!"specificity" %in% colnames(df)) {
    df$specificity <- df$neglog10_pvalue
  }
  if (!"pair" %in% colnames(df)) {
    df$pair <- paste(df$sender, df$receiver, sep = " -> ")
  }
  if (!"interaction_label" %in% colnames(df)) {
    df$interaction_label <- df$interaction_name
  }
  df_use <- df[
    !is.na(df$sender) &
      nzchar(df$sender) &
      !is.na(df$receiver) &
      nzchar(df$receiver), ,
    drop = FALSE
  ]
  if (nrow(df_use) == 0L) {
    return(data.frame())
  }
  pair_df <- aggregate_ccc_long(df_use, backend = "r")
  if (is.null(pair_df) || nrow(pair_df) == 0L) {
    return(data.frame())
  }
  p_df <- group_summary(
    df = df_use,
    group_cols = c("sender", "receiver"),
    value_col = "pvalue",
    out_col = "pvalue",
    fun = function(x) {
      x <- as.numeric(x)
      x <- x[is.finite(x) & x > 0]
      if (length(x) == 0L) {
        return(NA_real_)
      }
      min(x)
    }
  )
  pair_df <- merge(pair_df, p_df, by = c("sender", "receiver"), all.x = TRUE)
  pair_df$pair <- paste(pair_df$sender, pair_df$receiver, sep = " -> ")
  pair_df$specificity <- ifelse(
    is.finite(pair_df$pvalue) & pair_df$pvalue > 0,
    -log10(pair_df$pvalue),
    NA_real_
  )
  pair_df
}

interaction_plot_df <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  if (!all(c("significant", "neglog10_pvalue") %in% colnames(df))) {
    df <- ccc_mark_significance(df)
  }
  if (!"specificity" %in% colnames(df)) {
    df$specificity <- df$neglog10_pvalue
  }
  if (!"pair" %in% colnames(df)) {
    df$pair <- paste(df$sender, df$receiver, sep = " -> ")
  }
  if (!"interaction_label" %in% colnames(df)) {
    df$interaction_label <- df$interaction_name
  }
  df_use <- df[
    !is.na(df$sender) &
      nzchar(df$sender) &
      !is.na(df$receiver) &
      nzchar(df$receiver) &
      !is.na(df$interaction_label) &
      nzchar(df$interaction_label), ,
    drop = FALSE
  ]
  if (nrow(df_use) == 0L) {
    return(data.frame())
  }
  df_use$ligand[is.na(df_use$ligand)] <- ""
  df_use$receptor[is.na(df_use$receptor)] <- ""
  group_cols <- c(
    "sender",
    "receiver",
    "pair",
    "interaction_label",
    "ligand",
    "receptor"
  )
  out <- group_summary(
    df = transform(
      df_use,
      specificity = ifelse(is.finite(specificity), specificity, NA_real_)
    ),
    group_cols = group_cols,
    value_col = "score",
    out_col = "score",
    fun = function(x) {
      x <- as.numeric(x)
      if (all(is.na(x))) {
        return(NA_real_)
      }
      sum(x, na.rm = TRUE)
    }
  )
  s_df <- group_summary(
    df = transform(
      df_use,
      specificity = ifelse(is.finite(specificity), specificity, NA_real_)
    ),
    group_cols = group_cols,
    value_col = "specificity",
    out_col = "specificity",
    fun = function(x) {
      x <- as.numeric(x)
      if (all(is.na(x))) {
        return(NA_real_)
      }
      sum(x, na.rm = TRUE)
    }
  )
  out <- merge(out, s_df, by = group_cols, all = TRUE)
  p_df <- group_summary(
    df = df_use,
    group_cols = group_cols,
    value_col = "pvalue",
    out_col = "pvalue",
    fun = function(x) {
      x <- as.numeric(x)
      x <- x[is.finite(x) & x > 0]
      if (length(x) == 0L) {
        return(NA_real_)
      }
      min(x)
    }
  )
  out <- merge(
    out,
    p_df,
    by = group_cols,
    all.x = TRUE
  )
  count_df <- group_summary(
    df = df_use,
    group_cols = group_cols,
    value_col = "significant",
    out_col = "count",
    fun = function(x) {
      sum(as.numeric(x), na.rm = TRUE)
    }
  )
  out <- merge(
    out,
    count_df,
    by = group_cols,
    all.x = TRUE
  )
  out$count <- suppressWarnings(as.numeric(out$count))
  out$count[!is.finite(out$count)] <- 0
  out$specificity <- ifelse(
    is.finite(out$pvalue) & out$pvalue > 0,
    -log10(out$pvalue),
    out$specificity
  )
  out
}

top_pairs <- function(pair_df, top_n = 20, value_col = "sum") {
  if (
    is.null(pair_df) || nrow(pair_df) == 0L || !is.numeric(top_n) || top_n <= 0L
  ) {
    return(pair_df)
  }
  ord <- order(pair_df[[value_col]], decreasing = TRUE, na.last = TRUE)
  pair_df[utils::head(ord, top_n), , drop = FALSE]
}

top_interactions <- function(
  interaction_df,
  top_n = 20,
  value_col = "score"
) {
  if (
    is.null(interaction_df) ||
      nrow(interaction_df) == 0L ||
      !is.numeric(top_n) ||
      top_n <= 0L
  ) {
    return(interaction_df)
  }
  if (!value_col %in% colnames(interaction_df)) {
    value_col <- c("score", "sum", "mean", "max", "count")[
      c("score", "sum", "mean", "max", "count") %in% colnames(interaction_df)
    ][1]
  }
  if (is.null(value_col) || is.na(value_col)) {
    return(interaction_df)
  }
  label_col <- c("interaction_label", "interaction_display", "interaction_name")[
    c("interaction_label", "interaction_display", "interaction_name") %in% colnames(interaction_df)
  ][1]
  if (is.null(label_col) || is.na(label_col)) {
    ord <- order(interaction_df[[value_col]], decreasing = TRUE, na.last = TRUE)
    return(interaction_df[utils::head(ord, top_n), , drop = FALSE])
  }
  labels <- as.character(interaction_df[[label_col]])
  labels[is.na(labels) | !nzchar(labels)] <- "Unclassified"
  vals <- suppressWarnings(as.numeric(interaction_df[[value_col]]))
  vals[!is.finite(vals)] <- 0
  ranked <- tapply(vals, labels, sum, na.rm = TRUE)
  ranked <- sort(ranked, decreasing = TRUE)
  keep <- names(ranked)[seq_len(min(as.integer(top_n), length(ranked)))]
  out <- interaction_df[labels %in% keep, , drop = FALSE]
  out$.ccc_rank <- match(as.character(out[[label_col]]), keep)
  out <- out[order(out$.ccc_rank, -vals[labels %in% keep], na.last = TRUE), , drop = FALSE]
  out$.ccc_rank <- NULL
  rownames(out) <- NULL
  out
}

ccc_palettes <- function(
  palette = "Chinese",
  palcolor = NULL,
  value_palette = NULL,
  value_palcolor = NULL,
  cell_palette = NULL,
  cell_palcolor = NULL,
  link_palette = NULL,
  link_palcolor = NULL
) {
  base_cfg <- ccc_resolve_category_palette(palette = palette, palcolor = palcolor)
  cell_cfg <- ccc_resolve_category_palette(
    palette = cell_palette %||% palette,
    palcolor = cell_palcolor %||% palcolor
  )
  link_cfg <- ccc_resolve_category_palette(
    palette = link_palette %||% palette,
    palcolor = link_palcolor %||% palcolor
  )
  list(
    value_palette = ccc_resolve_value_palette(value_palette, fallback = "RdBu"),
    value_palcolor = value_palcolor,
    cell_palette = cell_cfg$palette %||% base_cfg$palette,
    cell_palcolor = cell_cfg$palcolor %||% base_cfg$palcolor,
    link_palette = link_cfg$palette %||% base_cfg$palette,
    link_palcolor = link_cfg$palcolor %||% base_cfg$palcolor
  )
}

ccc_assign_plot_score <- function(df, value = "score") {
  if (is.null(df) || nrow(df) == 0L) {
    return(df)
  }
  df <- ccc_mark_significance(df)
  if (!"pvalue" %in% colnames(df)) {
    df$pvalue <- if ("pval" %in% colnames(df)) df$pval else NA_real_
  }
  if (identical(value, "count")) {
    df$score <- as.numeric(df$significant)
    return(df)
  }
  score_col <- if (identical(value, "weight")) {
    NULL
  } else {
    value
  }
  if (is.null(score_col) || !score_col %in% colnames(df)) {
    score_col <- c("score", "prob", "means")[
      c("score", "prob", "means") %in% colnames(df)
    ][1]
  }
  if (is.null(score_col) || is.na(score_col)) {
    score_col <- "score"
    if (!"score" %in% colnames(df)) {
      df$score <- NA_real_
      return(df)
    }
  }
  df$score <- df[[score_col]]
  df
}

ccc_filter_table_context <- function(
  df,
  resource = NULL,
  condition = NULL,
  sample = NULL
) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    return(df)
  }
  if (!is.null(resource)) {
    available_resource <- if ("resource" %in% colnames(df)) {
      unique(as.character(df$resource))
    } else {
      character(0)
    }
    available_resource <- available_resource[
      !is.na(available_resource) & nzchar(available_resource)
    ]
    if (length(available_resource) == 0L) {
      log_message(
        "The selected CCC result does not provide {.arg resource} provenance",
        message_type = "error"
      )
    }
    df <- df[as.character(df$resource) %in% as.character(resource), , drop = FALSE]
  }
  if (!is.null(condition)) {
    available_condition <- if ("condition" %in% colnames(df)) {
      unique(as.character(df$condition))
    } else {
      character(0)
    }
    available_condition <- available_condition[
      !is.na(available_condition) & nzchar(available_condition)
    ]
    if (length(available_condition) == 0L) {
      log_message(
        "The selected CCC result does not provide {.arg condition} provenance",
        message_type = "error"
      )
    }
    df <- df[
      as.character(df$condition) %in% as.character(condition), ,
      drop = FALSE
    ]
  }
  if (!is.null(sample)) {
    sample_col <- c("sample", "context", "dataset")
    sample_col <- sample_col[sample_col %in% colnames(df)][1]
    available_sample <- if (!is.na(sample_col)) {
      unique(as.character(df[[sample_col]]))
    } else {
      character(0)
    }
    available_sample <- available_sample[
      !is.na(available_sample) & nzchar(available_sample)
    ]
    if (is.na(sample_col) || length(available_sample) == 0L) {
      log_message(
        "The selected CCC result does not provide sample or context provenance",
        message_type = "error"
      )
    }
    df <- df[as.character(df[[sample_col]]) %in% as.character(sample), , drop = FALSE]
  }
  df
}

ccc_prepare_filtered_object <- function(
  srt,
  method,
  resource = NULL,
  condition = NULL,
  sample = NULL
) {
  if (is.null(resource) && is.null(condition) && is.null(sample)) {
    return(srt)
  }
  method <- normalize_ccc_method(method)
  targets <- unique(c("CCC", if (!identical(method, "CCC")) method))
  targets <- targets[targets %in% names(srt@tools)]
  for (target in targets) {
    bundle <- srt@tools[[target]]
    for (field in c("long_table", "primary_table", "consensus_table")) {
      if (is.data.frame(bundle[[field]])) {
        bundle[[field]] <- ccc_filter_table_context(
          bundle[[field]],
          resource = resource,
          condition = condition,
          sample = sample
        )
      }
    }
    srt@tools[[target]] <- bundle
  }
  srt
}

ccc_combine_methods <- function(
  df,
  mode = c("separate", "support", "rank", "legacy")
) {
  mode <- match.arg(mode)
  df <- ccc_semantic_long_table(df)
  if (nrow(df) == 0L || identical(mode, "separate")) {
    return(df)
  }
  if (identical(mode, "legacy")) {
    log_message(
      "{.val combine_methods = 'legacy'} sums backend scores on incompatible scales and is deprecated",
      message_type = "warning"
    )
    return(df)
  }
  required <- c("sender", "receiver", "ligand", "receptor", "method")
  if (!all(required %in% colnames(df))) {
    log_message("Unified CCC data are missing method or interaction identifiers", message_type = "error")
  }
  has_interaction <- !is.na(df$ligand) & nzchar(trimws(as.character(df$ligand))) &
    !is.na(df$receptor) & nzchar(trimws(as.character(df$receptor)))
  df <- df[has_interaction, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(df)
  }
  methods <- unique(as.character(df$method))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  df$.ccc_method_percentile <- 1
  for (method_i in methods) {
    idx <- which(as.character(df$method) == method_i)
    raw_rank <- suppressWarnings(as.numeric(df$priority_rank[idx]))
    finite <- is.finite(raw_rank)
    if (any(finite)) {
      percentile <- (rank(raw_rank[finite], ties.method = "average") - 1) /
        max(1, sum(finite) - 1)
      df$.ccc_method_percentile[idx[finite]] <- percentile
    }
  }
  key_cols <- c("sender", "receiver", "ligand", "receptor")
  key <- do.call(paste, c(lapply(df[key_cols], as.character), sep = "\r"))
  split_idx <- split(seq_len(nrow(df)), key, drop = TRUE)
  rows <- lapply(split_idx, function(idx) {
    x <- df[idx, , drop = FALSE]
    row <- x[1, , drop = FALSE]
    support_methods <- unique(as.character(x$method))
    support_methods <- support_methods[!is.na(support_methods) & nzchar(support_methods)]
    row$method <- "CCC"
    row$support_methods <- paste(sort(support_methods), collapse = ";")
    row$support_count <- length(support_methods)
    row$support_fraction <- length(support_methods) / max(1L, length(methods))
    if (identical(mode, "support")) {
      row$score <- row$support_count
      row$priority_score <- row$support_fraction
      row$priority_rank <- 1 - row$support_fraction
      row$score_type <- "method_support_count"
      row$support_type <- "cross_method_support"
    } else {
      rank_values <- suppressWarnings(as.numeric(x$.ccc_method_percentile))
      ranks <- tapply(
        rank_values,
        as.character(x$method),
        function(value) {
          value <- value[is.finite(value)]
          if (length(value)) mean(value) else 1
        }
      )
      ranks <- as.numeric(ranks)
      missing_methods <- max(0L, length(methods) - length(ranks))
      ranks <- c(ranks, rep(1, missing_methods))
      mean_rank <- if (length(ranks)) mean(ranks) else 1
      row$score <- 1 - mean_rank
      row$priority_score <- 1 - mean_rank
      row$priority_rank <- mean_rank
      row$score_type <- "mean_within_method_percentile"
      row$support_type <- "scop_visualization_consensus"
    }
    row$.ccc_method_percentile <- NULL
    row$pvalue <- NA_real_
    row$pvalue_type <- "not_available"
    row$significant <- row$support_count > 0L
    row
  })
  ccc_bind_long_tables(rows)
}

ccc_prepare_combined_object <- function(
  srt,
  method,
  combine_methods = c("separate", "support", "rank", "legacy")
) {
  combine_methods <- match.arg(combine_methods)
  method <- normalize_ccc_method(method)
  if (!identical(method, "CCC") || identical(combine_methods, "separate")) {
    return(srt)
  }
  bundle <- srt@tools[["CCC"]]
  if (is.null(bundle)) {
    bundle <- ccc_build_unified_bundle(srt)
  }
  liana_methods <- tolower(as.character(
    srt@tools[["LIANA"]]$parameters$method %||% character(0)
  ))
  if (
    combine_methods %in% c("support", "rank") &&
      all(c("CellphoneDB", "LIANA") %in% (bundle$methods %||% character(0))) &&
      "cellphonedb" %in% liana_methods
  ) {
    log_message(
      paste0(
        "Standalone CellphoneDB and the LIANA consensus are not independent ",
        "evidence because LIANA includes its CellPhoneDB scoring method; ",
        "interpret this as backend agreement, not independent validation"
      ),
      message_type = "warning"
    )
  }
  bundle$long_table <- ccc_combine_methods(bundle$long_table, mode = combine_methods)
  bundle$pair_table <- aggregate_ccc_long(bundle$long_table, backend = "r")
  bundle$metadata$combine_methods <- combine_methods
  srt@tools[["CCC"]] <- bundle
  srt
}

ccc_plot_methods_separately <- function(call, srt, env) {
  unified <- ccc_semantic_long_table(srt@tools[["CCC"]]$long_table)
  methods <- unique(as.character(unified$method))
  methods <- methods[!is.na(methods) & nzchar(methods) & methods != "CCC"]
  plots <- lapply(methods, function(method) {
    method_table <- unified[as.character(unified$method) == method, , drop = FALSE]
    bundle <- srt@tools[[method]] %||% list(method = method)
    bundle$long_table <- method_table
    bundle$primary_table <- method_table
    bundle$pair_table <- aggregate_ccc_long(method_table, backend = "r")
    method_object <- srt
    method_object@tools[[method]] <- bundle
    next_call <- call
    next_call$srt <- NULL
    next_call$object <- method_object
    next_call$method <- method
    next_call$combine_methods <- "legacy"
    next_call$resource <- NULL
    next_call$condition <- NULL
    next_call$sample <- NULL
    eval(next_call, envir = env)
  })
  names(plots) <- methods
  plots
}


ccc_matrix_heatmap_plot <- function(
  mat,
  value_label = "Value",
  row_title = "Row",
  column_title = "Column",
  row_annotation_name = row_title,
  column_annotation_name = column_title,
  top_anno = "bar",
  right_anno = "cell",
  left_anno = "bar",
  bottom_anno = "cell",
  bar_value = "sum",
  add_text = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  x_text_angle = 90,
  border = TRUE,
  width = NULL,
  height = NULL,
  units = "inch",
  title = NULL,
  subtitle = NULL,
  value_palette = "RdBu",
  value_palcolor = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  symmetric = FALSE
) {
  check_r("circlize", verbose = FALSE)
  top_anno <- ccc_match_side_anno(top_anno)
  right_anno <- ccc_match_side_anno(right_anno)
  left_anno <- ccc_match_side_anno(left_anno)
  bottom_anno <- ccc_match_side_anno(bottom_anno)
  bar_value <- ccc_match_bar_value(bar_value)

  if (is.null(mat) || length(mat) == 0L || nrow(mat) == 0L || ncol(mat) == 0L) {
    log_message(
      "The heatmap matrix is empty after applying the current filters",
      message_type = "error"
    )
  }

  top_values <- ccc_matrix_group_values_from_matrix(mat, margin = "col")
  right_values <- ccc_matrix_group_values_from_matrix(mat, margin = "row")
  top_bar_vec <- if (ccc_side_has_bar(top_anno, bottom_anno)) {
    ccc_matrix_bar_stats_from_values(top_values, metrics = bar_value)
  } else {
    NULL
  }
  top_box_vec <- if (ccc_side_has_distribution(top_anno, bottom_anno)) {
    top_values
  } else {
    NULL
  }
  top_summary_vec <- if (ccc_side_has_summary(top_anno, bottom_anno)) {
    ccc_matrix_summary_from_values(top_values, metric = "mean")
  } else {
    NULL
  }
  right_bar_vec <- if (ccc_side_has_bar(left_anno, right_anno)) {
    ccc_matrix_bar_stats_from_values(right_values, metrics = bar_value)
  } else {
    NULL
  }
  right_box_vec <- if (ccc_side_has_distribution(left_anno, right_anno)) {
    right_values
  } else {
    NULL
  }
  right_summary_vec <- if (ccc_side_has_summary(left_anno, right_anno)) {
    ccc_matrix_summary_from_values(right_values, metric = "mean")
  } else {
    NULL
  }

  mat_vals <- mat[is.finite(mat)]
  if (length(mat_vals) == 0L) {
    mat_vals <- c(0, 1)
  }
  if (isTRUE(symmetric)) {
    max_abs <- max(abs(mat_vals), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- 1
    }
    val_range <- c(-max_abs, max_abs)
  } else {
    val_range <- range(mat_vals, na.rm = TRUE)
    if (val_range[1] == val_range[2]) {
      val_range <- val_range + c(-0.5, 0.5)
    }
  }
  fill_cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 100
  )
  col_fun <- circlize::colorRamp2(
    breaks = seq(val_range[1], val_range[2], length.out = length(fill_cols)),
    colors = fill_cols
  )

  cell_fun <- if (isTRUE(add_text)) {
    mat_local <- mat
    function(j, i, x, y, width, height, fill) {
      v <- mat_local[i, j]
      if (is.finite(v)) {
        grid::grid.text(
          sprintf("%.2g", v),
          x,
          y,
          gp = grid::gpar(
            fontsize = font.size * 0.7,
            col = if (!is.na(fill) && grDevices::col2rgb(fill)[1] < 128) {
              "white"
            } else {
              "grey20"
            }
          )
        )
      }
    }
  } else {
    NULL
  }

  body_size <- ccc_heatmap_body_size(
    nrow = nrow(mat),
    ncol = ncol(mat),
    width = width,
    height = height,
    units = units
  )

  col_ann <- if (!is.null(colnames(mat))) {
    col_nms <- colnames(mat)
    col_cols <- palette_colors(
      col_nms,
      palette = cell_palette,
      palcolor = cell_palcolor
    )
    ComplexHeatmap::HeatmapAnnotation(
      df = stats::setNames(
        data.frame(col_nms, row.names = col_nms, stringsAsFactors = FALSE),
        column_annotation_name
      ),
      col = stats::setNames(
        list(stats::setNames(col_cols, col_nms)),
        column_annotation_name
      ),
      annotation_name_gp = grid::gpar(fontsize = font.size * 0.8),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      which = "column",
      border = border,
      gp = grid::gpar(col = if (isTRUE(border)) "white" else NA)
    )
  } else {
    NULL
  }

  row_ann <- if (!is.null(rownames(mat))) {
    row_nms <- rownames(mat)
    row_cols <- palette_colors(
      row_nms,
      palette = cell_palette,
      palcolor = cell_palcolor
    )
    ComplexHeatmap::rowAnnotation(
      df = stats::setNames(
        data.frame(row_nms, row.names = row_nms, stringsAsFactors = FALSE),
        row_annotation_name
      ),
      col = stats::setNames(
        list(stats::setNames(row_cols, row_nms)),
        row_annotation_name
      ),
      annotation_name_gp = grid::gpar(fontsize = font.size * 0.8),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      border = border,
      gp = grid::gpar(col = if (isTRUE(border)) "white" else NA)
    )
  } else {
    NULL
  }

  top_names <- colnames(mat)
  right_names <- rownames(mat)
  top_cols <- if (!is.null(top_names)) {
    palette_colors(top_names, palette = cell_palette, palcolor = cell_palcolor)
  } else {
    NULL
  }
  right_cols <- if (!is.null(right_names)) {
    palette_colors(
      right_names,
      palette = cell_palette,
      palcolor = cell_palcolor
    )
  } else {
    NULL
  }

  top_bar_ann <- ccc_build_bar_annotation(
    stats_list = top_bar_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    reverse = FALSE
  )
  left_bar_ann <- ccc_build_bar_annotation(
    stats_list = right_bar_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    reverse = TRUE
  )
  right_bar_ann <- ccc_build_bar_annotation(
    stats_list = right_bar_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    reverse = FALSE
  )
  top_box_ann <- ccc_build_distribution_annotation(
    values = top_box_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    type = "box",
    reverse = FALSE
  )
  left_box_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "box",
    reverse = TRUE
  )
  right_box_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "box",
    reverse = FALSE
  )
  top_hist_ann <- ccc_build_distribution_annotation(
    values = top_box_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    type = "histogram",
    reverse = FALSE
  )
  left_hist_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "histogram",
    reverse = TRUE
  )
  right_hist_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "histogram",
    reverse = FALSE
  )
  top_density_ann <- ccc_build_distribution_annotation(
    values = top_box_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    type = "density",
    reverse = FALSE
  )
  left_density_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "density",
    reverse = TRUE
  )
  right_density_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "density",
    reverse = FALSE
  )
  top_violin_ann <- ccc_build_distribution_annotation(
    values = top_box_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    type = "violin",
    reverse = FALSE
  )
  left_violin_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "violin",
    reverse = TRUE
  )
  right_violin_ann <- ccc_build_distribution_annotation(
    values = right_box_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "violin",
    reverse = FALSE
  )
  top_point_ann <- ccc_build_summary_annotation(
    values = top_summary_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    type = "point",
    reverse = FALSE
  )
  left_point_ann <- ccc_build_summary_annotation(
    values = right_summary_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "point",
    reverse = TRUE
  )
  right_point_ann <- ccc_build_summary_annotation(
    values = right_summary_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "point",
    reverse = FALSE
  )
  top_line_ann <- ccc_build_summary_annotation(
    values = top_summary_vec,
    fill_cols = top_cols,
    which = "column",
    font.size = font.size,
    border = border,
    type = "line",
    reverse = FALSE
  )
  left_line_ann <- ccc_build_summary_annotation(
    values = right_summary_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "line",
    reverse = TRUE
  )
  right_line_ann <- ccc_build_summary_annotation(
    values = right_summary_vec,
    fill_cols = right_cols,
    which = "row",
    font.size = font.size,
    border = border,
    type = "line",
    reverse = FALSE
  )

  top_combined <- ccc_combine_side_annotations(
    side_anno = top_anno,
    annotations = list(
      bar = top_bar_ann,
      box = top_box_ann,
      point = top_point_ann,
      line = top_line_ann,
      histogram = top_hist_ann,
      density = top_density_ann,
      violin = top_violin_ann,
      cell = col_ann
    )
  )
  bottom_combined <- ccc_combine_side_annotations(
    side_anno = bottom_anno,
    annotations = list(
      bar = top_bar_ann,
      box = top_box_ann,
      point = top_point_ann,
      line = top_line_ann,
      histogram = top_hist_ann,
      density = top_density_ann,
      violin = top_violin_ann,
      cell = col_ann
    )
  )
  left_combined <- ccc_combine_side_annotations(
    side_anno = left_anno,
    annotations = list(
      bar = left_bar_ann,
      box = left_box_ann,
      point = left_point_ann,
      line = left_line_ann,
      histogram = left_hist_ann,
      density = left_density_ann,
      violin = left_violin_ann,
      cell = row_ann
    )
  )
  right_combined <- ccc_combine_side_annotations(
    side_anno = right_anno,
    annotations = list(
      bar = right_bar_ann,
      box = right_box_ann,
      point = right_point_ann,
      line = right_line_ann,
      histogram = right_hist_ann,
      density = right_density_ann,
      violin = right_violin_ann,
      cell = row_ann
    )
  )

  ht <- ComplexHeatmap::Heatmap(
    matrix = mat,
    col = col_fun,
    name = value_label,
    cell_fun = cell_fun,
    top_annotation = top_combined,
    bottom_annotation = bottom_combined,
    left_annotation = left_combined,
    right_annotation = right_combined,
    row_title = row_title,
    column_title = column_title,
    row_title_gp = grid::gpar(fontsize = font.size),
    column_title_gp = grid::gpar(fontsize = font.size),
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_names_gp = grid::gpar(fontsize = font.size),
    column_names_gp = grid::gpar(fontsize = font.size),
    column_names_rot = x_text_angle,
    border = border,
    rect_gp = grid::gpar(
      col = if (isTRUE(border)) "white" else NA,
      lwd = 1
    ),
    use_raster = FALSE,
    width = body_size$width,
    height = body_size$height,
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = font.size * 0.9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = font.size * 0.8),
      border = "black",
      grid_width = grid::unit(4, "mm")
    )
  )

  g_tree <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht,
      padding = grid::unit(c(4, 4, 12, 4), "mm")
    ),
    wrap = TRUE,
    wrap.grobs = TRUE
  )

  if (isTRUE(body_size$fixed)) {
    overall_size <- ccc_heatmap_capture_size(
      body_width = body_size$width_num,
      body_height = body_size$height_num,
      top_tracks = ccc_count_side_tracks(
        top_anno,
        track_counts = list(
          bar = if (is.null(top_bar_vec)) 0L else length(top_bar_vec),
          box = if (is.null(top_box_vec)) 0L else 1L,
          point = if (is.null(top_summary_vec)) 0L else 1L,
          line = if (is.null(top_summary_vec)) 0L else 1L,
          histogram = if (is.null(top_box_vec)) 0L else 1L,
          density = if (is.null(top_box_vec)) 0L else 1L,
          violin = if (is.null(top_box_vec)) 0L else 1L,
          cell = if (is.null(col_ann)) 0L else 1L
        )
      ),
      bottom_tracks = ccc_count_side_tracks(
        bottom_anno,
        track_counts = list(
          bar = if (is.null(top_bar_vec)) 0L else length(top_bar_vec),
          box = if (is.null(top_box_vec)) 0L else 1L,
          point = if (is.null(top_summary_vec)) 0L else 1L,
          line = if (is.null(top_summary_vec)) 0L else 1L,
          histogram = if (is.null(top_box_vec)) 0L else 1L,
          density = if (is.null(top_box_vec)) 0L else 1L,
          violin = if (is.null(top_box_vec)) 0L else 1L,
          cell = if (is.null(col_ann)) 0L else 1L
        )
      ),
      left_tracks = ccc_count_side_tracks(
        left_anno,
        track_counts = list(
          bar = if (is.null(right_bar_vec)) 0L else length(right_bar_vec),
          box = if (is.null(right_box_vec)) 0L else 1L,
          point = if (is.null(right_summary_vec)) 0L else 1L,
          line = if (is.null(right_summary_vec)) 0L else 1L,
          histogram = if (is.null(right_box_vec)) 0L else 1L,
          density = if (is.null(right_box_vec)) 0L else 1L,
          violin = if (is.null(right_box_vec)) 0L else 1L,
          cell = if (is.null(row_ann)) 0L else 1L
        )
      ),
      right_tracks = ccc_count_side_tracks(
        right_anno,
        track_counts = list(
          bar = if (is.null(right_bar_vec)) 0L else length(right_bar_vec),
          box = if (is.null(right_box_vec)) 0L else 1L,
          point = if (is.null(right_summary_vec)) 0L else 1L,
          line = if (is.null(right_summary_vec)) 0L else 1L,
          histogram = if (is.null(right_box_vec)) 0L else 1L,
          density = if (is.null(right_box_vec)) 0L else 1L,
          violin = if (is.null(right_box_vec)) 0L else 1L,
          cell = if (is.null(row_ann)) 0L else 1L
        )
      ),
      legend.position = legend.position,
      has_title = !is.null(title) || !is.null(subtitle)
    )
    p <- panel_fix_overall(
      g_tree,
      width = overall_size$width,
      height = overall_size$height,
      units = units
    )
  } else {
    p <- patchwork::wrap_plots(g_tree)
  }

  if (!is.null(title) || !is.null(subtitle)) {
    p <- p +
      ggplot2::labs(title = title, subtitle = subtitle) +
      apply_plot_theme(theme_use, theme_args) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = font.size * 1.2, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = font.size, hjust = 0.5)
      )
  }

  p
}


ccc_cellchat_role_heatmap_plot <- function(
  srt,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  signaling = NULL,
  pattern = "outgoing",
  top_anno = "bar",
  right_anno = "cell",
  left_anno = "bar",
  bottom_anno = "cell",
  bar_value = "sum",
  add_text = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  x_text_angle = 90,
  border = TRUE,
  width = NULL,
  height = NULL,
  units = "inch",
  title = NULL,
  subtitle = NULL,
  value_palette = "RdBu",
  value_palcolor = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  add_text_eff <- if (is.null(add_text)) FALSE else isTRUE(add_text)
  build_role_heatmap <- function(
    mat,
    plot_title,
    plot_subtitle = NULL,
    symmetric = FALSE
  ) {
    ccc_matrix_heatmap_plot(
      mat = mat,
      value_label = if (isTRUE(symmetric)) {
        "Difference"
      } else {
        "Relative strength"
      },
      row_title = "Pathway",
      column_title = "Cell type",
      row_annotation_name = "Pathway",
      column_annotation_name = "Cell type",
      top_anno = top_anno,
      right_anno = right_anno,
      left_anno = left_anno,
      bottom_anno = bottom_anno,
      bar_value = bar_value,
      add_text = add_text_eff,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      x_text_angle = x_text_angle,
      border = border,
      width = width,
      height = height,
      units = units,
      title = plot_title,
      subtitle = plot_subtitle,
      value_palette = value_palette,
      value_palcolor = value_palcolor,
      cell_palette = cell_palette,
      cell_palcolor = cell_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args,
      symmetric = symmetric
    )
  }
  if (isTRUE(use_cc_single_condition(srt, condition = condition))) {
    info <- get_dataset_object(srt, condition = condition, dataset = dataset)
    role <- ccc_cellchat_role_matrix(
      object = info$object,
      signaling = signaling,
      pattern = pattern,
      scale_rows = TRUE
    )
    return(build_role_heatmap(
      mat = role$scaled,
      plot_title = title %||% paste0(info$label, ": ", pattern),
      plot_subtitle = subtitle
    ))
  }

  cc_cmp <- ccc_cellchat_heatmap_comparison(
    srt = srt,
    condition = condition,
    comparison = comparison
  )
  object_names <- cc_cmp$object_names
  object_list <- cc_cmp$object_list
  if (is.null(signaling)) {
    signaling <- Reduce(
      union,
      lapply(object_list, function(obj) obj@netP$pathways)
    )
  }
  plots <- lapply(seq_along(object_list), function(i) {
    role <- ccc_cellchat_role_matrix(
      object = object_list[[i]],
      signaling = signaling,
      pattern = pattern,
      scale_rows = TRUE
    )
    build_role_heatmap(
      mat = role$scaled,
      plot_title = object_names[i]
    )
  })
  out <- Reduce(`+`, plots)
  if (!is.null(title) || !is.null(subtitle)) {
    out <- out +
      patchwork::plot_annotation(
        title = title %||% paste0("CellChat role heatmap: ", pattern),
        subtitle = subtitle
      ) &
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = font.size * 1.2, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = font.size, hjust = 0.5)
      )
  }
  out
}


ccc_cellchat_diff_heatmap_plot <- function(
  srt,
  condition = NULL,
  comparison = c(1, 2),
  signaling = NULL,
  pattern = "outgoing",
  top_anno = "bar",
  right_anno = "cell",
  left_anno = "bar",
  bottom_anno = "cell",
  bar_value = "sum",
  add_text = NULL,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  x_text_angle = 90,
  border = TRUE,
  width = NULL,
  height = NULL,
  units = "inch",
  title = NULL,
  subtitle = NULL,
  value_palette = "RdBu",
  value_palcolor = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  cc_cmp <- ccc_cellchat_heatmap_comparison(
    srt = srt,
    condition = condition,
    comparison = comparison,
    min_n = 2L
  )
  object_names <- cc_cmp$object_names[seq_len(2)]
  object_list <- cc_cmp$object_list[seq_len(2)]
  if (is.null(signaling)) {
    signaling <- Reduce(
      union,
      lapply(object_list, function(obj) obj@netP$pathways)
    )
  }
  mats <- lapply(object_list, function(obj) {
    ccc_cellchat_role_matrix(
      object = obj,
      signaling = signaling,
      pattern = pattern,
      scale_rows = TRUE
    )$scaled
  })
  diff_mat <- mats[[2]] - mats[[1]]
  diff_mat[!is.finite(diff_mat)] <- NA_real_
  ccc_matrix_heatmap_plot(
    mat = diff_mat,
    value_label = "Difference",
    row_title = "Pathway",
    column_title = "Cell type",
    row_annotation_name = "Pathway",
    column_annotation_name = "Cell type",
    top_anno = top_anno,
    right_anno = right_anno,
    left_anno = left_anno,
    bottom_anno = bottom_anno,
    bar_value = bar_value,
    add_text = if (is.null(add_text)) FALSE else isTRUE(add_text),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    x_text_angle = x_text_angle,
    border = border,
    width = width,
    height = height,
    units = units,
    title = title %||%
      paste0(object_names[2], " vs ", object_names[1], ": ", pattern),
    subtitle = subtitle,
    value_palette = value_palette,
    value_palcolor = value_palcolor,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    legend.position = legend.position,
    legend.direction = legend.direction,
    font.size = font.size,
    theme_use = theme_use,
    theme_args = theme_args,
    symmetric = TRUE
  )
}


ccc_match_side_anno <- function(side_anno) {
  if (is.null(side_anno)) {
    return(character(0))
  }
  side_anno <- unique(as.character(side_anno))
  side_anno <- side_anno[!is.na(side_anno)]
  side_anno <- side_anno[side_anno != ""]
  if (length(side_anno) == 0L) {
    return(character(0))
  }
  if ("none" %in% tolower(side_anno)) {
    return(character(0))
  }
  allowed <- c(
    "bar",
    "box",
    "point",
    "line",
    "histogram",
    "density",
    "violin",
    "cell"
  )
  bad <- setdiff(side_anno, allowed)
  if (length(bad) > 0L) {
    log_message(
      "{.arg side_anno} values must be drawn from {.val {allowed}}. Invalid values: {.val {bad}}",
      message_type = "error"
    )
  }
  side_anno
}


ccc_update_side_anno_legacy <- function(side_anno, show_bar = TRUE) {
  side_anno <- ccc_match_side_anno(side_anno)
  if (isTRUE(show_bar)) {
    return(unique(c("bar", side_anno)))
  }
  setdiff(side_anno, "bar")
}


ccc_side_has_bar <- function(...) {
  annos <- unlist(list(...), use.names = FALSE)
  "bar" %in% annos
}


ccc_side_has_distribution <- function(...) {
  annos <- unlist(list(...), use.names = FALSE)
  any(c("box", "histogram", "density", "violin") %in% annos)
}


ccc_side_has_summary <- function(...) {
  annos <- unlist(list(...), use.names = FALSE)
  any(c("point", "line") %in% annos)
}


ccc_combine_side_annotations <- function(
  side_anno,
  annotations = list()
) {
  side_anno <- ccc_match_side_anno(side_anno)
  ann_list <- list()
  for (nm in side_anno) {
    ann_cur <- annotations[[nm]]
    if (!is.null(ann_cur)) {
      ann_list[[length(ann_list) + 1L]] <- ann_cur
    }
  }
  if (length(ann_list) == 0L) {
    return(NULL)
  }
  out <- ann_list[[1]]
  if (length(ann_list) > 1L) {
    for (i in 2:length(ann_list)) {
      out <- c(out, ann_list[[i]])
    }
  }
  out
}


ccc_count_side_tracks <- function(side_anno, track_counts = list()) {
  side_anno <- ccc_match_side_anno(side_anno)
  n <- 0L
  for (nm in side_anno) {
    cur <- track_counts[[nm]] %||% 0L
    n <- n + as.integer(cur)[1]
  }
  n
}


ccc_match_bar_value <- function(bar_value) {
  allowed <- c("count", "sum", "mean", "max")
  out <- unique(as.character(bar_value %||% "sum"))
  bad <- setdiff(out, allowed)
  if (length(bad) > 0L) {
    log_message(
      "{.arg bar_value} must be one or more of {.val {allowed}}. Invalid values: {.val {bad}}",
      message_type = "error"
    )
  }
  out
}


ccc_bar_fun <- function(bar_value) {
  if (identical(bar_value, "count")) {
    return(function(x) {
      sum(is.finite(as.numeric(x)) & as.numeric(x) > 0, na.rm = TRUE)
    })
  }
  if (identical(bar_value, "mean")) {
    return(function(x) {
      x <- as.numeric(x)
      if (all(is.na(x))) {
        return(NA_real_)
      }
      mean(x, na.rm = TRUE)
    })
  }
  if (identical(bar_value, "max")) {
    return(function(x) {
      x <- as.numeric(x)
      if (all(is.na(x))) {
        return(NA_real_)
      }
      max(x, na.rm = TRUE)
    })
  }
  function(x) {
    sum(as.numeric(x), na.rm = TRUE)
  }
}


ccc_collect_bar_stats <- function(values, groups, ordered_levels, metrics) {
  out <- lapply(metrics, function(metric) {
    vec <- tapply(values, groups, ccc_bar_fun(metric))
    vec <- as.numeric(vec[ordered_levels])
    stats::setNames(vec, ordered_levels)
  })
  names(out) <- metrics
  out
}


ccc_collect_group_values <- function(values, groups, ordered_levels) {
  values <- as.numeric(values)
  groups <- as.character(groups)
  split_vals <- split(values, groups)
  out <- lapply(ordered_levels, function(level) {
    vec <- split_vals[[level]]
    vec <- vec[is.finite(vec)]
    if (length(vec) == 0L) {
      return(NA_real_)
    }
    vec
  })
  stats::setNames(out, ordered_levels)
}


ccc_collect_group_summary <- function(
  values,
  groups,
  ordered_levels,
  metric = "mean"
) {
  vec <- tapply(values, groups, ccc_bar_fun(metric))
  vec <- as.numeric(vec[ordered_levels])
  stats::setNames(vec, ordered_levels)
}


ccc_reverse_numeric_container <- function(x) {
  finite_vals <- unlist(x, use.names = FALSE)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  if (length(finite_vals) == 0L) {
    return(x)
  }
  rng <- range(finite_vals, na.rm = TRUE)
  rev_fun <- function(v) {
    v <- as.numeric(v)
    keep <- is.finite(v)
    out <- v
    out[keep] <- rng[2] - v[keep] + rng[1]
    out
  }
  if (is.list(x)) {
    return(lapply(x, rev_fun))
  }
  rev_fun(x)
}


ccc_axis_param <- function(which, reverse = FALSE) {
  axis_param <- list(gp = grid::gpar())
  if (identical(which, "row")) {
    axis_param$direction <- if (isTRUE(reverse)) "reverse" else "normal"
  }
  axis_param
}


ccc_build_annotation_wrapper <- function(
  anno,
  which = c("column", "row"),
  font.size = 10,
  border = TRUE,
  name = "Value"
) {
  which <- match.arg(which)
  if (is.null(anno)) {
    return(NULL)
  }
  if (identical(which, "column")) {
    return(do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(
        stats::setNames(list(anno), name),
        list(
          annotation_name_gp = grid::gpar(fontsize = font.size * 0.8),
          annotation_name_side = "left",
          show_annotation_name = TRUE,
          which = "column",
          border = border,
          gap = grid::unit(1, "mm")
        )
      )
    ))
  }
  do.call(
    ComplexHeatmap::rowAnnotation,
    c(
      stats::setNames(list(anno), name),
      list(
        annotation_name_gp = grid::gpar(fontsize = font.size * 0.8),
        annotation_name_side = "top",
        show_annotation_name = TRUE,
        border = border,
        gap = grid::unit(1, "mm")
      )
    )
  )
}


ccc_build_bar_annotation <- function(
  stats_list,
  fill_cols,
  which = c("column", "row"),
  font.size = 10,
  border = TRUE,
  reverse = FALSE
) {
  which <- match.arg(which)
  if (is.null(stats_list) || is.null(fill_cols)) {
    return(NULL)
  }
  wrapper_args <- c(
    ccc_build_bar_annotation_args(
      stats_list = stats_list,
      fill_cols = fill_cols,
      which = which,
      font.size = font.size,
      border = border,
      reverse = reverse
    ),
    list(
      annotation_name_gp = grid::gpar(fontsize = font.size * 0.8),
      annotation_name_side = if (identical(which, "column")) "left" else "top",
      show_annotation_name = TRUE,
      border = border,
      gap = grid::unit(1, "mm")
    )
  )
  if (identical(which, "column")) {
    wrapper_args$which <- "column"
    return(do.call(ComplexHeatmap::HeatmapAnnotation, wrapper_args))
  }
  do.call(ComplexHeatmap::rowAnnotation, wrapper_args)
}


ccc_build_distribution_annotation <- function(
  values,
  fill_cols,
  which = c("column", "row"),
  font.size = 10,
  border = TRUE,
  type = c("box", "histogram", "density", "violin"),
  reverse = FALSE
) {
  which <- match.arg(which)
  type <- match.arg(type)
  if (is.null(values) || is.null(fill_cols)) {
    return(NULL)
  }
  values_use <- values
  if (isTRUE(reverse) && type %in% c("histogram", "density", "violin")) {
    values_use <- ccc_reverse_numeric_container(values_use)
  }
  if (
    type %in%
      c("histogram", "density", "violin") &&
      !ccc_distribution_has_enough_points(values_use)
  ) {
    type <- "box"
  }
  anno <- switch(type,
    box = ComplexHeatmap::anno_boxplot(
      values_use,
      which = which,
      gp = grid::gpar(
        fill = fill_cols,
        col = if (isTRUE(border)) "grey25" else NA
      ),
      border = border,
      outline = FALSE,
      axis = FALSE,
      axis_param = ccc_axis_param(which = which, reverse = reverse),
      box_width = 0.7,
      width = if (identical(which, "row")) grid::unit(2, "cm") else NULL,
      height = if (identical(which, "column")) grid::unit(2, "cm") else NULL
    ),
    histogram = ComplexHeatmap::anno_histogram(
      values_use,
      which = which,
      gp = grid::gpar(
        fill = fill_cols,
        col = if (isTRUE(border)) "grey25" else NA
      ),
      border = border,
      axis = FALSE,
      width = if (identical(which, "row")) grid::unit(2, "cm") else NULL,
      height = if (identical(which, "column")) grid::unit(2, "cm") else NULL
    ),
    density = ComplexHeatmap::anno_density(
      values_use,
      which = which,
      type = "lines",
      gp = grid::gpar(col = fill_cols, fill = fill_cols, lwd = 1),
      border = border,
      axis = FALSE,
      width = if (identical(which, "row")) grid::unit(2, "cm") else NULL,
      height = if (identical(which, "column")) grid::unit(2, "cm") else NULL
    ),
    violin = ComplexHeatmap::anno_density(
      values_use,
      which = which,
      type = "violin",
      gp = grid::gpar(
        fill = fill_cols,
        col = if (isTRUE(border)) "grey25" else NA
      ),
      border = border,
      axis = FALSE,
      width = if (identical(which, "row")) grid::unit(2, "cm") else NULL,
      height = if (identical(which, "column")) grid::unit(2, "cm") else NULL
    )
  )
  ccc_build_annotation_wrapper(
    anno = anno,
    which = which,
    font.size = font.size,
    border = border,
    name = tools::toTitleCase(type)
  )
}


ccc_distribution_has_enough_points <- function(values) {
  values_list <- if (is.list(values)) {
    values
  } else if (is.matrix(values)) {
    lapply(seq_len(ncol(values)), function(i) values[, i])
  } else {
    list(values)
  }
  values_list <- lapply(values_list, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[is.finite(x)]
  })
  values_list <- Filter(length, values_list)
  if (length(values_list) == 0L) {
    return(FALSE)
  }
  all(vapply(
    values_list,
    function(x) {
      length(x) >= 2L && length(unique(x)) >= 2L
    },
    logical(1)
  ))
}


ccc_build_summary_annotation <- function(
  values,
  fill_cols,
  which = c("column", "row"),
  font.size = 10,
  border = TRUE,
  type = c("point", "line"),
  reverse = FALSE
) {
  which <- match.arg(which)
  type <- match.arg(type)
  if (is.null(values) || is.null(fill_cols)) {
    return(NULL)
  }
  anno <- switch(type,
    point = ComplexHeatmap::anno_points(
      values,
      which = which,
      gp = grid::gpar(col = fill_cols),
      pch = 16,
      size = grid::unit(1.8, "mm"),
      border = border,
      axis = FALSE,
      axis_param = ccc_axis_param(which = which, reverse = reverse),
      width = if (identical(which, "row")) grid::unit(2, "cm") else NULL,
      height = if (identical(which, "column")) grid::unit(2, "cm") else NULL
    ),
    line = ComplexHeatmap::anno_lines(
      values,
      which = which,
      gp = grid::gpar(col = "grey25", lwd = 1),
      add_points = TRUE,
      pt_gp = grid::gpar(col = fill_cols),
      pch = 16,
      size = grid::unit(1.6, "mm"),
      border = border,
      axis = FALSE,
      axis_param = ccc_axis_param(which = which, reverse = reverse),
      width = if (identical(which, "row")) grid::unit(2, "cm") else NULL,
      height = if (identical(which, "column")) grid::unit(2, "cm") else NULL
    )
  )
  ccc_build_annotation_wrapper(
    anno = anno,
    which = which,
    font.size = font.size,
    border = border,
    name = tools::toTitleCase(type)
  )
}


ccc_build_bar_annotation_args <- function(
  stats_list,
  fill_cols,
  which = c("column", "row"),
  font.size = 10,
  border = TRUE,
  reverse = FALSE
) {
  which <- match.arg(which)
  out <- list()
  for (metric in names(stats_list)) {
    values <- as.numeric(stats_list[[metric]])
    values[!is.finite(values)] <- 0
    anno_args <- list(
      values,
      gp = grid::gpar(
        fill = fill_cols,
        col = if (isTRUE(border)) "grey25" else NA,
        lwd = 0.6
      ),
      which = which,
      border = border,
      bar_width = 0.8
    )
    anno_args$axis_param <- list(
      gp = grid::gpar(fontsize = font.size * 0.7)
    )
    if (identical(which, "row")) {
      anno_args$axis_param$direction <- if (isTRUE(reverse)) {
        "reverse"
      } else {
        "normal"
      }
    }
    if (identical(which, "column")) {
      anno_args$height <- grid::unit(2, "cm")
    } else {
      anno_args$width <- grid::unit(2, "cm")
    }
    out[[tools::toTitleCase(metric)]] <- do.call(
      ComplexHeatmap::anno_barplot,
      anno_args
    )
  }
  out
}


ccc_heatmap_body_size <- function(
  nrow,
  ncol,
  width = NULL,
  height = NULL,
  units = "inch",
  default_cell_size = 0.35
) {
  nrow <- max(as.numeric(nrow), 1)
  ncol <- max(as.numeric(ncol), 1)
  width_num <- if (is.numeric(width) && length(width) > 0L) {
    as.numeric(width[1])
  } else {
    NA_real_
  }
  height_num <- if (is.numeric(height) && length(height) > 0L) {
    as.numeric(height[1])
  } else {
    NA_real_
  }
  fixed <- !all(is.na(c(width_num, height_num)))
  if (is.na(width_num) && is.na(height_num)) {
    width_num <- ncol * default_cell_size
    height_num <- nrow * default_cell_size
  } else if (is.na(width_num)) {
    width_num <- height_num * ncol / nrow
  } else if (is.na(height_num)) {
    height_num <- width_num * nrow / ncol
  }
  list(
    width = grid::unit(width_num, units),
    height = grid::unit(height_num, units),
    width_num = width_num,
    height_num = height_num,
    fixed = fixed
  )
}


ccc_heatmap_capture_size <- function(
  body_width,
  body_height,
  top_tracks = 0,
  bottom_tracks = 0,
  left_tracks = 0,
  right_tracks = 0,
  legend.position = "right",
  has_title = FALSE
) {
  extra_width <- 0.8 + left_tracks * 0.85 + right_tracks * 0.85
  extra_height <- 0.8 + top_tracks * 0.45 + bottom_tracks * 0.35
  if (legend.position %in% c("right", "left")) {
    extra_width <- extra_width + 1
  }
  if (legend.position %in% c("top", "bottom")) {
    extra_height <- extra_height + 0.8
  }
  if (isTRUE(has_title)) {
    extra_height <- extra_height + 0.4
  }
  list(
    width = body_width + extra_width,
    height = body_height + extra_height
  )
}
