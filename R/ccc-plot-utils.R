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
