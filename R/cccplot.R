get_cc_obj <- function(x) {
  if (!inherits(x, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  x@tools[["CellChat"]]
}

cc_names <- function(srt) {
  store <- get_cc_obj(srt)
  if (is.null(store) || is.null(store$results)) {
    return(character(0))
  }
  names(store$results)
}

comparison_cc_names <- function(srt) {
  store <- get_cc_obj(srt)
  if (is.null(store) || is.null(store$comparisons)) {
    return(character(0))
  }
  names(store$comparisons)
}

resolve_single_cc_condition <- function(srt, condition = NULL) {
  result_names <- cc_names(srt)
  if (!is.null(condition)) {
    if (!condition %in% result_names) {
      log_message(
        "Condition {.val {condition}} not found in CellChat results",
        message_type = "error"
      )
    }
    return(condition)
  }
  if (length(result_names) == 1L) {
    return(result_names[1])
  }
  if ("ALL" %in% result_names) {
    return("ALL")
  }
  log_message(
    "Multiple CellChat results found. Please specify {.arg condition}",
    message_type = "error"
  )
}

use_cc_single_condition <- function(srt, condition = NULL) {
  result_names <- cc_names(srt)
  cmp_names <- comparison_cc_names(srt)
  if (!is.null(condition)) {
    if (condition %in% result_names) {
      return(TRUE)
    }
    if (condition %in% cmp_names) {
      return(FALSE)
    }
    log_message(
      "{.arg condition} must be one of CellChat result names or comparison names",
      message_type = "error"
    )
  }
  if (length(result_names) == 1L && length(cmp_names) == 0L) {
    return(TRUE)
  }
  if ("ALL" %in% result_names && length(cmp_names) == 0L) {
    return(TRUE)
  }
  if (length(cmp_names) == 1L && length(result_names) == 0L) {
    return(FALSE)
  }
  log_message(
    "The requested plot is ambiguous because both single-condition results and/or comparison results are available. Please specify {.arg condition}",
    message_type = "error"
  )
}

resolve_group_index_single <- function(object, group.use = NULL) {
  if (is.null(group.use) || length(group.use) == 0L) {
    return(NULL)
  }

  groups <- character(0)
  if (inherits(object, "CellChat")) {
    groups <- tryCatch(
      {
        if (!is.null(object@idents)) {
          unique(as.character(object@idents))
        } else {
          character(0)
        }
      },
      error = function(e) character(0)
    )
  }
  if (length(groups) == 0L) {
    groups <- unique(as.character(group.use))
  }

  group.use_chr <- as.character(group.use)
  if (length(group.use_chr) == 1L && identical(group.use_chr, "all")) {
    return(groups)
  }

  group.use_num <- suppressWarnings(as.integer(group.use_chr))
  if (
    length(group.use_chr) > 0L &&
      all(!is.na(group.use_num)) &&
      all(group.use_chr == as.character(group.use_num))
  ) {
    if (any(group.use_num < 1L) || any(group.use_num > length(groups))) {
      log_message(
        "{.arg group.use} indices are out of range for available groups",
        message_type = "error"
      )
    }
    return(groups[group.use_num])
  }

  missing <- setdiff(group.use_chr, groups)
  if (length(missing) > 0L) {
    log_message(
      "Unknown group labels in {.arg group.use}: {.val {missing}}",
      message_type = "error"
    )
  }

  groups[groups %in% group.use_chr]
}

get_single_cc_obj <- function(srt, condition = NULL) {
  store <- get_cc_obj(srt)
  condition <- resolve_single_cc_condition(srt, condition = condition)
  store$results[[condition]]$cellchat_object
}

.cc_get_cmp <- function(srt, condition = NULL) {
  store <- get_cc_obj(srt)
  cmp_names <- comparison_cc_names(srt)
  if (is.null(condition)) {
    if (length(cmp_names) == 1L) {
      condition <- cmp_names[1]
    } else if (length(cmp_names) == 0L) {
      log_message(
        "No CellChat comparisons found",
        message_type = "error"
      )
    } else {
      log_message(
        "Multiple CellChat comparisons found. Please specify {.arg condition}",
        message_type = "error"
      )
    }
  }
  if (!condition %in% cmp_names) {
    log_message(
      "Comparison {.val {condition}} not found in CellChat results",
      message_type = "error"
    )
  }
  store$comparisons[[condition]]
}

get_dataset_object <- function(srt, condition = NULL, dataset = 1) {
  store <- get_cc_obj(srt)
  result_names <- cc_names(srt)
  cmp_names <- comparison_cc_names(srt)

  if (!is.null(condition) && condition %in% result_names) {
    return(list(
      object = store$results[[condition]]$cellchat_object,
      seurat_object = store$results[[condition]]$seurat_object,
      label = condition,
      source = "single"
    ))
  }

  if (!is.null(condition) && condition %in% cmp_names) {
    cmp <- .cc_get_cmp(srt, condition = condition)
    ds_name <- .cc_pick_dataset_name(cmp, dataset)
    seu_object <- NULL
    if (
      !is.null(store$results[[ds_name]]) &&
        !is.null(store$results[[ds_name]]$seurat_object)
    ) {
      seu_object <- store$results[[ds_name]]$seurat_object
    }
    return(list(
      object = cmp$object.list[[ds_name]],
      seurat_object = seu_object,
      label = ds_name,
      source = "comparison"
    ))
  }

  if (is.null(condition) && length(result_names) == 1L) {
    return(list(
      object = store$results[[result_names[1]]]$cellchat_object,
      seurat_object = store$results[[result_names[1]]]$seurat_object,
      label = result_names[1],
      source = "single"
    ))
  }

  if (is.null(condition) && "ALL" %in% result_names) {
    return(list(
      object = store$results[["ALL"]]$cellchat_object,
      seurat_object = store$results[["ALL"]]$seurat_object,
      label = "ALL",
      source = "single"
    ))
  }

  if (is.null(condition) && length(cmp_names) == 1L) {
    cmp <- .cc_get_cmp(srt, condition = cmp_names[1])
    ds_name <- .cc_pick_dataset_name(cmp, dataset)
    seu_object <- NULL
    if (
      !is.null(store$results[[ds_name]]) &&
        !is.null(store$results[[ds_name]]$seurat_object)
    ) {
      seu_object <- store$results[[ds_name]]$seurat_object
    }
    return(list(
      object = cmp$object.list[[ds_name]],
      seurat_object = seu_object,
      label = ds_name,
      source = "comparison"
    ))
  }

  log_message(
    "Unable to determine which CellChat object to plot. Please specify {.arg condition}",
    message_type = "error"
  )
}

.cc_resolve_dataset_index <- function(cmp, comparison = c(1, 2)) {
  ds_names <- names(cmp$object.list)
  if (is.numeric(comparison)) {
    idx <- as.integer(comparison)
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(ds_names))) {
      log_message(
        "comparison indices out of range. Available datasets: {.val {seq_along(ds_names)}}",
        message_type = "error"
      )
    }
    return(idx)
  }
  idx <- match(as.character(comparison), ds_names)
  if (any(is.na(idx))) {
    log_message(
      "comparison names not found: {.val {as.character(comparison)[is.na(idx)]}}. Available datasets: {.val {ds_names}}",
      message_type = "error"
    )
  }
  idx
}

.cc_pick_dataset_name <- function(cmp, dataset = 1) {
  nm <- names(cmp$object.list)
  if (is.character(dataset)) {
    if (!dataset %in% nm) {
      log_message(
        "dataset {.val {dataset}} not found. Available: {.val {nm}}",
        message_type = "error"
      )
    }
    return(dataset)
  }
  if (length(dataset) != 1L || dataset < 1L || dataset > length(nm)) {
    log_message(
      "dataset index out of range. Available: {.val {seq_along(nm)}}",
      message_type = "error"
    )
  }
  nm[dataset]
}

cmp_cc_label <- function(cmp, comp_idx = c(1, 2)) {
  nm <- names(cmp$object.list)[comp_idx]
  if (length(nm) <= 1L) {
    return(nm[1])
  }
  paste0(nm[1], "_vs_", nm[2])
}

cc_theme <- function(theme_use = "theme_scop", theme_args = list()) {
  if (is.null(theme_use)) {
    return(NULL)
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  theme_fun <- if (is.character(theme_use)) {
    tryCatch(
      get(theme_use, mode = "function", inherits = TRUE),
      error = function(e) NULL
    )
  } else if (is.function(theme_use)) {
    theme_use
  } else {
    NULL
  }

  if (is.null(theme_fun)) {
    theme_fun <- ggplot2::theme_bw
  }

  do.call(theme_fun, theme_args)
}

finalize_cc_plot <- function(
  p,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  aspect.ratio = NULL,
  font.size = NULL
) {
  if (!inherits(p, c("gg", "ggplot"))) {
    return(p)
  }
  theme_obj <- cc_theme(theme_use = theme_use, theme_args = theme_args)

  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = legend.title,
      fill = legend.title,
      size = legend.title
    ) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  if (!is.null(theme_obj)) {
    p <- p + theme_obj
  }
  if (!is.null(aspect.ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect.ratio)
  }
  if (!is.null(font.size)) {
    p <- p +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = font.size),
        axis.title = ggplot2::element_text(size = font.size + 1),
        plot.title = ggplot2::element_text(size = font.size + 2, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = font.size),
        legend.text = ggplot2::element_text(size = font.size),
        legend.title = ggplot2::element_text(size = font.size + 1)
      )
  }
  p
}

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
  theme_obj <- cc_theme(theme_use = theme_use, theme_args = theme_args)
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

subset_cc_table <- function(
  object,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05,
  dataset = NULL
) {
  if (!is.null(pairLR.use) && !is.data.frame(pairLR.use)) {
    pairLR.use <- data.frame(
      interaction_name = as.character(pairLR.use),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(pairLR.use)) {
    signaling <- NULL
  }
  df <- CellChat::subsetCommunication(
    object = object,
    slot.name = slot.name,
    sources.use = sources.use,
    targets.use = targets.use,
    signaling = signaling,
    pairLR.use = pairLR.use,
    thresh = thresh
  )
  df <- as.data.frame(df)
  if (!is.null(dataset) && nrow(df) > 0L) {
    df$dataset <- dataset
  }
  df
}

prepare_cc_bubble_data <- function(
  object,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05,
  dataset = NULL,
  top_n = 10,
  remove.isolate = TRUE
) {
  df <- subset_cc_table(
    object = object,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sources.use,
    targets.use = targets.use,
    thresh = thresh,
    dataset = dataset
  )
  if (is.null(df) || nrow(df) == 0L) {
    return(df)
  }

  if (isTRUE(remove.isolate)) {
    df <- df[df$prob > 0, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    return(df)
  }

  if (
    !is.null(top_n) &&
      is.numeric(top_n) &&
      length(top_n) == 1L &&
      top_n > 0L &&
      is.null(pairLR.use)
  ) {
    if (length(unique(df$interaction_name)) > top_n) {
      score_interaction <- sort(
        tapply(df$prob, as.character(df$interaction_name), sum, na.rm = TRUE),
        decreasing = TRUE
      )
      keep_interaction <- names(score_interaction)[seq_len(min(
        top_n,
        length(score_interaction)
      ))]
      df <- df[df$interaction_name %in% keep_interaction, , drop = FALSE]
    }

    max_pairs <- max(15L, as.integer(top_n) * 2L)
    pair_vec <- paste(df$source, "->", df$target)
    if (
      length(unique(pair_vec)) > max_pairs &&
        is.null(sources.use) &&
        is.null(targets.use)
    ) {
      score_pair <- sort(
        tapply(df$prob, pair_vec, sum, na.rm = TRUE),
        decreasing = TRUE
      )
      keep_pair <- names(score_pair)[seq_len(min(
        max_pairs,
        length(score_pair)
      ))]
      df <- df[paste(df$source, "->", df$target) %in% keep_pair, , drop = FALSE]
    }
  }

  df
}

custom_cc_bubble_plot <- function(
  df,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  color.by = c("prob", "pval"),
  bubble_size.range = c(1.5, 8),
  palette = "RdBu",
  palcolor = NULL,
  font.size = 10,
  angle.x = 45,
  hjust.x = 1,
  vjust.x = 1,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  aspect.ratio = NULL,
  remove.isolate = TRUE
) {
  color.by <- match.arg(color.by)
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No communication records available for bubble plot",
      message_type = "error"
    )
  }

  req_cols <- c("source", "target", "interaction_name", "prob", "pval")
  miss <- setdiff(req_cols, colnames(df))
  if (length(miss) > 0L) {
    log_message(
      "Missing required columns for bubble plot: {.val {miss}}",
      message_type = "error"
    )
  }

  if (isTRUE(remove.isolate)) {
    df <- df[df$prob > 0, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No non-zero communication records remain after filtering",
      message_type = "error"
    )
  }

  df$pair <- paste(df$source, "->", df$target)
  if ("dataset" %in% colnames(df) && length(unique(df$dataset)) > 1L) {
    df$pair <- paste0(df$pair, " [", df$dataset, "]")
  }
  df$interaction_plot <- as.character(df$interaction_name)

  pair_order <- sort(
    tapply(df$prob, df$pair, sum, na.rm = TRUE),
    decreasing = TRUE
  )
  interaction_order <- sort(
    tapply(df$prob, df$interaction_plot, sum, na.rm = TRUE),
    decreasing = TRUE
  )
  df$pair <- factor(df$pair, levels = names(pair_order))
  df$interaction_plot <- factor(
    df$interaction_plot,
    levels = rev(names(interaction_order))
  )

  cols <- palette_colors(palette = palette, palcolor = palcolor, n = 9)
  fill_var <- if (identical(color.by, "pval")) {
    df$fill_val <- -log10(df$pval + 1e-300)
    "fill_val"
  } else {
    df$fill_val <- df$prob
    "fill_val"
  }
  fill_label <- legend.title %||%
    if (identical(color.by, "pval")) "-log10(pval)" else "prob"

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = pair,
      y = interaction_plot,
      size = prob,
      fill = .data[[fill_var]]
    )
  ) +
    ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.2) +
    ggplot2::scale_size_area(
      name = "prob",
      max_size = max(bubble_size.range),
      n.breaks = 4,
      guide = ggplot2::guide_legend(
        override.aes = list(fill = "grey30", shape = 21, color = "grey20"),
        order = 2
      )
    ) +
    ggplot2::scale_fill_gradientn(
      name = fill_label,
      colours = cols,
      n.breaks = 4,
      guide = ggplot2::guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        title.hjust = 0,
        order = 1
      )
    ) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = angle.x,
        hjust = hjust.x,
        vjust = vjust.x
      ),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.3,
        colour = "grey80",
        linetype = 2
      ),
      panel.grid.minor = ggplot2::element_blank()
    )

  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.title = legend.title,
    theme_use = theme_use,
    theme_args = theme_args,
    aspect.ratio = aspect.ratio,
    font.size = font.size
  )
}

detect_method <- function(srt, method = NULL) {
  alias_map <- c(
    "CellPhoneDB" = "CellphoneDB",
    "CellphoneDB" = "CellphoneDB",
    "Liana" = "LIANA",
    "liana" = "LIANA",
    "NicheNet" = "Nichenetr",
    "MultiNicheNet" = "MultiNichenetr"
  )
  alias_map_lower <- c(
    "cellchat" = "CellChat",
    "cellphonedb" = "CellphoneDB",
    "cellphone_db" = "CellphoneDB",
    "cellphone db" = "CellphoneDB",
    "liana" = "LIANA",
    "nichenet" = "Nichenetr",
    "nichenetr" = "Nichenetr",
    "multinichenet" = "MultiNichenetr",
    "multinichenetr" = "MultiNichenetr"
  )
  if (!is.null(method)) {
    method_chr <- as.character(method)[1]
    if (is.na(method_chr) || !nzchar(method_chr)) {
      log_message(
        "{.arg method} must be a non-empty string or NULL to auto-detect",
        message_type = "error"
      )
    }
    method_chr <- trimws(method_chr)
    if (method_chr %in% names(alias_map)) {
      return(unname(alias_map[[method_chr]]))
    }
    key <- tolower(method_chr)
    if (key %in% names(alias_map_lower)) {
      return(unname(alias_map_lower[[key]]))
    }
    return(method_chr)
  }
  candidates <- c("CellphoneDB", "Nichenetr", "MultiNichenetr", "LIANA", "CellChat")
  available <- candidates[candidates %in% names(srt@tools)]
  if (length(available) == 1L) {
    return(available[1])
  }
  if (length(available) == 0L) {
    log_message(
      "No cell-cell communication results were found in {.cls Seurat}",
      message_type = "error"
    )
  }
  log_message(
    "Multiple CCC methods are available. Please specify {.arg method}. Candidates: {.val {available}}",
    message_type = "error"
  )
}

get_bundle <- function(srt, method) {
  x <- srt@tools[[method]]
  if (is.null(x)) {
    log_message(
      "{.pkg {method}} results not found in {.cls Seurat}",
      message_type = "error"
    )
  }
  x
}

get_group_by <- function(srt, method) {
  method <- detect_method(srt = srt, method = method)
  if (identical(method, "CellChat")) {
    store <- get_cc_obj(srt)
    return(store$parameters$group.by %||% NULL)
  }
  bundle <- get_bundle(srt, method = method)
  bundle$parameters$group.by %||% NULL
}

ccc_pair_table <- function(
  srt,
  method,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05
) {
  method <- detect_method(srt = srt, method = method)
  if (identical(method, "CellChat")) {
    df <- extract_long_table(
      srt = srt,
      condition = condition,
      dataset = dataset,
      slot.name = slot.name,
      signaling = signaling,
      pairLR.use = pairLR.use,
      sources.use = sources.use,
      targets.use = targets.use,
      thresh = thresh
    )
  } else {
    bundle <- get_bundle(srt, method = method)
    df <- bundle$long_table %||% data.frame()
  }
  df <- standardize_long_df(df)
  df <- ccc_mark_significance(df, thresh = thresh)
  aggregate_ccc_long(df)
}

extract_long_table <- function(
  srt,
  condition = NULL,
  dataset = 1,
  slot.name = "net",
  signaling = NULL,
  pairLR.use = NULL,
  sources.use = NULL,
  targets.use = NULL,
  thresh = 0.05
) {
  info <- get_dataset_object(srt, condition = condition, dataset = dataset)
  df <- subset_cc_table(
    object = info$object,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sources.use,
    targets.use = targets.use,
    thresh = thresh,
    dataset = info$label
  )
  df <- standardize_long_df(df)
  df <- ccc_mark_significance(df, thresh = thresh)
  df$method <- "CellChat"
  df
}

standardize_long_df <- function(df) {
  df <- standardize_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }

  rename_map <- list(
    sender = c("sender", "source"),
    receiver = c("receiver", "target"),
    interaction_name = c("interaction_name", "interacting_pair", "interaction"),
    ligand = c("ligand", "ligand_complex", "gene_a", "partner_a", "from"),
    receptor = c("receptor", "receptor_complex", "gene_b", "partner_b", "to"),
    pathway_name = c("pathway_name", "classification", "signaling"),
    score = c("score", "prob", "means", "mean", "prioritization_score", "LRscore", "lrscore", "magnitude"),
    pvalue = c("pvalue", "pval", "pvalues", "cellphone_pvals", "aggregate_rank", "specificity_rank", "magnitude_rank")
  )
  out <- df
  for (nm in names(rename_map)) {
    col <- ccc_pick_col(out, rename_map[[nm]])
    if (!is.null(col) && !nm %in% colnames(out)) {
      colnames(out)[match(col, colnames(out))] <- nm
    }
  }
  for (nm in c(
    "sender",
    "receiver",
    "interaction_name",
    "ligand",
    "receptor",
    "pathway_name"
  )) {
    if (!nm %in% colnames(out)) {
      out[[nm]] <- NA_character_
    }
  }
  interaction_label_col <- ccc_pick_col(
    out,
    c("interaction_label", "interaction_name_2", "interaction_name", "interacting_pair")
  )
  if (!is.null(interaction_label_col) && !"interaction_label" %in% colnames(out)) {
    out[["interaction_label"]] <- out[[interaction_label_col]]
  }
  if (!"interaction_label" %in% colnames(out)) {
    out[["interaction_label"]] <- out[["interaction_name"]]
  }
  out[["interaction_label"]] <- ccc_display_interaction(out[["interaction_label"]])
  if (!"classification" %in% colnames(out)) {
    out[["classification"]] <- out[["pathway_name"]]
  }
  out[["classification"]] <- as.character(out[["classification"]])
  out[["classification"]][is.na(out[["classification"]]) | !nzchar(out[["classification"]])] <- "Unclassified"
  out[["pathway_name"]] <- out[["classification"]]
  pair_lr_col <- ccc_pick_col(
    out,
    c("pair_lr", "pairLR", "interacting_pair", "interaction_name_2", "interaction_name")
  )
  if (!is.null(pair_lr_col) && !"pair_lr" %in% colnames(out)) {
    out[["pair_lr"]] <- out[[pair_lr_col]]
  }
  if (!"pair_lr" %in% colnames(out)) {
    out[["pair_lr"]] <- paste(out[["ligand"]], out[["receptor"]], sep = "-")
  }
  if (!"score" %in% colnames(out)) {
    out$score <- NA_real_
  }
  if (!"pvalue" %in% colnames(out)) {
    out$pvalue <- NA_real_
  }
  if (!"ligand_display" %in% colnames(out)) {
    out$ligand_display <- ccc_display_gene(out$ligand)
  }
  if (!"receptor_display" %in% colnames(out)) {
    out$receptor_display <- ccc_display_gene(out$receptor)
  }
  if (!"interaction_display" %in% colnames(out)) {
    out$interaction_display <- out$interaction_label
  }
  out <- ccc_mark_significance(out)
  out
}

ccc_display_gene <- function(x) {
  x <- ccc_clean_identifier(x, drop_receptor_suffix = TRUE)
  x[is.na(x)] <- ""
  x
}

ccc_display_interaction <- function(x) {
  x <- ccc_clean_identifier(x, drop_receptor_suffix = FALSE)
  x[is.na(x)] <- ""
  x <- gsub("_", " - ", x, fixed = TRUE)
  x <- gsub("\\s*-\\s*", " - ", x)
  x
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

ccc_clean_identifier <- function(x, drop_receptor_suffix = FALSE) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- sub("^(complex:|simple:)", "", x, ignore.case = TRUE)
  x <- sub("(?:_|\\s)complex$", "", x, ignore.case = TRUE)
  if (isTRUE(drop_receptor_suffix)) {
    x <- sub("^integrin[_\\s]+", "", x, ignore.case = TRUE)
    x <- sub("_receptor_inhibitor$", "", x, ignore.case = TRUE)
    x <- sub("_receptor$", "", x, ignore.case = TRUE)
    x <- sub("_ligand$", "", x, ignore.case = TRUE)
  }
  x
}

ccc_identifier_key <- function(x) {
  x <- ccc_clean_identifier(x, drop_receptor_suffix = FALSE)
  x <- tolower(trimws(x))
  x[x == "nan"] <- ""
  gsub("[^[:alnum:]]+", "", x)
}

ccc_row_identifier_candidates <- function(df, columns = character(), include_lr = TRUE) {
  if (is.null(df) || nrow(df) == 0L) {
    return(list())
  }
  out <- list()
  for (nm in intersect(columns, colnames(df))) {
    out[[nm]] <- as.character(df[[nm]])
  }
  if (isTRUE(include_lr) && all(c("ligand", "receptor") %in% colnames(df))) {
    ligand <- ccc_display_gene(df$ligand)
    receptor <- ccc_display_gene(df$receptor)
    out[["lr_dash"]] <- paste(ligand, receptor, sep = "-")
    out[["lr_spaced"]] <- paste(ligand, receptor, sep = " - ")
    out[["lr_arrow"]] <- paste(ligand, receptor, sep = " -> ")
    out[["lr_underscore"]] <- paste(ligand, receptor, sep = "_")
  }
  out
}

ccc_match_identifiers <- function(
  df,
  values,
  columns,
  include_lr = TRUE
) {
  values <- ccc_normalize_use_arg(values)
  if (is.null(values)) {
    return(rep(TRUE, nrow(df)))
  }
  candidates <- ccc_row_identifier_candidates(
    df = df,
    columns = columns,
    include_lr = include_lr
  )
  if (length(candidates) == 0L) {
    return(rep(FALSE, nrow(df)))
  }
  values_key <- ccc_identifier_key(values)
  keep <- rep(FALSE, nrow(df))
  for (candidate in candidates) {
    candidate_chr <- as.character(candidate)
    candidate_key <- ccc_identifier_key(candidate_chr)
    keep <- keep |
      candidate_chr %in% values |
      candidate_key %in% values_key
  }
  keep[is.na(keep)] <- FALSE
  keep
}

ccc_filter_bubble_df <- function(
  df,
  interaction.use = NULL,
  pairLR.use = NULL
) {
  if (is.null(df) || nrow(df) == 0L) {
    return(df)
  }
  if (!is.null(interaction.use)) {
    keep <- ccc_match_identifiers(
      df = df,
      values = interaction.use,
      columns = c(
        "interaction_name",
        "interaction_name_2",
        "interaction_label",
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
        "interaction_name",
        "interaction_name_2"
      ),
      include_lr = TRUE
    )
    df <- df[keep, , drop = FALSE]
  }
  df
}

ccc_panel_grid <- function(n_panels, ncols = NULL, nrows = NULL, context = "plot") {
  n_panels <- max(as.integer(n_panels %||% 0L), 1L)
  has_ncols <- !is.null(ncols)
  has_nrows <- !is.null(nrows)
  if (!is.null(ncols)) {
    if (!is.numeric(ncols) || length(ncols) != 1L || is.na(ncols) || ncols <= 0) {
      log_message("{.arg ncols} must be a positive integer", message_type = "error")
    }
    ncols <- as.integer(ncols)
  }
  if (!is.null(nrows)) {
    if (!is.numeric(nrows) || length(nrows) != 1L || is.na(nrows) || nrows <= 0) {
      log_message("{.arg nrows} must be a positive integer", message_type = "error")
    }
    nrows <- as.integer(nrows)
  }
  if (!is.null(ncols) && !is.null(nrows) && ncols * nrows < n_panels) {
    log_message(
      "{.arg ncols * nrows} must be at least the number of panels for {context}",
      message_type = "error"
    )
  }
  if (is.null(ncols) && is.null(nrows)) {
    ncols <- n_panels
    nrows <- 1L
  } else if (is.null(ncols)) {
    ncols <- ceiling(n_panels / nrows)
  } else if (is.null(nrows)) {
    nrows <- ceiling(n_panels / ncols)
  }
  if (!isTRUE(has_ncols)) {
    ncols <- min(as.integer(ncols), n_panels)
  }
  if (!isTRUE(has_nrows)) {
    nrows <- ceiling(n_panels / max(as.integer(ncols), 1L))
  }
  ncols <- as.integer(ncols)
  nrows <- as.integer(nrows)
  list(ncols = ncols, nrows = nrows)
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
        class_key %in% signaling_key,
      ,
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
  pair_df <- aggregate_ccc_long(df_use)
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

scale_var <- function(df, color.by = "score", agg_value = "sum") {
  if (identical(color.by, "pvalue") || identical(color.by, "specificity")) {
    var <- "specificity"
    label <- "-log10(pvalue)"
  } else if (agg_value %in% colnames(df)) {
    var <- agg_value
    label <- agg_value
  } else {
    var <- "score"
    label <- "score"
  }
  list(var = var, label = label)
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

ccc_scatter_df <- function(pair_df) {
  if (is.null(pair_df) || nrow(pair_df) == 0L) {
    return(data.frame())
  }
  value_col <- c("sum", "mean", "max", "count")
  value_col <- value_col[value_col %in% colnames(pair_df)][1]
  value_col <- value_col %||% "score"
  if (!value_col %in% colnames(pair_df)) {
    pair_df[[value_col]] <- 0
  }

  sender_df <- stats::aggregate(
    pair_df[[value_col]],
    by = list(cell = pair_df$sender),
    FUN = sum,
    na.rm = TRUE
  )
  colnames(sender_df)[2] <- "outgoing"

  receiver_df <- stats::aggregate(
    pair_df[[value_col]],
    by = list(cell = pair_df$receiver),
    FUN = sum,
    na.rm = TRUE
  )
  colnames(receiver_df)[2] <- "incoming"

  degree_df <- stats::aggregate(
    rep(1, nrow(pair_df)),
    by = list(cell = pair_df$sender),
    FUN = sum,
    na.rm = TRUE
  )
  colnames(degree_df)[2] <- "outgoing_links"

  degree_in_df <- stats::aggregate(
    rep(1, nrow(pair_df)),
    by = list(cell = pair_df$receiver),
    FUN = sum,
    na.rm = TRUE
  )
  colnames(degree_in_df)[2] <- "incoming_links"

  out <- merge(sender_df, receiver_df, by = "cell", all = TRUE)
  out <- merge(out, degree_df, by = "cell", all = TRUE)
  out <- merge(out, degree_in_df, by = "cell", all = TRUE)
  out[is.na(out)] <- 0
  out$total_strength <- out$outgoing + out$incoming
  out$total_links <- out$outgoing_links + out$incoming_links
  out
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

ccc_circle_value_col <- function(value = "weight") {
  if (identical(value, "count")) {
    return("count")
  }
  if (value %in% c("sum", "mean", "max")) {
    return(value)
  }
  "sum"
}

ccc_cellchat_net_matrix <- function(object, measure = "count") {
  mat <- object@net[[measure]] %||% NULL
  if (is.null(mat)) {
    log_message(
      "CellChat network slot {.val {measure}} is not available",
      message_type = "error"
    )
  }
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) && !is.null(colnames(mat))) {
    rownames(mat) <- colnames(mat)
  }
  if (is.null(colnames(mat)) && !is.null(rownames(mat))) {
    colnames(mat) <- rownames(mat)
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    lev <- levels(object@idents)
    if (length(lev) == nrow(mat) && length(lev) == ncol(mat)) {
      rownames(mat) <- lev
      colnames(mat) <- lev
    } else {
      rn <- rownames(mat) %||% paste0("row", seq_len(nrow(mat)))
      cn <- colnames(mat) %||% paste0("col", seq_len(ncol(mat)))
      rownames(mat) <- rn
      colnames(mat) <- cn
    }
  }
  storage.mode(mat) <- "numeric"
  mat
}

ccc_align_named_matrix <- function(mat, row_levels, col_levels = row_levels) {
  out <- matrix(
    0,
    nrow = length(row_levels),
    ncol = length(col_levels),
    dimnames = list(row_levels, col_levels)
  )
  if (is.null(mat) || length(mat) == 0L) {
    return(out)
  }
  rn <- intersect(rownames(mat), row_levels)
  cn <- intersect(colnames(mat), col_levels)
  if (length(rn) > 0L && length(cn) > 0L) {
    out[rn, cn] <- mat[rn, cn, drop = FALSE]
  }
  out
}

ccc_cellchat_diff_network_data <- function(
  srt,
  condition = NULL,
  comparison = c(1, 2),
  measure = "count",
  sender.use = NULL,
  receiver.use = NULL,
  top_n = 20,
  edge_threshold = 0
) {
  cmp <- .cc_get_cmp(srt = srt, condition = condition)
  comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)
  if (length(comp_idx) < 2L) {
    log_message(
      "{.arg comparison} must contain at least two datasets for {.val plot_type = 'diff_network'}",
      message_type = "error"
    )
  }
  comp_idx <- comp_idx[seq_len(2)]
  object_names <- names(cmp$object.list)[comp_idx]
  object_list <- cmp$object.list[object_names]
  mats <- lapply(object_list, ccc_cellchat_net_matrix, measure = measure)
  node_levels <- Reduce(
    union,
    lapply(mats, function(mat) unique(c(rownames(mat), colnames(mat))))
  )
  mats <- lapply(mats, function(mat) {
    ccc_align_named_matrix(mat, row_levels = node_levels, col_levels = node_levels)
  })
  diff_mat <- mats[[2]] - mats[[1]]
  diff_mat[!is.finite(diff_mat)] <- 0

  if (!is.null(sender.use)) {
    keep_rows <- intersect(as.character(sender.use), rownames(diff_mat))
    diff_mat <- diff_mat[keep_rows, , drop = FALSE]
  }
  if (!is.null(receiver.use)) {
    keep_cols <- intersect(as.character(receiver.use), colnames(diff_mat))
    diff_mat <- diff_mat[, keep_cols, drop = FALSE]
  }
  if (nrow(diff_mat) == 0L || ncol(diff_mat) == 0L) {
    log_message(
      "No sender-receiver groups remain after filtering for {.val plot_type = 'diff_network'}",
      message_type = "error"
    )
  }

  edge_df <- expand.grid(
    sender = rownames(diff_mat),
    receiver = colnames(diff_mat),
    stringsAsFactors = FALSE
  )
  edge_df$diff <- as.numeric(diff_mat)
  edge_df$abs_diff <- abs(edge_df$diff)
  edge_df <- edge_df[is.finite(edge_df$diff), , drop = FALSE]
  edge_df <- edge_df[edge_df$abs_diff > edge_threshold, , drop = FALSE]
  if (nrow(edge_df) == 0L) {
    log_message(
      "No differential CCC edges remain after thresholding",
      message_type = "error"
    )
  }
  if (is.numeric(top_n) && length(top_n) == 1L && is.finite(top_n) && top_n > 0L) {
    ord <- order(edge_df$abs_diff, decreasing = TRUE, na.last = TRUE)
    edge_df <- edge_df[utils::head(ord, top_n), , drop = FALSE]
  }
  rownames(edge_df) <- NULL

  node_levels <- unique(c(as.character(edge_df$sender), as.character(edge_df$receiver)))
  node_levels <- node_levels[!is.na(node_levels) & nzchar(node_levels)]
  node_df <- data.frame(
    node = node_levels,
    outgoing = rowSums(diff_mat[node_levels, , drop = FALSE], na.rm = TRUE),
    incoming = colSums(diff_mat[, node_levels, drop = FALSE], na.rm = TRUE),
    total_abs = rowSums(abs(diff_mat[node_levels, , drop = FALSE]), na.rm = TRUE) +
      colSums(abs(diff_mat[, node_levels, drop = FALSE]), na.rm = TRUE) -
      diag(abs(diff_mat[node_levels, node_levels, drop = FALSE])),
    balance = rowSums(diff_mat[node_levels, , drop = FALSE], na.rm = TRUE) +
      colSums(diff_mat[, node_levels, drop = FALSE], na.rm = TRUE) -
      diag(diff_mat[node_levels, node_levels, drop = FALSE]),
    stringsAsFactors = FALSE
  )
  list(
    edge_df = edge_df,
    node_df = node_df,
    object_names = object_names,
    cmp = cmp,
    diff_mat = diff_mat
  )
}

ccc_diff_edge_palette <- function(edge_color = NULL, negative_label, positive_label) {
  default_cols <- c("#2b6cb0", "#c53030")
  if (is.null(edge_color) || length(edge_color) == 0L) {
    return(stats::setNames(default_cols, c(negative_label, positive_label)))
  }
  edge_color <- as.character(edge_color)
  if (length(edge_color) == 1L) {
    return(stats::setNames(rep(edge_color, 2), c(negative_label, positive_label)))
  }
  if (!is.null(names(edge_color)) && any(nzchar(names(edge_color)))) {
    nm <- tolower(names(edge_color))
    neg_idx <- match(TRUE, nm %in% c("negative", "decrease", "down", tolower(negative_label)))
    pos_idx <- match(TRUE, nm %in% c("positive", "increase", "up", tolower(positive_label)))
    cols <- default_cols
    if (!is.na(neg_idx)) {
      cols[1] <- edge_color[neg_idx]
    }
    if (!is.na(pos_idx)) {
      cols[2] <- edge_color[pos_idx]
    }
    return(stats::setNames(cols, c(negative_label, positive_label)))
  }
  stats::setNames(edge_color[seq_len(2)], c(negative_label, positive_label))
}

ccc_diff_network_plot <- function(
  srt,
  condition = NULL,
  comparison = c(1, 2),
  measure = "count",
  sender.use = NULL,
  receiver.use = NULL,
  top_n = 20,
  edge_threshold = 0,
  edge_size = c(0.5, 1.8),
  edge_color = NULL,
  link_curvature = 0.2,
  link_alpha = 0.65,
  directed = FALSE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 5,
  node_alpha = 0.95,
  layout = "circle",
  title = NULL,
  subtitle = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  check_r("igraph", verbose = FALSE)
  dat <- ccc_cellchat_diff_network_data(
    srt = srt,
    condition = condition,
    comparison = comparison,
    measure = measure,
    sender.use = sender.use,
    receiver.use = receiver.use,
    top_n = top_n,
    edge_threshold = edge_threshold
  )
  edge_df <- dat$edge_df
  object_names <- dat$object_names
  negative_label <- paste0("Higher in ", object_names[1])
  positive_label <- paste0("Higher in ", object_names[2])
  node_levels <- unique(c(as.character(edge_df$sender), as.character(edge_df$receiver)))
  node_levels <- node_levels[!is.na(node_levels) & nzchar(node_levels)]
  if (length(node_levels) == 0L) {
    log_message(
      "No cell groups remain for differential circle plotting",
      message_type = "error"
    )
  }

  net_diff <- matrix(
    0,
    nrow = length(node_levels),
    ncol = length(node_levels),
    dimnames = list(node_levels, node_levels)
  )
  snd <- as.character(edge_df$sender)
  rcv <- as.character(edge_df$receiver)
  val <- as.numeric(edge_df$diff)
  snd_idx <- match(snd, node_levels)
  rcv_idx <- match(rcv, node_levels)
  keep <- !is.na(snd_idx) & !is.na(rcv_idx) & is.finite(val)
  if (any(keep)) {
    net_diff[cbind(snd_idx[keep], rcv_idx[keep])] <- val[keep]
  }
  net_abs <- abs(net_diff)
  node_weight <- rowSums(net_abs, na.rm = TRUE) + colSums(net_abs, na.rm = TRUE)
  if (!any(is.finite(node_weight)) || max(node_weight, na.rm = TRUE) <= 0) {
    node_weight <- rep(1, length(node_levels))
  }

  node_cols <- palette_colors(
    node_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  edge_cols <- ccc_diff_edge_palette(
    edge_color = edge_color,
    negative_label = negative_label,
    positive_label = positive_label
  )

  g <- igraph::graph_from_adjacency_matrix(
    net_abs,
    mode = "directed",
    weighted = TRUE,
    diag = TRUE
  )
  coords <- igraph::layout_in_circle(g)
  coords_scale <- if (nrow(coords) != 1L) scale(coords) else coords
  edge_start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  loop_angle <- ifelse(
    coords_scale[igraph::V(g), 1] > 0,
    -atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1]),
    pi - atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1])
  )

  vertex_size_max <- if (length(unique(node_weight)) == 1L) 5 else 15
  vertex_size <- node_weight / max(node_weight, na.rm = TRUE) * vertex_size_max + 5
  igraph::V(g)$size <- vertex_size
  igraph::V(g)$color <- grDevices::adjustcolor(
    node_cols[igraph::V(g)$name],
    alpha.f = node_alpha
  )
  igraph::V(g)$frame.color <- node_cols[igraph::V(g)$name]
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$label.cex <- max(0.8, font.size / 10)

  edge_weight_max <- max(igraph::E(g)$weight, na.rm = TRUE)
  if (!is.finite(edge_weight_max) || edge_weight_max <= 0) {
    edge_weight_max <- 1
  }
  edge_width_max <- max(edge_size) * 4
  igraph::E(g)$width <- 0.3 + igraph::E(g)$weight / edge_weight_max * edge_width_max
  edge_sign <- net_diff[cbind(
    igraph::V(g)$name[edge_start[, 1]],
    igraph::V(g)$name[edge_start[, 2]]
  )]
  edge_base_col <- ifelse(
    edge_sign >= 0,
    unname(edge_cols[positive_label]),
    unname(edge_cols[negative_label])
  )
  igraph::E(g)$color <- grDevices::adjustcolor(edge_base_col, alpha.f = link_alpha)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  if (sum(edge_start[, 2] == edge_start[, 1]) != 0) {
    loop_idx <- which(edge_start[, 2] == edge_start[, 1])
    igraph::E(g)$loop.angle[loop_idx] <- loop_angle[edge_start[loop_idx, 1]]
  }
  igraph::E(g)$arrow.size <- if (isTRUE(directed)) 0.35 else 0

  radian_rescale <- function(x, start = 0, direction = 1) {
    rotate <- function(y) (y + start) %% (2 * pi) * direction
    rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label_locs <- radian_rescale(
    x = seq_len(length(igraph::V(g))),
    direction = -1,
    start = 0
  )
  label_dist <- vertex_size / max(vertex_size) + 2

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(old_par), silent = TRUE), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, if (is.null(title) && is.null(subtitle)) 1.8 else 3, 0.5))
  plot(
    g,
    edge.curved = link_curvature,
    vertex.shape = "circle",
    layout = coords_scale,
    margin = 0.2,
    vertex.label.dist = label_dist,
    vertex.label.degree = label_locs,
    vertex.label.family = "Helvetica",
    edge.label.family = "Helvetica"
  )
  graphics::title(
    main = title %||% paste0(object_names[2], " vs ", object_names[1]),
    sub = subtitle
  )
  grDevices::recordPlot()
}

ccc_scatter_plot <- function(
  srt,
  method,
  pair_df,
  condition = NULL,
  dataset = 1,
  signaling = NULL,
  title = NULL,
  subtitle = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  link_palette = "Dark2",
  link_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  dots <- list(...)
  dot_alpha <- dots[["dot.alpha"]] %||% 0.6
  dot_size <- dots[["dot.size"]] %||% c(2, 6)
  label_size <- dots[["label.size"]] %||% 3
  do_label <- if (is.null(dots[["do.label"]])) TRUE else isTRUE(dots[["do.label"]])
  show_legend <- if (is.null(dots[["show.legend"]])) TRUE else isTRUE(dots[["show.legend"]])
  show_axes <- if (is.null(dots[["show.axes"]])) TRUE else isTRUE(dots[["show.axes"]])
  xlabel <- dots[["xlabel"]] %||% "Outgoing interaction strength"
  ylabel <- dots[["ylabel"]] %||% "Incoming interaction strength"
  weight_minmax <- dots[["weight.MinMax"]] %||% NULL

  plot_df <- ccc_scatter_df(pair_df)
  if (is.null(plot_df) || nrow(plot_df) == 0L) {
    log_message(
      "No aggregated sender-receiver interactions are available for role scatter plotting",
      message_type = "error"
    )
  }

  cell_levels <- unique(as.character(plot_df$cell))
  cell_cols <- palette_colors(
    cell_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  plot_df$cell <- factor(plot_df$cell, levels = cell_levels)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = outgoing,
      y = incoming,
      size = total_links,
      fill = cell
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      color = "grey20",
      stroke = 0.8,
      alpha = dot_alpha
    ) +
    ggplot2::scale_fill_manual(values = cell_cols[cell_levels], drop = FALSE) +
    ggplot2::scale_size_continuous(
      range = dot_size,
      limits = weight_minmax
    ) +
    ggplot2::labs(
      x = xlabel,
      y = ylabel,
      fill = "Cell type",
      size = "Link count"
    )
  if (isTRUE(do_label)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = cell),
      show.legend = FALSE,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.size = 0.2,
      segment.alpha = 0.5,
      seed = 1,
      size = label_size,
      color = "grey15"
    )
  }

  p <- p +
    ggplot2::labs(title = title, subtitle = subtitle) +
    cc_theme(theme_use = theme_use, theme_args = theme_args) +
    ggplot2::theme(
      text = ggplot2::element_text(size = font.size),
      legend.key.height = grid::unit(0.15, "in"),
      plot.title = ggplot2::element_text(size = font.size, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = font.size),
      legend.position = legend.position,
      legend.direction = legend.direction,
      axis.line.x = ggplot2::element_line(linewidth = 0.25),
      axis.line.y = ggplot2::element_line(linewidth = 0.25)
    )
  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!isTRUE(show_axes)) {
    p <- p + ggplot2::theme_void(base_size = font.size)
  }
  p
}

ccc_circle_plot <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "aggregation",
  top_n = 20,
  value = "weight",
  edge_threshold = 0,
  edge_size = c(0.5, 1.8),
  node_size = 5,
  node_alpha = 0.9,
  link_alpha = 0.6,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  link_palette = "Dark2",
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  label = FALSE,
  label.size = 4,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1
) {
  check_r("igraph", verbose = FALSE)
  value_col <- ccc_circle_value_col(value)
  plot_df <- ccc_network_df(
    pair_df = pair_df,
    interaction_df = interaction_df,
    display_by = display_by,
    top_n = top_n,
    value_col = value_col,
    edge_threshold = edge_threshold
  )
  if (is.null(plot_df) || nrow(plot_df) == 0L) {
    log_message(
      "No CCC records are available for circle plotting",
      message_type = "error"
    )
  }

  plot_df <- plot_df[is.finite(plot_df$weight), , drop = FALSE]
  if (nrow(plot_df) == 0L) {
    log_message(
      "No finite CCC edge values are available for circle plotting",
      message_type = "error"
    )
  }

  node_levels <- unique(c(
    as.character(plot_df$sender),
    as.character(plot_df$receiver)
  ))
  node_levels <- node_levels[!is.na(node_levels) & nzchar(node_levels)]
  net <- matrix(
    0,
    nrow = length(node_levels),
    ncol = length(node_levels),
    dimnames = list(node_levels, node_levels)
  )
  snd <- as.character(plot_df$sender)
  rcv <- as.character(plot_df$receiver)
  val <- as.numeric(plot_df$weight)
  snd_idx <- match(snd, node_levels)
  rcv_idx <- match(rcv, node_levels)
  keep <- !is.na(snd_idx) & !is.na(rcv_idx) & is.finite(val)
  if (any(keep)) {
    net[cbind(snd_idx[keep], rcv_idx[keep])] <- val[keep]
  }

  cell_cols <- palette_colors(
    node_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  link_cols <- palette_colors(
    unique(as.character(plot_df$sender)),
    palette = link_palette,
    palcolor = link_palcolor,
    NA_keep = TRUE
  )
  vertex_weight <- rowSums(net, na.rm = TRUE) + colSums(net, na.rm = TRUE)
  if (!any(is.finite(vertex_weight)) || max(vertex_weight, na.rm = TRUE) <= 0) {
    vertex_weight <- rep(1, length(node_levels))
  }
  vertex_size_max <- if (length(unique(vertex_weight)) == 1L) 5 else 15
  vertex_weight_max <- max(vertex_weight, na.rm = TRUE)
  vertex_size <- vertex_weight / vertex_weight_max * vertex_size_max + 5

  g <- igraph::graph_from_adjacency_matrix(
    net,
    mode = "directed",
    weighted = TRUE,
    diag = TRUE
  )
  edge_start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- igraph::layout_in_circle(g)
  coords_scale <- if (nrow(coords) != 1L) scale(coords) else coords
  loop_angle <- ifelse(
    coords_scale[igraph::V(g), 1] > 0,
    -atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1]),
    pi - atan(coords_scale[igraph::V(g), 2] / coords_scale[igraph::V(g), 1])
  )

  igraph::V(g)$size <- vertex_size
  igraph::V(g)$color <- cell_cols[igraph::V(g)$name]
  igraph::V(g)$frame.color <- cell_cols[igraph::V(g)$name]
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$label.cex <- max(0.8, font.size / 10)

  edge_weight_max <- max(igraph::E(g)$weight, na.rm = TRUE)
  if (!is.finite(edge_weight_max) || edge_weight_max <= 0) {
    edge_weight_max <- 1
  }
  edge_width_max <- max(edge_size) * 4
  igraph::E(g)$width <- 0.3 + igraph::E(g)$weight / edge_weight_max * edge_width_max
  igraph::E(g)$color <- grDevices::adjustcolor(
    cell_cols[igraph::V(g)$name[edge_start[, 1]]],
    alpha.f = link_alpha
  )
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  if (sum(edge_start[, 2] == edge_start[, 1]) != 0) {
    loop_idx <- which(edge_start[, 2] == edge_start[, 1])
    igraph::E(g)$loop.angle[loop_idx] <- loop_angle[edge_start[loop_idx, 1]]
  }

  radian_rescale <- function(x, start = 0, direction = 1) {
    rotate <- function(y) (y + start) %% (2 * pi) * direction
    rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label_locs <- radian_rescale(
    x = seq_len(length(igraph::V(g))),
    direction = -1,
    start = 0
  )
  label_dist <- vertex_size / max(vertex_size) + 2

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(old_par), silent = TRUE), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, if (is.null(title)) 0.5 else 2, 0.5))
  plot(
    g,
    edge.curved = 0.2,
    vertex.shape = "circle",
    layout = coords_scale,
    margin = 0.2,
    vertex.label.dist = label_dist,
    vertex.label.degree = label_locs,
    vertex.label.family = "Helvetica",
    edge.label.family = "Helvetica"
  )
  if (!is.null(title) || !is.null(subtitle)) {
    graphics::title(main = title, sub = subtitle)
  }
  grDevices::recordPlot()
}

ccc_network_df <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "aggregation",
  top_n = 20,
  value_col = "sum",
  edge_threshold = 0
) {
  edge_df <- if (identical(display_by, "interaction")) {
    tmp <- top_interactions(
      interaction_df,
      top_n = top_n,
      value_col = "score"
    )
    if (is.null(tmp) || nrow(tmp) == 0L) {
      data.frame()
    } else {
      stats::aggregate(
        score ~ sender + receiver,
        data = tmp,
        FUN = sum,
        na.rm = TRUE
      )
    }
  } else {
    pair_df
  }
  if (is.null(edge_df) || nrow(edge_df) == 0L) {
    return(data.frame())
  }
  if ("score" %in% colnames(edge_df) && !value_col %in% colnames(edge_df)) {
    edge_df[[value_col]] <- edge_df$score
  }
  edge_df <- edge_df[edge_df[[value_col]] >= edge_threshold, , drop = FALSE]
  edge_df$weight <- edge_df[[value_col]]
  edge_df
}

ccc_heatmap_plot <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "aggregation",
  top_n = 20,
  edge_value = "sum",
  color.by = "score",
  x_text_angle = 90,
  facet_by = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  title = NULL,
  subtitle = NULL,
  value_palette = "RdBu",
  value_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  fill_cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 9
  )
  if (identical(display_by, "interaction")) {
    plot_df <- top_interactions(
      interaction_df,
      top_n = top_n,
      value_col = "score"
    )
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No interaction-level CCC records are available for heatmap plotting",
        message_type = "error"
      )
    }
    sc <- scale_var(plot_df, color.by = color.by, agg_value = edge_value)
    if (is.null(facet_by)) {
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = interaction_label,
          y = pair,
          fill = .data[[sc$var]]
        )
      ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::labs(x = "Interaction", y = "Pair", fill = sc$label)
    } else if (identical(facet_by, "sender")) {
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = interaction_label,
          y = receiver,
          fill = .data[[sc$var]]
        )
      ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::facet_wrap(~sender, scales = "free_y") +
        ggplot2::labs(x = "Interaction", y = "Receiver", fill = sc$label)
    } else {
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = interaction_label,
          y = sender,
          fill = .data[[sc$var]]
        )
      ) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::facet_wrap(~receiver, scales = "free_y") +
        ggplot2::labs(x = "Interaction", y = "Sender", fill = sc$label)
    }
  } else {
    plot_df <- top_pairs(pair_df, top_n = top_n, value_col = edge_value)
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No aggregated sender-receiver interactions are available for heatmap plotting",
        message_type = "error"
      )
    }
    sc <- scale_var(plot_df, color.by = color.by, agg_value = edge_value)
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = sender,
        y = receiver,
        fill = .data[[sc$var]]
      )
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(
        ggplot2::aes(label = round(.data[[edge_value]], 2)),
        size = 3
      ) +
      ggplot2::labs(x = "Sender", y = "Receiver", fill = sc$label)
  }
  p <- p +
    ggplot2::scale_fill_gradientn(colours = fill_cols) +
    ggplot2::theme(
      axis.text.x = if (isTRUE(show_column_names)) {
        ggplot2::element_text(angle = x_text_angle, hjust = 1)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = if (isTRUE(show_row_names)) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      }
    )
  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}

ccc_flow_plot_df <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "interaction",
  top_n = 20,
  edge_value = "sum",
  edge_threshold = 0
) {
  if (identical(display_by, "interaction")) {
    plot_df <- top_interactions(
      interaction_df,
      top_n = top_n,
      value_col = "score"
    )
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      return(list(
        nodes = data.frame(),
        edges = data.frame(),
        breaks = c(1, 2, 3),
        labels = c("Sender", "Interaction", "Receiver")
      ))
    }
    plot_df <- plot_df[
      is.finite(plot_df$score) & plot_df$score > edge_threshold, ,
      drop = FALSE
    ]
    if (nrow(plot_df) == 0L) {
      return(list(
        nodes = data.frame(),
        edges = data.frame(),
        breaks = c(1, 2, 3),
        labels = c("Sender", "Interaction", "Receiver")
      ))
    }
    plot_df$weight <- plot_df$score

    branch_limit <- max(3L, min(5L, ceiling(sqrt(max(as.integer(top_n), 1L)))))
    sender_rank <- group_summary(
      df = plot_df,
      group_cols = c("interaction_label", "sender"),
      value_col = "weight",
      out_col = "weight",
      fun = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    receiver_rank <- group_summary(
      df = plot_df,
      group_cols = c("interaction_label", "receiver"),
      value_col = "weight",
      out_col = "weight",
      fun = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    sender_rank <- do.call(
      rbind,
      lapply(split(sender_rank, sender_rank$interaction_label), function(x) {
        x <- x[order(x$weight, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
        utils::head(x, branch_limit)
      })
    )
    receiver_rank <- do.call(
      rbind,
      lapply(split(receiver_rank, receiver_rank$interaction_label), function(x) {
        x <- x[order(x$weight, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
        utils::head(x, branch_limit)
      })
    )
    sender_keep <- paste(sender_rank$interaction_label, sender_rank$sender, sep = "\r")
    receiver_keep <- paste(receiver_rank$interaction_label, receiver_rank$receiver, sep = "\r")
    plot_df <- plot_df[
      paste(plot_df$interaction_label, plot_df$sender, sep = "\r") %in% sender_keep &
        paste(plot_df$interaction_label, plot_df$receiver, sep = "\r") %in% receiver_keep,
      ,
      drop = FALSE
    ]
    if (nrow(plot_df) == 0L) {
      return(list(
        nodes = data.frame(),
        edges = data.frame(),
        breaks = c(1, 2, 3),
        labels = c("Sender", "Interaction", "Receiver")
      ))
    }

    sender_levels <- ccc_rank_flow_nodes(plot_df$sender, plot_df$weight)
    interaction_levels <- ccc_rank_flow_nodes(
      plot_df$interaction_label,
      plot_df$weight
    )
    receiver_levels <- ccc_rank_flow_nodes(plot_df$receiver, plot_df$weight)

    sender_nodes <- ccc_make_flow_nodes(
      sender_levels,
      x = 1,
      column = "sender",
      plot_df = plot_df,
      label_col = "sender"
    )
    interaction_nodes <- ccc_make_flow_nodes(
      interaction_levels,
      x = 2,
      column = "interaction",
      plot_df = plot_df,
      label_col = "interaction_label"
    )
    receiver_nodes <- ccc_make_flow_nodes(
      receiver_levels,
      x = 3,
      column = "receiver",
      plot_df = plot_df,
      label_col = "receiver"
    )
    node_df <- rbind(sender_nodes, interaction_nodes, receiver_nodes)

    left_edges <- group_summary(
      df = plot_df,
      group_cols = c("sender", "interaction_label"),
      value_col = "weight",
      out_col = "weight",
      fun = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    left_edges$edge_id <- paste0("L", seq_len(nrow(left_edges)))
    left_edges$from_id <- paste0("sender::", left_edges$sender)
    left_edges$to_id <- paste0("interaction::", left_edges$interaction_label)
    left_edges$edge_group <- left_edges$sender
    left_edges$edge_label <- left_edges$interaction_label
    left_edges <- left_edges[, c(
      "edge_id",
      "from_id",
      "to_id",
      "weight",
      "edge_group",
      "edge_label"
    ), drop = FALSE]

    right_edges <- group_summary(
      df = plot_df,
      group_cols = c("interaction_label", "receiver"),
      value_col = "weight",
      out_col = "weight",
      fun = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    sender_lookup <- group_summary(
      df = plot_df,
      group_cols = c("interaction_label", "receiver", "sender"),
      value_col = "weight",
      out_col = "weight",
      fun = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    sender_lookup <- sender_lookup[
      order(sender_lookup$interaction_label, sender_lookup$receiver, -sender_lookup$weight),
      ,
      drop = FALSE
    ]
    sender_lookup <- sender_lookup[!duplicated(paste(sender_lookup$interaction_label, sender_lookup$receiver, sep = "\r")), , drop = FALSE]
    sender_key <- paste(sender_lookup$interaction_label, sender_lookup$receiver, sep = "\r")
    right_key <- paste(right_edges$interaction_label, right_edges$receiver, sep = "\r")
    right_edges$edge_id <- paste0("R", seq_len(nrow(right_edges)))
    right_edges$from_id <- paste0("interaction::", right_edges$interaction_label)
    right_edges$to_id <- paste0("receiver::", right_edges$receiver)
    right_edges$edge_group <- sender_lookup$sender[match(right_key, sender_key)]
    right_edges$edge_group[is.na(right_edges$edge_group)] <- right_edges$receiver[is.na(right_edges$edge_group)]
    right_edges$edge_label <- right_edges$interaction_label
    right_edges <- right_edges[, c(
      "edge_id",
      "from_id",
      "to_id",
      "weight",
      "edge_group",
      "edge_label"
    ), drop = FALSE]
    edge_df <- rbind(left_edges, right_edges)
    breaks <- c(1, 2, 3)
    labels <- c("Sender", "Interaction", "Receiver")
  } else {
    plot_df <- top_pairs(pair_df, top_n = top_n, value_col = edge_value)
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      return(list(
        nodes = data.frame(),
        edges = data.frame(),
        breaks = c(1, 2),
        labels = c("Sender", "Receiver")
      ))
    }
    plot_df <- plot_df[
      is.finite(plot_df[[edge_value]]) & plot_df[[edge_value]] > edge_threshold, ,
      drop = FALSE
    ]
    if (nrow(plot_df) == 0L) {
      return(list(
        nodes = data.frame(),
        edges = data.frame(),
        breaks = c(1, 2),
        labels = c("Sender", "Receiver")
      ))
    }
    plot_df$weight <- plot_df[[edge_value]]

    sender_levels <- ccc_rank_flow_nodes(plot_df$sender, plot_df$weight)
    receiver_levels <- ccc_rank_flow_nodes(plot_df$receiver, plot_df$weight)
    sender_nodes <- ccc_make_flow_nodes(
      sender_levels,
      x = 1,
      column = "sender",
      plot_df = plot_df,
      label_col = "sender"
    )
    receiver_nodes <- ccc_make_flow_nodes(
      receiver_levels,
      x = 2,
      column = "receiver",
      plot_df = plot_df,
      label_col = "receiver"
    )
    node_df <- rbind(sender_nodes, receiver_nodes)
    edge_df <- data.frame(
      edge_id = paste0("A", seq_len(nrow(plot_df))),
      from_id = paste0("sender::", plot_df$sender),
      to_id = paste0("receiver::", plot_df$receiver),
      weight = plot_df$weight,
      edge_group = plot_df$sender,
      edge_label = plot_df$pair,
      stringsAsFactors = FALSE
    )
    breaks <- c(1, 2)
    labels <- c("Sender", "Receiver")
  }

  if (nrow(node_df) == 0L || nrow(edge_df) == 0L) {
    return(list(
      nodes = data.frame(),
      edges = data.frame(),
      breaks = breaks,
      labels = labels
    ))
  }

  edge_df <- merge(
    edge_df,
    node_df[, c("node_id", "x", "y")],
    by.x = "from_id",
    by.y = "node_id",
    all.x = TRUE
  )
  edge_df <- merge(
    edge_df,
    node_df[, c("node_id", "x", "y")],
    by.x = "to_id",
    by.y = "node_id",
    all.x = TRUE,
    suffixes = c("_from", "_to")
  )
  edge_df <- edge_df[
    !is.na(edge_df$x_from) & !is.na(edge_df$x_to), ,
    drop = FALSE
  ]
  list(nodes = node_df, edges = edge_df, breaks = breaks, labels = labels)
}

ccc_rank_flow_nodes <- function(labels, weights) {
  labels <- as.character(labels)
  keep <- !is.na(labels) & nzchar(labels) & is.finite(weights)
  labels <- labels[keep]
  weights <- weights[keep]
  if (length(labels) == 0L) {
    return(character(0))
  }
  ord <- stats::aggregate(
    x = weights,
    by = list(label = labels),
    FUN = sum,
    na.rm = TRUE
  )
  ord$label[order(ord$x, decreasing = TRUE)]
}

ccc_make_flow_nodes <- function(levels, x, column, plot_df, label_col) {
  if (length(levels) == 0L) {
    return(data.frame())
  }
  agg <- stats::aggregate(
    x = plot_df$weight,
    by = list(label = plot_df[[label_col]]),
    FUN = sum,
    na.rm = TRUE
  )
  agg <- agg[match(levels, agg$label), , drop = FALSE]
  data.frame(
    node_id = paste0(column, "::", levels),
    label = levels,
    column = column,
    x = x,
    y = rev(seq_along(levels)),
    weight = agg$x,
    stringsAsFactors = FALSE
  )
}

ccc_sigmoid_curve_df <- function(edge_df, n = 80) {
  if (is.null(edge_df) || nrow(edge_df) == 0L) {
    return(data.frame())
  }
  pieces <- lapply(seq_len(nrow(edge_df)), function(i) {
    row <- edge_df[i, , drop = FALSE]
    t <- seq(0, 1, length.out = n)
    s <- stats::plogis(seq(-6, 6, length.out = n))
    s <- (s - min(s)) / diff(range(s))
    data.frame(
      edge_id = row$edge_id,
      x = row$x_from + (row$x_to - row$x_from) * t,
      y = row$y_from + (row$y_to - row$y_from) * s,
      weight = row$weight,
      edge_group = row$edge_group,
      edge_colour = row$edge_colour,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

ccc_network_label_layer <- function(
  data,
  mapping,
  label = FALSE,
  label_size = 4,
  label_fg = "white",
  label_bg = "black",
  label_bg_r = 0.1,
  repel = FALSE,
  direction = "both",
  seed = 11,
  ...
) {
  if (is.null(data) || nrow(data) == 0L) {
    return(NULL)
  }
  args <- c(
    list(
      data = data,
      mapping = mapping,
      size = label_size,
      inherit.aes = FALSE,
      show.legend = FALSE
    ),
    list(...)
  )
  if (isTRUE(label) || isTRUE(repel)) {
    args <- c(args, list(
      color = label_fg,
      bg.color = label_bg,
      bg.r = label_bg_r,
      box.padding = 0.2,
      point.padding = 0.15,
      segment.alpha = 0,
      segment.color = NA,
      min.segment.length = 0,
      direction = direction,
      seed = seed
    ))
    return(do.call(ggrepel::geom_text_repel, args))
  }
  args$color <- "grey15"
  do.call(ggplot2::geom_text, args)
}

ccc_flow_network_plot <- function(
  pair_df,
  interaction_df = NULL,
  plot_type = c("arrow", "sigmoid"),
  display_by = "interaction",
  top_n = 20,
  edge_value = "sum",
  edge_threshold = 0,
  edge_size = c(0.2, 1),
  edge_color = NULL,
  link_curvature = 0.2,
  link_alpha = 0.6,
  directed = FALSE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 6,
  node_alpha = 0.9,
  title = NULL,
  subtitle = NULL,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  label = FALSE,
  label.size = 4,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1,
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  plot_type <- match.arg(plot_type)
  flow <- ccc_flow_plot_df(
    pair_df = pair_df,
    interaction_df = interaction_df,
    display_by = display_by,
    top_n = top_n,
    edge_value = edge_value,
    edge_threshold = edge_threshold
  )
  node_df <- flow$nodes
  edge_df <- flow$edges
  if (
    is.null(node_df) ||
      nrow(node_df) == 0L ||
      is.null(edge_df) ||
      nrow(edge_df) == 0L
  ) {
    log_message(
      paste0(
        "No ",
        if (identical(display_by, "interaction")) {
          "interaction-level"
        } else {
          "aggregated"
        },
        " CCC records are available for ",
        plot_type,
        " plotting"
      ),
      message_type = "error"
    )
  }

  celltype_levels <- unique(node_df$label[
    node_df$column %in% c("sender", "receiver")
  ])
  celltype_cols <- palette_colors(
    celltype_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  node_df$fill <- ifelse(
    node_df$column %in% c("sender", "receiver"),
    unname(celltype_cols[node_df$label]),
    "white"
  )
  node_df$border <- ifelse(
    node_df$column == "interaction",
    "grey55",
    unname(celltype_cols[node_df$label])
  )
  node_df$text_colour <- ifelse(
    node_df$column == "interaction",
    "grey15",
    "grey10"
  )
  if (
    length(stats::na.omit(node_df$weight)) <= 1L ||
      diff(range(node_df$weight, na.rm = TRUE)) == 0
  ) {
    node_df$size_scaled <- node_size
  } else {
    node_df$size_scaled <- scales::rescale(
      node_df$weight,
      to = c(node_size * 0.8, node_size * 1.5),
      from = range(node_df$weight, na.rm = TRUE)
    )
  }
  if (all(!is.finite(node_df$size_scaled))) {
    node_df$size_scaled <- node_size
  }

  if (is.null(edge_color) || length(edge_color) == 0L) {
    edge_df$edge_colour <- unname(celltype_cols[edge_df$edge_group])
  } else if (!is.null(names(edge_color)) && any(nzchar(names(edge_color)))) {
    edge_df$edge_colour <- unname(edge_color[as.character(edge_df$edge_group)])
  } else if (length(edge_color) == 1L) {
    edge_df$edge_colour <- rep(edge_color, nrow(edge_df))
  } else if (length(edge_color) == nrow(edge_df)) {
    edge_df$edge_colour <- as.character(edge_color)
  } else if (length(edge_color) == length(unique(edge_df$edge_group))) {
    edge_map <- stats::setNames(
      as.character(edge_color),
      unique(as.character(edge_df$edge_group))
    )
    edge_df$edge_colour <- unname(edge_map[as.character(edge_df$edge_group)])
  } else {
    log_message(
      paste0(
        "{.arg edge_color} must be length 1, match the number of edges, ",
        "or be a named vector keyed by sender group"
      ),
      message_type = "error"
    )
  }
  edge_df$edge_colour[is.na(edge_df$edge_colour)] <- "#6C757D"
  edge_df$curvature <- ifelse(
    edge_df$x_from < edge_df$x_to,
    link_curvature,
    -link_curvature
  )

  label_sender <- node_df[node_df$column == "sender", , drop = FALSE]
  label_interaction <- node_df[node_df$column == "interaction", , drop = FALSE]
  label_receiver <- node_df[node_df$column == "receiver", , drop = FALSE]
  cell_node_df <- node_df[node_df$column %in% c("sender", "receiver"), , drop = FALSE]
  interaction_node_df <- node_df[node_df$column == "interaction", , drop = FALSE]
  label_sender$x_label <- label_sender$x - 0.06
  label_receiver$x_label <- label_receiver$x + 0.06
  label_interaction$y_label <- label_interaction$y + 0.08

  p <- ggplot2::ggplot()
  if (identical(plot_type, "sigmoid")) {
    curve_df <- ccc_sigmoid_curve_df(edge_df)
    p <- p +
      ggplot2::geom_path(
        data = curve_df,
        ggplot2::aes(
          x = x,
          y = y,
          group = edge_id,
          linewidth = weight,
          color = edge_colour
        ),
        alpha = link_alpha,
        lineend = "round",
        show.legend = FALSE
      )
  } else if (identical(plot_type, "arrow")) {
    p <- p +
      ggplot2::geom_segment(
        data = edge_df,
        ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = weight,
          color = edge_colour
        ),
        alpha = link_alpha,
        lineend = "round",
        arrow = NULL,
        show.legend = FALSE
      )
  } else {
    p <- p +
      ggplot2::geom_curve(
        data = edge_df,
        ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = weight,
          color = edge_colour
        ),
        curvature = if (identical(display_by, "interaction")) {
          0.12
        } else {
          link_curvature
        },
        alpha = link_alpha,
        lineend = "round",
        arrow = if (isTRUE(directed) || identical(display_by, "interaction")) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE
      )
  }

  p <- p +
    ggplot2::scale_color_identity() +
    ggplot2::geom_point(
      data = interaction_node_df,
      ggplot2::aes(
        x = x,
        y = y,
        size = size_scaled
      ),
      shape = 22,
      fill = interaction_node_df$fill,
      color = interaction_node_df$border,
      stroke = 0.9,
      alpha = node_alpha,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = cell_node_df,
      ggplot2::aes(
        x = x,
        y = y,
        fill = label,
        size = size_scaled
      ),
      shape = 21,
      color = "grey20",
      stroke = 0.8,
      alpha = node_alpha,
      show.legend = FALSE
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_linewidth_continuous(range = edge_size, guide = "none") +
    ccc_network_label_layer(
      data = label_sender,
      mapping = ggplot2::aes(x = x_label, y = y, label = label),
      label = label,
      label_size = label.size,
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = TRUE,
      direction = "y",
      hjust = 1
    ) +
    ccc_network_label_layer(
      data = label_interaction,
      mapping = ggplot2::aes(x = x, y = y_label, label = label),
      label = label,
      label_size = label.size,
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = isTRUE(label) && identical(display_by, "interaction"),
      direction = "y",
      hjust = 0.5,
      vjust = 0,
      lineheight = 0.9
    ) +
    ccc_network_label_layer(
      data = label_receiver,
      mapping = ggplot2::aes(x = x_label, y = y, label = label),
      label = label,
      label_size = label.size,
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = TRUE,
      direction = "y",
      hjust = 0
    ) +
    ggplot2::geom_point(
      data = data.frame(
        x = NA_real_,
        y = NA_real_,
        label = celltype_levels,
        stringsAsFactors = FALSE
      ),
      ggplot2::aes(x = x, y = y, fill = label),
      shape = 21,
      size = max(3, node_size * 0.75),
      color = "grey20",
      stroke = 0.8,
      show.legend = TRUE,
      na.rm = TRUE,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = celltype_cols,
      name = legend.title %||% "Cell type",
      breaks = names(celltype_cols),
      guide = ggplot2::guide_legend(
        override.aes = list(shape = 21, size = 4, color = "grey20", stroke = 0.8)
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = flow$breaks,
      labels = flow$labels,
      expand = ggplot2::expansion(mult = c(0.18, 0.18))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL)

  p <- finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )

  p +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 28, 5.5, 28)
    )
}

ccc_dot_plot <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "aggregation",
  top_n = 20,
  edge_value = "sum",
  color.by = "score",
  x_text_angle = 90,
  facet_by = NULL,
  title = NULL,
  subtitle = NULL,
  value_palette = "RdBu",
  value_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (identical(display_by, "interaction")) {
    plot_df <- top_interactions(
      interaction_df,
      top_n = top_n,
      value_col = "score"
    )
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No interaction-level CCC records are available for dot plotting",
        message_type = "error"
      )
    }
    sc <- scale_var(plot_df, color.by = color.by, agg_value = edge_value)
    if (is.null(facet_by)) {
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = interaction_label,
          y = pair,
          size = score,
          fill = .data[[sc$var]]
        )
      ) +
        ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.2) +
        ggplot2::labs(
          x = "Interaction",
          y = "Pair",
          size = "score",
          fill = sc$label
        )
    } else if (identical(facet_by, "sender")) {
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = interaction_label,
          y = receiver,
          size = score,
          fill = .data[[sc$var]]
        )
      ) +
        ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.2) +
        ggplot2::facet_wrap(~sender, scales = "free_y") +
        ggplot2::labs(
          x = "Interaction",
          y = "Receiver",
          size = "score",
          fill = sc$label
        )
    } else {
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = interaction_label,
          y = sender,
          size = score,
          fill = .data[[sc$var]]
        )
      ) +
        ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.2) +
        ggplot2::facet_wrap(~receiver, scales = "free_y") +
        ggplot2::labs(
          x = "Interaction",
          y = "Sender",
          size = "score",
          fill = sc$label
        )
    }
  } else {
    plot_df <- top_pairs(pair_df, top_n = top_n, value_col = edge_value)
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No aggregated sender-receiver interactions are available for dot plotting",
        message_type = "error"
      )
    }
    sc <- scale_var(plot_df, color.by = color.by, agg_value = edge_value)
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = sender,
        y = receiver,
        size = count,
        fill = .data[[sc$var]]
      )
    ) +
      ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.2) +
      ggplot2::labs(
        x = "Sender",
        y = "Receiver",
        size = "count",
        fill = sc$label
      ) +
      ggplot2::coord_equal()
  }
  fill_cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 9
  )
  p <- p +
    ggplot2::scale_fill_gradientn(colours = fill_cols) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = 1)
    )
  p <- p +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = font.size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = font.size * 1.15,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(size = font.size),
      legend.position = legend.position,
      legend.direction = legend.direction,
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_blank()
    )
  p
}

ccc_source_target_dot_plot <- function(
  interaction_df,
  top_n = 20,
  color.by = "score",
  interaction.use = NULL,
  pairLR.use = NULL,
  ncols = NULL,
  nrows = NULL,
  title = NULL,
  subtitle = NULL,
  value_palette = "Reds",
  value_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (is.null(interaction_df) || nrow(interaction_df) == 0L) {
    log_message(
      "No interaction-level CCC records are available for dot plotting",
      message_type = "error"
    )
  }
  score_col <- if (identical(color.by, "pvalue")) {
    c("neglog10_pvalue", "specificity")[
      c("neglog10_pvalue", "specificity") %in% colnames(interaction_df)
    ][1]
  } else {
    "score"
  }
  score_col <- score_col %||% "score"
  top_n_use <- if (!is.null(interaction.use) || !is.null(pairLR.use)) {
    NULL
  } else {
    top_n
  }
  plot_df <- top_interactions(
    interaction_df,
    top_n = top_n_use,
    value_col = score_col
  )
  plot_df <- plot_df[
    is.finite(suppressWarnings(as.numeric(plot_df$score))) &
      suppressWarnings(as.numeric(plot_df$score)) > 0,
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0L) {
    log_message(
      "No dotplot records remain after filtering",
      message_type = "error"
    )
  }
  plot_df$target <- as.character(plot_df$receiver)
  plot_df$source <- as.character(plot_df$sender)
  plot_df$interaction_plot <- as.character(plot_df$interaction_label)
  metric <- suppressWarnings(as.numeric(plot_df[[score_col]]))
  metric[!is.finite(metric)] <- 0
  plot_df$plot_metric <- metric

  source_order <- names(sort(tapply(plot_df$plot_metric, plot_df$source, sum), decreasing = TRUE))
  target_order <- names(sort(tapply(plot_df$plot_metric, plot_df$target, sum), decreasing = TRUE))
  interaction_order <- names(sort(tapply(plot_df$plot_metric, plot_df$interaction_plot, sum), decreasing = FALSE))
  plot_df$source <- factor(plot_df$source, levels = source_order)
  plot_df$target <- factor(plot_df$target, levels = target_order)
  plot_df$interaction_plot <- factor(plot_df$interaction_plot, levels = interaction_order)

  grid <- ccc_panel_grid(
    n_panels = length(source_order),
    ncols = ncols,
    nrows = nrows,
    context = "source-target dot plot"
  )

  fill_cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 9
  )
  fill_label <- if (identical(color.by, "pvalue")) {
    "-log10(p-value)"
  } else {
    "Communication score"
  }
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = target,
      y = interaction_plot,
      size = score,
      fill = plot_metric
    )
  ) +
    ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.25, alpha = 0.9) +
    ggplot2::facet_wrap(
      ~source,
      ncol = grid$ncols,
      nrow = grid$nrows,
      scales = "free_x"
    ) +
    ggplot2::scale_fill_gradientn(colours = fill_cols, name = fill_label) +
    ggplot2::scale_size_area(name = "Communication score", max_size = 8) +
    ggplot2::labs(x = "Target", y = "Interactions (Ligand -> Receptor)")

  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

ccc_source_target_tile_plot <- function(
  interaction_df,
  top_n = 20,
  color.by = "score",
  interaction.use = NULL,
  pairLR.use = NULL,
  title = NULL,
  subtitle = NULL,
  value_palette = "Reds",
  value_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (is.null(interaction_df) || nrow(interaction_df) == 0L) {
    log_message(
      "No interaction-level CCC records are available for tile plotting",
      message_type = "error"
    )
  }
  score_col <- if (identical(color.by, "pvalue")) {
    c("neglog10_pvalue", "specificity")[
      c("neglog10_pvalue", "specificity") %in% colnames(interaction_df)
    ][1]
  } else {
    "score"
  }
  score_col <- score_col %||% "score"
  top_n_use <- if (!is.null(interaction.use) || !is.null(pairLR.use)) {
    NULL
  } else {
    top_n
  }
  data <- top_interactions(interaction_df, top_n = top_n_use, value_col = score_col)
  if (nrow(data) == 0L) {
    log_message("No tile plot records remain after filtering", message_type = "error")
  }
  data$ligand_display <- ifelse(
    is.na(data$ligand) | !nzchar(data$ligand),
    data$interaction_label,
    data$ligand
  )
  data$receptor_display <- ifelse(
    is.na(data$receptor) | !nzchar(data$receptor),
    data$interaction_label,
    data$receptor
  )
  metric <- suppressWarnings(as.numeric(data[[score_col]]))
  metric[!is.finite(metric)] <- 0
  data$plot_metric <- metric
  interaction_order <- names(sort(tapply(data$plot_metric, data$interaction_label, sum), decreasing = FALSE))

  left <- data.frame(
    type = "Source ligand",
    cell_type = data$sender,
    gene = data$ligand_display,
    interaction_label = data$interaction_label,
    plot_metric = data$plot_metric,
    stringsAsFactors = FALSE
  )
  right <- data.frame(
    type = "Target receptor",
    cell_type = data$receiver,
    gene = data$receptor_display,
    interaction_label = data$interaction_label,
    plot_metric = data$plot_metric,
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(left, right)
  plot_df$interaction_label <- factor(plot_df$interaction_label, levels = interaction_order)
  fill_cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 9
  )
  fill_label <- if (identical(color.by, "pvalue")) {
    "-log10(p-value)"
  } else {
    "Communication score"
  }
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = cell_type, y = interaction_label, fill = plot_metric)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::facet_wrap(~type, nrow = 1, scales = "free_x") +
    ggplot2::scale_fill_gradientn(colours = fill_cols, name = fill_label) +
    ggplot2::labs(x = "Cell type", y = NULL)

  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
}

ccc_sample_dot_plot <- function(
  df,
  sample_key = "sample",
  top_n = 20,
  top_n_pairs = NULL,
  color.by = "score",
  interaction.use = NULL,
  pairLR.use = NULL,
  ncols = NULL,
  nrows = NULL,
  title = NULL,
  subtitle = NULL,
  value_palette = "Reds",
  value_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (is.null(df) || nrow(df) == 0L) {
    log_message("No communication records are available for sample dot plotting", message_type = "error")
  }
  sample_col <- if (sample_key %in% colnames(df)) {
    sample_key
  } else if ("dataset" %in% colnames(df)) {
    "dataset"
  } else {
    NULL
  }
  if (is.null(sample_col)) {
    log_message(
      "{.arg sample_key} was not found in the CCC result table",
      message_type = "error"
    )
  }
  df <- prepare_plot_df(df)
  score_col <- if (identical(color.by, "pvalue")) {
    c("neglog10_pvalue", "specificity")[
      c("neglog10_pvalue", "specificity") %in% colnames(df)
    ][1]
  } else {
    "score"
  }
  score_col <- score_col %||% "score"
  score_by_interaction <- sort(
    tapply(df[[score_col]], df$interaction_label, sum, na.rm = TRUE),
    decreasing = TRUE
  )
  if (!is.null(interaction.use) || !is.null(pairLR.use) || is.null(top_n) || !is.numeric(top_n) || top_n <= 0L) {
    keep_interactions <- names(score_by_interaction)
  } else {
    keep_interactions <- names(score_by_interaction)[seq_len(min(top_n, length(score_by_interaction)))]
  }
  df <- df[df$interaction_label %in% keep_interactions, , drop = FALSE]
  df$sample <- as.character(df[[sample_col]])
  df$pair <- paste(df$sender, df$receiver, sep = " -> ")
  if (!is.null(top_n_pairs) && is.numeric(top_n_pairs) && top_n_pairs > 0L) {
    df$.ccc_metric <- suppressWarnings(as.numeric(df[[score_col]]))
    df$.ccc_metric[!is.finite(df$.ccc_metric)] <- 0
    df <- do.call(rbind, lapply(split(df, df$interaction_label), function(df_i) {
      pair_scores <- sort(tapply(df_i$.ccc_metric, df_i$pair, sum, na.rm = TRUE), decreasing = TRUE)
      keep_pairs <- names(pair_scores)[seq_len(min(top_n_pairs, length(pair_scores)))]
      df_i[df_i$pair %in% keep_pairs, , drop = FALSE]
    }))
    rownames(df) <- NULL
    df$.ccc_metric <- NULL
  }
  if (nrow(df) == 0L) {
    log_message("No sample dotplot records remain after filtering", message_type = "error")
  }
  df$interaction_label <- factor(df$interaction_label, levels = keep_interactions)
  grid <- ccc_panel_grid(
    n_panels = length(keep_interactions),
    ncols = ncols %||% min(3L, length(keep_interactions)),
    nrows = nrows,
    context = "sample dot plot"
  )
  fill_cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 9
  )
  fill_label <- if (identical(color.by, "pvalue")) {
    "-log10(p-value)"
  } else {
    "Communication score"
  }
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = sample, y = pair, size = score, fill = .data[[score_col]])
  ) +
    ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.25, alpha = 0.9) +
    ggplot2::facet_wrap(
      ~interaction_label,
      ncol = grid$ncols,
      nrow = grid$nrows,
      scales = "free_y"
    ) +
    ggplot2::scale_fill_gradientn(colours = fill_cols, name = fill_label) +
    ggplot2::scale_size_area(name = "Communication score", max_size = 8) +
    ggplot2::labs(x = sample_col, y = "Sender -> receiver")

  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

ccc_chord_plot <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "aggregation",
  top_n = 20,
  edge_value = "sum",
  edge_threshold = 0,
  link_alpha = 0.6,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  font.size = 10,
  reduce = TRUE,
  max.groups = 8,
  small.gap = 1,
  big.gap = 8,
  lab.cex = 0.6
) {
  check_r("circlize", verbose = FALSE)
  display_by <- match.arg(display_by, c("aggregation", "interaction"))
  cell_lab_cex <- max(lab.cex, 0.72)
  pair_alpha <- min(link_alpha, 0.5)

  pair_plot <- if (identical(display_by, "interaction")) {
    interaction_plot <- top_interactions(
      interaction_df %||% data.frame(),
      top_n = top_n,
      value_col = "score"
    )
    if (is.null(interaction_plot) || nrow(interaction_plot) == 0L) {
      log_message(
        "No interaction-level CCC records are available for chord plotting",
        message_type = "error"
      )
    }
    if (!"pair" %in% colnames(interaction_plot)) {
      interaction_plot$pair <- paste(
        interaction_plot$sender,
        interaction_plot$receiver,
        sep = " -> "
      )
    }
    stats::aggregate(
      score ~ sender + receiver + pair,
      data = interaction_plot,
      FUN = sum,
      na.rm = TRUE
    )
  } else {
    top_pairs(
      pair_df,
      top_n = top_n,
      value_col = edge_value
    )
  }

  if (is.null(pair_plot) || nrow(pair_plot) == 0L) {
    log_message(
      "No CCC records are available for chord plotting",
      message_type = "error"
    )
  }
  if (!"pair" %in% colnames(pair_plot)) {
    pair_plot$pair <- paste(pair_plot$sender, pair_plot$receiver, sep = " -> ")
  }
  pair_value_col <- if (identical(display_by, "interaction")) {
    "score"
  } else if (edge_value %in% colnames(pair_plot)) {
    edge_value
  } else {
    c("sum", "score", "count")[c("sum", "score", "count") %in% colnames(pair_plot)][1]
  }
  if (is.na(pair_value_col) || !nzchar(pair_value_col)) {
    log_message(
      "No CCC records are available for chord plotting",
      message_type = "error"
    )
  }

  pair_plot[[pair_value_col]] <- suppressWarnings(as.numeric(pair_plot[[pair_value_col]]))
  pair_plot <- pair_plot[
    is.finite(pair_plot[[pair_value_col]]) &
      pair_plot[[pair_value_col]] > 0 &
      pair_plot[[pair_value_col]] >= edge_threshold, ,
    drop = FALSE
  ]
  if (nrow(pair_plot) == 0L) {
    log_message(
      "No CCC records remain for chord plotting after filtering",
      message_type = "error"
    )
  }

  strength_df <- rbind(
    data.frame(
      cell = pair_plot$sender,
      strength = pair_plot[[pair_value_col]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      cell = pair_plot$receiver,
      strength = pair_plot[[pair_value_col]],
      stringsAsFactors = FALSE
    )
  )
  strength_df <- stats::aggregate(
    strength ~ cell,
    data = strength_df,
    FUN = sum,
    na.rm = TRUE
  )
  strength_df <- strength_df[
    order(strength_df$strength, decreasing = TRUE, na.last = TRUE), ,
    drop = FALSE
  ]
  if (
    isTRUE(reduce) &&
      is.numeric(max.groups) &&
      length(max.groups) == 1L &&
      is.finite(max.groups) &&
      nrow(strength_df) > max.groups
  ) {
    keep_cells <- utils::head(strength_df$cell, max.groups)
    pair_plot <- pair_plot[
      pair_plot$sender %in% keep_cells &
        pair_plot$receiver %in% keep_cells, ,
      drop = FALSE
    ]
    if (nrow(pair_plot) == 0L) {
      log_message(
        "No CCC records remain for chord plotting after cell-group reduction",
        message_type = "error"
      )
    }
    strength_df <- strength_df[strength_df$cell %in% keep_cells, , drop = FALSE]
  }

  cell_order <- unique(c(
    strength_df$cell,
    as.character(pair_plot$sender),
    as.character(pair_plot$receiver)
  ))
  cell_order <- cell_order[!is.na(cell_order) & nzchar(cell_order)]
  if (length(cell_order) == 0L) {
    log_message(
      "No cell groups are available for chord plotting",
      message_type = "error"
    )
  }

  cell_cols <- palette_colors(
    cell_order,
    palette = cell_palette,
    palcolor = cell_palcolor
  )
  names(cell_cols) <- cell_order

  chord_df <- pair_plot[, c("sender", "receiver", pair_value_col), drop = FALSE]
  colnames(chord_df) <- c("source", "target", "prob")
  edge_col <- unname(cell_cols[as.character(chord_df$source)])
  edge_col[is.na(edge_col)] <- "#7CAAB0"
  grid.col <- cell_cols[cell_order]
  names(grid.col) <- cell_order

  preallocate_height <- tryCatch(
    max(
      graphics::strwidth(
        cell_order,
        cex = max(cell_lab_cex, 0.5)
      ),
      na.rm = TRUE
    ),
    error = function(e) {
      # Example devices may not have initialized a base plot yet.
      0.045
    }
  )
  if (!is.finite(preallocate_height) || preallocate_height <= 0) {
    preallocate_height <- 0.045
  }
  preallocate_height <- min(max(preallocate_height, 0.045), 0.085)

  old_par <- graphics::par(no.readonly = TRUE)
  circlize::circos.clear()
  on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
  on.exit(try(graphics::par(old_par), silent = TRUE), add = TRUE)

  graphics::par(mar = c(1, 1, if (!is.null(title) || !is.null(subtitle)) 3 else 1, 1))
  circlize::circos.par(
    start.degree = 90,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01),
    points.overflow.warning = FALSE
  )
  circlize::chordDiagram(
    x = chord_df,
    order = cell_order,
    grid.col = grid.col,
    col = grDevices::adjustcolor(edge_col, alpha.f = pair_alpha),
    transparency = max(0, min(1, 1 - pair_alpha)),
    annotationTrack = "grid",
    annotationTrackHeight = c(0.03),
    preAllocateTracks = list(track.height = preallocate_height),
    small.gap = small.gap,
    big.gap = big.gap,
    link.sort = TRUE,
    link.decreasing = FALSE,
    link.largest.ontop = TRUE,
    directional = 1,
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow",
    link.border = NA,
    reduce = 1e-5
  )
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim_i <- circlize::get.cell.meta.data("xlim")
      ylim_i <- circlize::get.cell.meta.data("ylim")
      sector_i <- circlize::get.cell.meta.data("sector.index")
      circlize::circos.text(
        mean(xlim_i),
        ylim_i[1],
        sector_i,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = cell_lab_cex,
        col = cell_cols[sector_i]
      )
    },
    bg.border = NA
  )

  if (!is.null(title) || !is.null(subtitle)) {
    graphics::title(
      main = title,
      sub = subtitle,
      cex.main = font.size / 10 * 1.15,
      cex.sub = font.size / 10
    )
  }
  grDevices::recordPlot()
}

ccc_gene_chord_plot <- function(
  df,
  top_n = 20,
  edge_threshold = 0,
  link_alpha = 0.6,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  font.size = 10,
  small.gap = 1,
  big.gap = 8,
  lab.cex = 0.55
) {
  check_r("circlize", verbose = FALSE)
  df <- prepare_plot_df(df)
  df <- top_interactions(df, top_n = top_n, value_col = "score")
  score <- suppressWarnings(as.numeric(df$score))
  df <- df[
    is.finite(score) &
      score > 0 &
      score >= edge_threshold &
      !is.na(df$sender) &
      nzchar(df$sender) &
      !is.na(df$receiver) &
      nzchar(df$receiver),
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0L) {
    log_message(
      "No CCC records remain for gene-level chord plotting",
      message_type = "error"
    )
  }

  ligand <- ccc_display_gene(df$ligand)
  receptor <- ccc_display_gene(df$receptor)
  ligand[!nzchar(ligand)] <- as.character(df$interaction_label[!nzchar(ligand)])
  receptor[!nzchar(receptor)] <- as.character(df$interaction_label[!nzchar(receptor)])
  source_node <- paste(as.character(df$sender), ligand, sep = "\n")
  target_node <- paste(as.character(df$receiver), receptor, sep = "\n")
  edge_df <- data.frame(
    source = source_node,
    target = target_node,
    sender = as.character(df$sender),
    receiver = as.character(df$receiver),
    score = suppressWarnings(as.numeric(df$score)),
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[
    !is.na(edge_df$source) &
      nzchar(edge_df$source) &
      !is.na(edge_df$target) &
      nzchar(edge_df$target),
    ,
    drop = FALSE
  ]
  edge_df <- group_summary(
    df = edge_df,
    group_cols = c("source", "target", "sender", "receiver"),
    value_col = "score",
    out_col = "score",
    fun = function(x) sum(as.numeric(x), na.rm = TRUE)
  )
  edge_df <- edge_df[
    is.finite(edge_df$score) & edge_df$score > 0,
    ,
    drop = FALSE
  ]
  if (nrow(edge_df) == 0L) {
    log_message(
      "No CCC records remain for gene-level chord plotting",
      message_type = "error"
    )
  }

  strength_df <- rbind(
    data.frame(node = edge_df$source, cell = edge_df$sender, strength = edge_df$score, stringsAsFactors = FALSE),
    data.frame(node = edge_df$target, cell = edge_df$receiver, strength = edge_df$score, stringsAsFactors = FALSE)
  )
  strength_df <- stats::aggregate(
    strength ~ node + cell,
    data = strength_df,
    FUN = sum,
    na.rm = TRUE
  )
  strength_df <- strength_df[
    order(strength_df$strength, decreasing = TRUE, na.last = TRUE),
    ,
    drop = FALSE
  ]
  node_order <- unique(strength_df$node)
  cell_levels <- unique(as.character(strength_df$cell))
  cell_cols <- palette_colors(
    cell_levels,
    palette = cell_palette,
    palcolor = cell_palcolor
  )
  names(cell_cols) <- cell_levels
  node_cell <- strength_df$cell[match(node_order, strength_df$node)]
  grid.col <- unname(cell_cols[node_cell])
  names(grid.col) <- node_order
  edge_col <- unname(cell_cols[edge_df$sender])
  edge_col[is.na(edge_col)] <- "#7CAAB0"

  chord_df <- edge_df[, c("source", "target", "score"), drop = FALSE]
  colnames(chord_df) <- c("source", "target", "prob")

  old_par <- graphics::par(no.readonly = TRUE)
  circlize::circos.clear()
  on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
  on.exit(try(graphics::par(old_par), silent = TRUE), add = TRUE)

  graphics::par(mar = c(1, 1, if (!is.null(title) || !is.null(subtitle)) 3 else 1, 1))
  circlize::circos.par(
    start.degree = 90,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01),
    points.overflow.warning = FALSE
  )
  circlize::chordDiagram(
    x = chord_df,
    order = node_order,
    grid.col = grid.col,
    col = grDevices::adjustcolor(edge_col, alpha.f = link_alpha),
    transparency = max(0, min(1, 1 - link_alpha)),
    annotationTrack = "grid",
    annotationTrackHeight = c(0.03),
    preAllocateTracks = list(track.height = 0.11),
    small.gap = small.gap,
    big.gap = big.gap,
    link.sort = TRUE,
    link.decreasing = FALSE,
    link.largest.ontop = TRUE,
    directional = 1,
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow",
    link.border = NA,
    reduce = 1e-5
  )
  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim_i <- circlize::get.cell.meta.data("xlim")
      ylim_i <- circlize::get.cell.meta.data("ylim")
      sector_i <- circlize::get.cell.meta.data("sector.index")
      circlize::circos.text(
        mean(xlim_i),
        ylim_i[1],
        sector_i,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = lab.cex,
        col = grid.col[sector_i]
      )
    },
    bg.border = NA
  )

  if (!is.null(title) || !is.null(subtitle)) {
    graphics::title(
      main = title,
      sub = subtitle,
      cex.main = font.size / 10 * 1.15,
      cex.sub = font.size / 10
    )
  }
  grDevices::recordPlot()
}

ccc_sankey_plot <- function(
  pair_df,
  interaction_df = NULL,
  display_by = "aggregation",
  top_n = 20,
  edge_value = "sum",
  min_receiver_flow = 0,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  check_r("thisplot", verbose = FALSE)
  plot_theme_use <- if (identical(theme_use, "theme_scop")) "theme_this" else theme_use
  min_receiver_flow <- suppressWarnings(as.numeric(min_receiver_flow)[1])
  if (!is.finite(min_receiver_flow) || min_receiver_flow < 0) {
    min_receiver_flow <- 0
  }

  if (identical(display_by, "interaction")) {
    plot_df <- top_interactions(
      interaction_df,
      top_n = top_n,
      value_col = "score"
    )
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No interaction-level CCC records are available for sankey plotting",
        message_type = "error"
      )
    }
    if (min_receiver_flow > 0) {
      receiver_weight <- tapply(
        suppressWarnings(as.numeric(plot_df$score)),
        as.character(plot_df$receiver),
        sum,
        na.rm = TRUE
      )
      keep_receivers <- names(receiver_weight)[receiver_weight >= min_receiver_flow]
      plot_df <- plot_df[
        as.character(plot_df$receiver) %in% keep_receivers,
        ,
        drop = FALSE
      ]
      if (nrow(plot_df) == 0L) {
        log_message(
          "No interaction-level communication flows remain after applying {.arg min_receiver_flow}",
          message_type = "error"
        )
      }
    }
    meta_data <- plot_df[, c("sender", "interaction_label", "receiver"), drop = FALSE]
    colnames(meta_data) <- c("Sender", "Interaction", "Receiver")
    stat_by <- c("Sender", "Interaction", "Receiver")
    ylab_use <- "Interaction count"
  } else {
    plot_df <- top_pairs(pair_df, top_n = top_n, value_col = edge_value)
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No aggregated CCC records are available for sankey plotting",
        message_type = "error"
      )
    }
    receiver_value_col <- edge_value
    if (!receiver_value_col %in% colnames(plot_df)) {
      receiver_value_col <- if ("count" %in% colnames(plot_df)) "count" else "sum"
    }
    if (min_receiver_flow > 0 && receiver_value_col %in% colnames(plot_df)) {
      receiver_weight <- tapply(
        suppressWarnings(as.numeric(plot_df[[receiver_value_col]])),
        as.character(plot_df$receiver),
        sum,
        na.rm = TRUE
      )
      keep_receivers <- names(receiver_weight)[receiver_weight >= min_receiver_flow]
      plot_df <- plot_df[
        as.character(plot_df$receiver) %in% keep_receivers,
        ,
        drop = FALSE
      ]
      if (nrow(plot_df) == 0L) {
        log_message(
          "No communication flows remain after applying {.arg min_receiver_flow}",
          message_type = "error"
        )
      }
    }
    if (!"count" %in% colnames(plot_df)) {
      plot_df$count <- 1
    }
    if (!identical(edge_value, "count")) {
      log_message(
        paste0(
          "{.fn thisplot::StatPlot} sankey is count-based. ",
          "For {.fn CCCStatPlot} with {.code plot_type = 'sankey'}, ",
          "{.arg edge_value} is used to rank/filter pairs, but flow width is shown by interaction count."
        ),
        message_type = "warning"
      )
    }
    repeat_n <- pmax(1L, round(as.numeric(plot_df$count)))
    meta_data <- plot_df[rep(seq_len(nrow(plot_df)), repeat_n), c("sender", "receiver"), drop = FALSE]
    rownames(meta_data) <- NULL
    colnames(meta_data) <- c("Sender", "Receiver")
    stat_by <- c("Sender", "Receiver")
    ylab_use <- "Interaction count"
  }

  StatPlot(
    meta.data = meta_data,
    stat.by = stat_by,
    plot_type = "sankey",
    stat_type = "count",
    palette = cell_palette,
    palcolor = cell_palcolor,
    title = title,
    subtitle = subtitle,
    ylab = ylab_use,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = plot_theme_use,
    theme_args = theme_args,
    combine = TRUE
  )
}

ccc_distribution_plot <- function(
  interaction_df,
  plot_type = c("box", "violin"),
  top_n = 20,
  interaction.use = NULL,
  pairLR.use = NULL,
  facet_by = NULL,
  x_text_angle = 90,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  plot_type <- match.arg(plot_type)
  top_n_use <- if (!is.null(interaction.use) || !is.null(pairLR.use)) {
    NULL
  } else {
    top_n
  }
  plot_df <- top_interactions(
    interaction_df,
    top_n = top_n_use,
    value_col = "score"
  )
  if (is.null(plot_df) || nrow(plot_df) == 0L) {
    log_message(
      "No interaction-level CCC records are available for distribution plotting",
      message_type = "error"
    )
  }
  facet_by <- facet_by %||% "sender"
  if (!facet_by %in% c("sender", "receiver", "pair")) {
    facet_by <- "sender"
  }
  group_var <- if (identical(facet_by, "sender")) {
    "receiver"
  } else if (identical(facet_by, "receiver")) {
    "sender"
  } else {
    "pair"
  }
  fill_levels <- unique(plot_df[[group_var]])
  fill_cols <- palette_colors(
    fill_levels,
    palette = cell_palette,
    palcolor = cell_palcolor
  )

  if (identical(plot_type, "violin")) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data[[group_var]],
        y = score,
        fill = .data[[group_var]]
      )
    ) +
      ggplot2::geom_violin(scale = "width", trim = FALSE, alpha = 0.8) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.5) +
      ggplot2::labs(x = group_var, y = "score", fill = group_var)
    if (!identical(facet_by, "pair")) {
      p <- p +
        ggplot2::facet_wrap(
          stats::as.formula(paste("~", facet_by)),
          scales = "free_x"
        )
    }
  } else {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data[[group_var]],
        y = score,
        fill = .data[[group_var]]
      )
    ) +
      ggplot2::geom_boxplot(alpha = 0.85, outlier.size = 0.3) +
      ggplot2::labs(x = group_var, y = "score", fill = group_var)
    if (!identical(facet_by, "pair")) {
      p <- p +
        ggplot2::facet_wrap(
          stats::as.formula(paste("~", facet_by)),
          scales = "free_x"
        )
    }
  }
  p <- p +
    ggplot2::scale_fill_manual(values = fill_cols[fill_levels], drop = FALSE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_text_angle, hjust = 1)
    )
  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}

ccc_dim_network_plot <- function(
  srt,
  pair_df,
  method,
  group.by = NULL,
  reduction = NULL,
  dims = c(1, 2),
  edge_value = "sum",
  edge_threshold = 0,
  edge_size = c(0.2, 1),
  edge_color = NULL,
  edge_alpha = 0.6,
  edge_line = "curved",
  edge_curvature = 0.2,
  directed = FALSE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 4,
  node_alpha = 0.9,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  dots <- list(...)
  label_top <- isTRUE(dots[["label"]])
  label_insitu <- isTRUE(dots[["label_insitu"]])
  label_repel <- isTRUE(dots[["label_repel"]])
  label_size <- dots[["label.size"]] %||% 4
  label_fg <- dots[["label.fg"]] %||% "white"
  label_bg <- dots[["label.bg"]] %||% "black"
  label_bg_r <- dots[["label.bg.r"]] %||% 0.1
  label_repulsion <- dots[["label_repulsion"]] %||% 20
  label_point_size <- dots[["label_point_size"]] %||% 1
  label_point_color <- dots[["label_point_color"]] %||% "black"
  label_segment_color <- dots[["label_segment_color"]] %||% "black"
  lineages <- dots[["lineages"]] %||% NULL
  lineages_trim <- dots[["lineages_trim"]] %||% c(0.01, 0.99)
  lineages_span <- dots[["lineages_span"]] %||% 0.75
  lineages_palette <- dots[["lineages_palette"]] %||% "Dark2"
  lineages_palcolor <- dots[["lineages_palcolor"]] %||% NULL
  lineages_arrow <- dots[["lineages_arrow"]] %||%
    grid::arrow(length = grid::unit(0.1, "inches"))
  lineages_linewidth <- dots[["lineages_linewidth"]] %||% 1
  lineages_line_bg <- dots[["lineages_line_bg"]] %||% "white"
  lineages_line_bg_stroke <- dots[["lineages_line_bg_stroke"]] %||% 0.5
  lineages_whiskers <- dots[["lineages_whiskers"]] %||% FALSE
  lineages_whiskers_linewidth <- dots[["lineages_whiskers_linewidth"]] %||% 0.5
  lineages_whiskers_alpha <- dots[["lineages_whiskers_alpha"]] %||% 0.5
  group.by <- group.by %||% get_group_by(srt = srt, method = method)
  if (is.null(group.by) || length(group.by) != 1L || !nzchar(group.by)) {
    log_message(
      "{.arg group.by} must be a single metadata column for {.val plot_type = 'embedding_network'}",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.val {group.by}} is not in the meta.data of srt object",
      message_type = "error"
    )
  }
  if ("split.by" %in% names(dots) && !is.null(dots[["split.by"]])) {
    log_message(
      "{.val plot_type = 'embedding_network'} does not support {.arg split.by}",
      message_type = "error"
    )
  }
  if ("combine" %in% names(dots) && !isTRUE(dots[["combine"]])) {
    log_message(
      "{.val plot_type = 'embedding_network'} requires {.arg combine = TRUE}",
      message_type = "error"
    )
  }
  reduction_use <- if (is.null(reduction)) {
    DefaultReduction(srt)
  } else {
    DefaultReduction(srt, pattern = reduction)
  }
  protected_args <- c(
    "srt",
    "group.by",
    "reduction",
    "dims",
    "palette",
    "palcolor",
    "cell_palette",
    "cell_palcolor",
    "link_palette",
    "link_palcolor",
    "title",
    "subtitle",
    "legend.position",
    "legend.direction",
    "legend.title",
    "theme_use",
    "theme_args",
    "combine",
    "lineages",
    "lineages_trim",
    "lineages_span",
    "lineages_palette",
    "lineages_palcolor",
    "lineages_arrow",
    "lineages_linewidth",
    "lineages_line_bg",
    "lineages_line_bg_stroke",
    "lineages_whiskers",
    "lineages_whiskers_linewidth",
    "lineages_whiskers_alpha",
    "label",
    "label_insitu",
    "label_repel",
    "label.size",
    "label.fg",
    "label.bg",
    "label.bg.r",
    "label_repulsion",
    "label_point_size",
    "label_point_color",
    "label_segment_color"
  )
  dots <- dots[!names(dots) %in% protected_args]

  base_plot <- do.call(
    CellDimPlot,
    c(
      list(
        srt = srt,
        group.by = group.by,
        reduction = reduction,
        dims = dims,
        palette = cell_palette,
        palcolor = cell_palcolor,
        title = title,
        subtitle = subtitle,
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.title = legend.title,
        theme_use = theme_use,
        theme_args = theme_args,
        combine = TRUE
      ),
      dots
    )
  )

  plot_data <- ccc_dim_network_plot_data(
    srt = srt,
    group.by = group.by,
    reduction = reduction_use,
    dims = dims,
    cells = dots[["cells"]] %||% NULL,
    show_na = dots[["show_na"]] %||% FALSE
  )

  overlay <- ccc_dim_network_layers(
    plot_data = plot_data,
    pair_df = pair_df,
    levels = levels(srt@meta.data[[group.by]]),
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    link_palette = link_palette,
    link_palcolor = link_palcolor,
    edge_value = edge_value,
    edge_threshold = edge_threshold,
    edge_size = edge_size,
    edge_color = edge_color,
    edge_alpha = edge_alpha,
    edge_line = edge_line,
    edge_curvature = edge_curvature,
    directed = directed,
    arrow_type = arrow_type,
    arrow_angle = arrow_angle,
    arrow_length = arrow_length,
    node_size = node_size,
    node_alpha = node_alpha
  )
  if (is.null(overlay)) {
    log_message(
      "No CCC edges are available for {.val plot_type = 'embedding_network'}",
      message_type = "error"
    )
  }
  lineage_layer <- NULL
  if (!is.null(lineages)) {
    lineage_layer <- ccc_dim_network_lineage_layer(
      srt = srt,
      lineages = lineages,
      reduction = reduction_use,
      dims = dims,
      cells = dots[["cells"]] %||% NULL,
      trim = lineages_trim,
      span = lineages_span,
      palette = lineages_palette,
      palcolor = lineages_palcolor,
      lineages_arrow = lineages_arrow,
      linewidth = lineages_linewidth,
      line_bg = lineages_line_bg,
      line_bg_stroke = lineages_line_bg_stroke,
      whiskers = lineages_whiskers,
      whiskers_linewidth = lineages_whiskers_linewidth,
      whiskers_alpha = lineages_whiskers_alpha
    )
  }
  label_layer <- NULL
  if (isTRUE(label_top)) {
    label_layer <- ccc_dim_network_label_layer(
      plot_data = plot_data,
      label_insitu = label_insitu,
      label_repel = label_repel,
      label_size = label_size,
      label_fg = label_fg,
      label_bg = label_bg,
      label_bg_r = label_bg_r,
      label_repulsion = label_repulsion,
      label_point_size = label_point_size,
      label_point_color = label_point_color,
      label_segment_color = label_segment_color
    )
  }

  suppressWarnings(base_plot + lineage_layer + overlay + label_layer)
}

ccc_dim_network_plot_data <- function(
  srt,
  group.by,
  reduction,
  dims = c(1, 2),
  cells = NULL,
  show_na = FALSE
) {
  emb <- Seurat::Embeddings(srt, reduction = reduction)
  if (max(dims) > ncol(emb)) {
    log_message(
      "{.arg dims} exceeds the available dimensions in reduction {.val {reduction}}",
      message_type = "error"
    )
  }

  cells_use <- rownames(emb)
  if (!is.null(cells)) {
    cells_use <- intersect(cells_use, cells)
  }
  emb <- emb[cells_use, dims, drop = FALSE]
  meta <- srt@meta.data[cells_use, group.by, drop = FALSE]
  colnames(meta) <- "group.by"

  plot_data <- data.frame(
    x = emb[, 1],
    y = emb[, 2],
    group.by = meta[, "group.by"],
    stringsAsFactors = FALSE
  )

  if (isTRUE(show_na) && any(is.na(plot_data$group.by))) {
    plot_data$group.by <- as.character(plot_data$group.by)
    plot_data$group.by[is.na(plot_data$group.by)] <- "NA"
  }
  plot_data
}

ccc_dim_network_label_layer <- function(
  plot_data,
  label_insitu = FALSE,
  label_repel = FALSE,
  label_size = 4,
  label_fg = "white",
  label_bg = "black",
  label_bg_r = 0.1,
  label_repulsion = 20,
  label_point_size = 1,
  label_point_color = "black",
  label_segment_color = "black"
) {
  if (
    is.null(plot_data) ||
      nrow(plot_data) == 0L ||
      !"group.by" %in% colnames(plot_data)
  ) {
    return(NULL)
  }

  label_df <- stats::aggregate(
    plot_data[, c("x", "y"), drop = FALSE],
    by = list(label = plot_data[["group.by"]]),
    FUN = stats::median
  )
  label_df <- label_df[!is.na(label_df$label), , drop = FALSE]
  if (nrow(label_df) == 0L) {
    return(NULL)
  }
  if (isFALSE(label_insitu)) {
    label_df$label <- as.character(label_df$label)
  }

  if (isTRUE(label_repel)) {
    list(
      ggplot2::geom_point(
        data = label_df,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        color = label_point_color,
        size = label_point_size,
        inherit.aes = FALSE,
        show.legend = FALSE
      ),
      ggrepel::geom_text_repel(
        data = label_df,
        mapping = ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          label = .data[["label"]]
        ),
        fontface = "bold",
        min.segment.length = 0,
        segment.color = label_segment_color,
        point.size = label_point_size,
        max.overlaps = 100,
        force = label_repulsion,
        color = label_fg,
        bg.color = label_bg,
        bg.r = label_bg_r,
        size = label_size,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    )
  } else {
    list(
      ggrepel::geom_text_repel(
        data = label_df,
        mapping = ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          label = .data[["label"]]
        ),
        fontface = "bold",
        min.segment.length = 0,
        segment.color = label_segment_color,
        point.size = NA,
        max.overlaps = 100,
        force = 0,
        color = label_fg,
        bg.color = label_bg,
        bg.r = label_bg_r,
        size = label_size,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    )
  }
}

ccc_dim_network_lineage_layer <- function(
  srt,
  lineages,
  reduction,
  dims = c(1, 2),
  cells = NULL,
  trim = c(0.01, 0.99),
  span = 0.75,
  palette = "Dark2",
  palcolor = NULL,
  lineages_arrow = grid::arrow(length = grid::unit(0.1, "inches")),
  linewidth = 1,
  line_bg = "white",
  line_bg_stroke = 0.5,
  whiskers = FALSE,
  whiskers_linewidth = 0.5,
  whiskers_alpha = 0.5
) {
  lineages_layers <- LineagePlot(
    srt = srt,
    lineages = lineages,
    reduction = reduction,
    dims = dims,
    cells = cells,
    trim = trim,
    span = span,
    palette = palette,
    palcolor = palcolor,
    lineages_arrow = lineages_arrow,
    linewidth = linewidth,
    line_bg = line_bg,
    line_bg_stroke = line_bg_stroke,
    whiskers = whiskers,
    whiskers_linewidth = whiskers_linewidth,
    whiskers_alpha = whiskers_alpha,
    return_layer = TRUE
  )
  c(list(ggnewscale::new_scale_color()), lineages_layers$curve_layer)
}

ccc_dim_network_layers <- function(
  plot_data,
  pair_df,
  levels = NULL,
  cell_palette = "RdBu",
  cell_palcolor = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  edge_value = "sum",
  edge_threshold = 0,
  edge_size = c(0.2, 1),
  edge_color = NULL,
  edge_alpha = 0.6,
  edge_line = "curved",
  edge_curvature = 0.2,
  directed = FALSE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 4,
  node_alpha = 0.9
) {
  ccc_group_key <- function(x) {
    x <- as.character(x)
    x <- tolower(trimws(x))
    gsub("[^[:alnum:]]+", "", x)
  }
  ccc_pick_edge_metric <- function(df, requested) {
    candidates <- unique(c(requested, "sum", "mean", "max", "count", "score"))
    candidates <- candidates[candidates %in% colnames(df)]
    if (length(candidates) == 0L) {
      return(NULL)
    }
    for (nm in candidates) {
      vals <- suppressWarnings(as.numeric(df[[nm]]))
      if (any(is.finite(vals))) {
        return(nm)
      }
    }
    NULL
  }

  edge_color_is_missing <- is.null(edge_color) ||
    length(edge_color) == 0L ||
    (length(edge_color) == 1L && is.na(edge_color))
  if (
    is.null(plot_data) ||
      nrow(plot_data) == 0L ||
      is.null(pair_df) ||
      nrow(pair_df) == 0L
  ) {
    return(NULL)
  }

  node_df <- stats::aggregate(
    plot_data[, c("x", "y"), drop = FALSE],
    by = list(group = plot_data[["group.by"]]),
    FUN = stats::median
  )
  node_df <- node_df[!is.na(node_df$group), , drop = FALSE]
  if (nrow(node_df) == 0L) {
    return(NULL)
  }
  node_df$group_key <- ccc_group_key(node_df$group)
  node_df <- node_df[!duplicated(node_df$group_key), , drop = FALSE]

  pair_df$sender_key <- ccc_group_key(pair_df$sender)
  pair_df$receiver_key <- ccc_group_key(pair_df$receiver)
  edge_df <- pair_df[
    pair_df$sender_key %in% node_df$group_key &
      pair_df$receiver_key %in% node_df$group_key, ,
    drop = FALSE
  ]
  if (nrow(edge_df) == 0L) {
    return(NULL)
  }

  edge_metric <- ccc_pick_edge_metric(edge_df, edge_value)
  if (is.null(edge_metric)) {
    return(NULL)
  }
  edge_df$weight <- abs(suppressWarnings(as.numeric(edge_df[[edge_metric]])))
  edge_df <- edge_df[is.finite(edge_df$weight), , drop = FALSE]
  if (nrow(edge_df) == 0L) {
    return(NULL)
  }

  node_weight_df <- rbind(
    data.frame(group = as.character(edge_df$sender), weight = edge_df$weight, stringsAsFactors = FALSE),
    data.frame(group = as.character(edge_df$receiver), weight = edge_df$weight, stringsAsFactors = FALSE)
  )
  node_weight_df <- stats::aggregate(
    weight ~ group,
    data = node_weight_df,
    FUN = sum,
    na.rm = TRUE
  )
  node_df <- merge(node_df, node_weight_df, by = "group", all.x = TRUE)
  node_df$weight[!is.finite(node_df$weight)] <- 0

  if (!isTRUE(directed)) {
    edge_df <- stats::aggregate(
      weight ~ sender + receiver,
      data = transform(
        edge_df,
        sender = pmin(as.character(sender), as.character(receiver)),
        receiver = pmax(as.character(sender), as.character(receiver))
      ),
      FUN = sum,
      na.rm = TRUE
    )
  }
  edge_df$sender_key <- ccc_group_key(edge_df$sender)
  edge_df$receiver_key <- ccc_group_key(edge_df$receiver)

  edge_df <- edge_df[
    edge_df$weight >= edge_threshold, ,
    drop = FALSE
  ]
  if (nrow(edge_df) == 0L) {
    return(NULL)
  }

  node_from <- node_df[, c("group_key", "group", "x", "y"), drop = FALSE]
  colnames(node_from) <- c("sender_key", "sender_display", "x_from", "y_from")
  node_to <- node_df[, c("group_key", "group", "x", "y"), drop = FALSE]
  colnames(node_to) <- c("receiver_key", "receiver_display", "x_to", "y_to")

  edge_df <- merge(
    edge_df,
    node_from,
    by = "sender_key",
    all.x = TRUE
  )
  edge_df <- merge(
    edge_df,
    node_to,
    by = "receiver_key",
    all.x = TRUE
  )
  edge_df <- edge_df[
    !is.na(edge_df$x_from) & !is.na(edge_df$x_to), ,
    drop = FALSE
  ]
  if (nrow(edge_df) == 0L) {
    return(NULL)
  }
  self_edge_df <- edge_df[edge_df$sender == edge_df$receiver, , drop = FALSE]
  edge_df <- edge_df[edge_df$sender != edge_df$receiver, , drop = FALSE]

  group_levels <- levels %||% unique(as.character(node_df$group))
  group_levels <- unique(c(group_levels, as.character(node_df$group)))
  colors <- palette_colors(
    group_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  edge_cols <- palette_colors(
    unique(edge_df$sender_display),
    palette = link_palette,
    palcolor = link_palcolor,
    NA_keep = TRUE
  )

  x_span <- diff(range(node_df$x, na.rm = TRUE))
  y_span <- diff(range(node_df$y, na.rm = TRUE))
  loop_dx <- if (is.finite(x_span) && x_span > 0) x_span * 0.06 else 0.25
  loop_dy <- if (is.finite(y_span) && y_span > 0) y_span * 0.06 else 0.25
  if (nrow(self_edge_df) > 0L) {
    self_edge_df$x_loop_from <- self_edge_df$x_from - loop_dx
    self_edge_df$y_loop_from <- self_edge_df$y_from + loop_dy * 0.4
    self_edge_df$x_loop_to <- self_edge_df$x_from + loop_dx
    self_edge_df$y_loop_to <- self_edge_df$y_from + loop_dy * 0.4
  }

  edge_geom <- if (identical(edge_line, "straight")) {
    geoms <- list()
    if (nrow(edge_df) > 0L && isTRUE(edge_color_is_missing)) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_segment(
        data = edge_df,
        mapping = ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = weight,
          color = sender_display
        ),
        alpha = edge_alpha,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    } else if (nrow(edge_df) > 0L) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_segment(
        data = edge_df,
        mapping = ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = weight
        ),
        color = edge_color,
        alpha = edge_alpha,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
    if (nrow(self_edge_df) > 0L && isTRUE(edge_color_is_missing)) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_curve(
        data = self_edge_df,
        mapping = ggplot2::aes(
          x = x_loop_from,
          y = y_loop_from,
          xend = x_loop_to,
          yend = y_loop_to,
          linewidth = weight,
          color = sender_display
        ),
        alpha = edge_alpha,
        curvature = 1.2,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    } else if (nrow(self_edge_df) > 0L) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_curve(
        data = self_edge_df,
        mapping = ggplot2::aes(
          x = x_loop_from,
          y = y_loop_from,
          xend = x_loop_to,
          yend = y_loop_to,
          linewidth = weight
        ),
        color = edge_color,
        alpha = edge_alpha,
        curvature = 1.2,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
    geoms
  } else {
    geoms <- list()
    if (nrow(edge_df) > 0L && isTRUE(edge_color_is_missing)) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_curve(
        data = edge_df,
        mapping = ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = weight,
          color = sender_display
        ),
        alpha = edge_alpha,
        curvature = edge_curvature,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    } else if (nrow(edge_df) > 0L) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_curve(
        data = edge_df,
        mapping = ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = weight
        ),
        color = edge_color,
        alpha = edge_alpha,
        curvature = edge_curvature,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
    if (nrow(self_edge_df) > 0L && isTRUE(edge_color_is_missing)) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_curve(
        data = self_edge_df,
        mapping = ggplot2::aes(
          x = x_loop_from,
          y = y_loop_from,
          xend = x_loop_to,
          yend = y_loop_to,
          linewidth = weight,
          color = sender_display
        ),
        alpha = edge_alpha,
        curvature = 1.2,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    } else if (nrow(self_edge_df) > 0L) {
      geoms[[length(geoms) + 1L]] <- ggplot2::geom_curve(
        data = self_edge_df,
        mapping = ggplot2::aes(
          x = x_loop_from,
          y = y_loop_from,
          xend = x_loop_to,
          yend = y_loop_to,
          linewidth = weight
        ),
        color = edge_color,
        alpha = edge_alpha,
        curvature = 1.2,
        arrow = if (isTRUE(directed)) {
          grid::arrow(
            type = arrow_type,
            angle = arrow_angle,
            length = arrow_length
          )
        } else {
          NULL
        },
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
    geoms
  }

  edge_color_scale <- if (isTRUE(edge_color_is_missing)) {
    ggplot2::scale_color_manual(
      values = edge_cols,
      drop = FALSE,
      guide = "none"
    )
  } else {
    NULL
  }

  layer_list <- c(
    list(ggnewscale::new_scale_color()),
    edge_geom,
    list(edge_color_scale),
    list(
      ggnewscale::new_scale_fill(),
      ggplot2::scale_size_continuous(
        range = c(node_size, node_size * 2.6),
        guide = "none"
      ),
      ggplot2::scale_linewidth_continuous(range = edge_size, guide = "none"),
      ggplot2::geom_point(
        data = node_df,
        mapping = ggplot2::aes(x = x, y = y, fill = group, size = weight),
        shape = 21,
        alpha = node_alpha,
        color = "grey20",
        stroke = 0.8,
        show.legend = FALSE,
        inherit.aes = FALSE
      ),
      ggplot2::scale_fill_manual(
        values = colors[group_levels],
        drop = FALSE,
        guide = "none"
      )
    )
  )
  layer_list[!vapply(layer_list, is.null, logical(1))]
}

ccc_bar_plot <- function(
  df,
  top_n = 20,
  title = NULL,
  subtitle = NULL,
  link_palette = "RdBu",
  link_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No interaction records are available for bar plotting",
      message_type = "error"
    )
  }
  agg <- stats::aggregate(
    x = df$score,
    by = list(
      interaction_name = df$interaction_name,
      sender = df$sender,
      receiver = df$receiver
    ),
    FUN = sum,
    na.rm = TRUE
  )
  agg$pair <- paste(agg$sender, "->", agg$receiver)
  agg <- agg[order(agg$x, decreasing = TRUE), , drop = FALSE]
  agg <- utils::head(agg, top_n)
  cols <- palette_colors(
    unique(agg$pair),
    palette = link_palette,
    palcolor = link_palcolor
  )
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = stats::reorder(interaction_name, x), y = x, fill = pair)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = "Aggregated score", fill = "Pair")
  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}

ccc_standardize_ligand_target_df <- function(
  df,
  top_n = 20,
  sender.use = NULL,
  receiver.use = NULL,
  context_df = NULL,
  sender_default = NULL,
  receiver_default = NULL
) {
  df <- standardize_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No ligand-target table is available for plotting",
      message_type = "error"
    )
  }
  df <- rename_by_candidates(
    df,
    "ligand",
    c("ligand", "from", "test_ligand", "ligand_oi", "ligand_source")
  )
  df <- rename_by_candidates(
    df,
    "target",
    c("target", "to", "gene", "target_gene", "target_genes", "gene_oi")
  )
  df <- rename_by_candidates(
    df,
    "sender",
    c("sender", "sender_celltype", "celltype_sender", "source_celltype")
  )
  df <- rename_by_candidates(
    df,
    "receiver",
    c("receiver", "receiver_celltype", "celltype_receiver", "target_celltype")
  )
  df <- rename_by_candidates(
    df,
    "weight",
    c(
      "weight",
      "score",
      "regulatory_potential",
      "pearson",
      "aupr_corrected",
      "aupr",
      "activity",
      "prioritization_score"
    )
  )

  if (!"weight" %in% colnames(df)) {
    numeric_cols <- setdiff(
      colnames(df)[vapply(df, is.numeric, logical(1))],
      c("ligand", "target")
    )
    if (length(numeric_cols) > 0L) {
      colnames(df)[match(numeric_cols[1], colnames(df))] <- "weight"
    }
  }

  if (!all(c("ligand", "target", "weight") %in% colnames(df))) {
    log_message(
      "Ligand-target table must contain ligand, target and weight columns",
      message_type = "error"
    )
  }

  sender_default <- unique(stats::na.omit(as.character(sender_default)))
  receiver_default <- unique(stats::na.omit(as.character(receiver_default)))
  if (!"sender" %in% colnames(df) && length(sender_default) == 1L) {
    df$sender <- sender_default
  }
  if (!"receiver" %in% colnames(df) && length(receiver_default) == 1L) {
    df$receiver <- receiver_default
  }

  context_df <- standardize_df(context_df)
  if (nrow(context_df) > 0L) {
    context_df <- standardize_long_df(context_df)
  }

  sender_requested <- !is.null(sender.use) && length(sender.use) > 0L
  receiver_requested <- !is.null(receiver.use) && length(receiver.use) > 0L
  direct_sender <- "sender" %in% colnames(df)
  direct_receiver <- "receiver" %in% colnames(df)

  if (sender_requested && direct_sender) {
    df <- df[df$sender %in% sender.use, , drop = FALSE]
  }
  if (receiver_requested && direct_receiver) {
    df <- df[df$receiver %in% receiver.use, , drop = FALSE]
  }

  if (
    (sender_requested && !direct_sender) ||
      (receiver_requested && !direct_receiver)
  ) {
    if (
      nrow(context_df) > 0L &&
        all(c("ligand", "sender", "receiver") %in% colnames(context_df))
    ) {
      context_keep <- context_df
      if (sender_requested) {
        context_keep <- context_keep[
          context_keep$sender %in% sender.use, ,
          drop = FALSE
        ]
      }
      if (receiver_requested) {
        context_keep <- context_keep[
          context_keep$receiver %in% receiver.use, ,
          drop = FALSE
        ]
      }
      keep_ligands <- unique(as.character(context_keep$ligand))
      keep_ligands <- keep_ligands[!is.na(keep_ligands) & nzchar(keep_ligands)]
      if (length(keep_ligands) > 0L) {
        df <- df[df$ligand %in% keep_ligands, , drop = FALSE]
      } else {
        df <- df[0, , drop = FALSE]
      }
    } else {
      log_message(
        paste0(
          "{.arg sender.use}/{.arg receiver.use} were provided, but the ",
          "ligand-target table does not contain sender/receiver columns and ",
          "no compatible long-table context is available for filtering."
        ),
        message_type = "warning"
      )
    }
  }

  if (nrow(df) == 0L) {
    log_message(
      "No ligand-target records remain after applying the current filters",
      message_type = "error"
    )
  }

  df <- df[order(df$weight, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(
    df,
    max(top_n, 1L) * max(3L, min(10L, length(unique(df$ligand))))
  )
  ligand_levels <- unique(df$ligand)
  target_levels <- unique(df$target)
  if (length(ligand_levels) > top_n) {
    ligand_levels <- ligand_levels[seq_len(top_n)]
    df <- df[df$ligand %in% ligand_levels, , drop = FALSE]
    target_levels <- unique(df$target)
  }
  df$weight <- suppressWarnings(as.numeric(df$weight))
  df$ligand <- factor(df$ligand, levels = rev(ligand_levels))
  df$target <- factor(df$target, levels = target_levels)
  df
}

ccc_ligand_target_plot <- function(
  df,
  top_n = 20,
  title = NULL,
  subtitle = NULL,
  value_palette = "RdBu",
  value_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  df <- ccc_standardize_ligand_target_df(df = df, top_n = top_n)

  cols <- palette_colors(
    palette = value_palette,
    palcolor = value_palcolor,
    n = 9
  )
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = target, y = ligand, fill = weight)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colours = cols) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "Target", y = "Ligand", fill = "weight")
  finalize_cc_plot(
    p,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}
