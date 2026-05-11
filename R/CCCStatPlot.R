#' @title CCC statistical distribution and summary plots
#'
#' @md
#' @param srt A `Seurat` object.
#' @param method Communication result type to use.
#' @param plot_type Plot type. One of:
#' - `"bar"` — horizontal bar chart of top pairs or interactions according to
#'   `display_by`.
#' - `"sankey"` — alluvial/sankey flow diagram.
#' - `"box"` / `"violin"` — distribution of interaction scores across sender-receiver pairs.
#' - `"comparison"` — comparison bars at overall or celltype level.
#' - `"lr_contribution"` — ligand-receptor contribution bar plot.
#' - `"gene"` —  pathway-related ligand/receptor gene expression panel.
#' - `"ranknet"` — pathway ranking comparison plot.
#' - `"scatter"` — outgoing vs. incoming signaling strength scatter.
#' - `"role_change"` — signaling change scatter for one cell identity.
#' @param condition Result name or comparison name.
#' @param dataset Dataset index or name.
#' @param comparison Comparison indices or names.
#' @param display_by Whether to summarize by `"aggregation"` or `"interaction"`.
#' @param sender.use Sender cell types to keep.
#' @param receiver.use Receiver cell types to keep.
#' @param ligand.use Ligands to keep.
#' @param receptor.use Receptors to keep.
#' @param interaction.use Interaction names to keep.
#' @param signaling Signaling pathway to focus on.
#' @param pairLR.use Specific ligand-receptor pair(s) to keep.
#' @param slot.name CellChat slot name.
#' @param thresh Significance threshold used when extracting communication results.
#' @param measure Summary measure for CellChat objects.
#' @param pattern Pattern used for pathway role plots.
#' @param compare_by Comparison mode for CellChat summary plots.
#' @param value Value column or summary statistic to use.
#' @param top_n Number of top records to retain.
#' @param x_text_angle Rotation angle for x-axis labels.
#' @param facet_by Faceting variable for interaction-level plots.
#' @param min_receiver_flow For `"sankey"`: minimum total receiver-side flow
#'   retained after top-N ranking. Useful when many small receiver nodes make
#'   the right side unreadable.
#' @param edge_value Aggregation statistic for network edges.
#' @param edge_threshold Minimum edge value to keep.
#' @param link_alpha Alpha used for network edges.
#' @param palette Main palette name.
#' @param palcolor Main custom palette colors.
#' @param cell_palette Cell annotation palette name.
#' @param cell_palcolor Custom cell annotation colors.
#' @param link_palette Link palette name.
#' @param link_palcolor Custom link palette colors.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param legend.position Legend position.
#' @param legend.direction Legend direction.
#' @param font.size Base font size.
#' @param theme_use Theme function used for styling.
#' @param theme_args Arguments passed to the theme function.
#' @param grid_major Whether to show major panel grid lines for applicable statistical panels.
#' Default is `TRUE`.
#' @param grid_major_colour Color of major panel grid lines.
#' @param grid_major_linetype Linetype of major panel grid lines.
#' @param grid_major_linewidth Line width of major panel grid lines.
#' @param verbose Whether to print messages.
#' @param combine Whether to combine multiple panels.
#' @param nrow Number of rows in combined layout.
#' @param ncol Number of columns in combined layout.
#' @param ... Additional plot-specific options.
#' @param stat_type For `"bar"`: what to summarize per interaction. One of
#'   `"score"` (total aggregated score) or `"count"` (number of significant
#'   interactions).
#'
#' @return A ggplot or recorded base plot object.
#' @export
#'
#' @examples
#' if (requireNamespace("CellChat", quietly = TRUE)) {
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#'
#' pc1 <- Seurat::Embeddings(pancreas_sub, "Standardpca")[, 1]
#' ct <- as.character(pancreas_sub$CellType)
#' ct_medians <- tapply(pc1, ct, median)
#' pancreas_sub$Condition <- ifelse(
#'   pc1 > ct_medians[ct],
#'   "ConditionA",
#'   "ConditionB"
#' )
#'
#' pancreas_sub <- RunCellChat(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group_column = "Condition",
#'   group_cmp = list(c("ConditionA", "ConditionB")),
#'   species = "Mus_musculus"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "sankey",
#'   display_by = "aggregation",
#'   top_n = 20
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "sankey",
#'   display_by = "interaction",
#'   top_n = 20
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "box",
#'   facet_by = "sender",
#'   top_n = 200
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "violin",
#'   facet_by = "receiver",
#'   top_n = 200
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "bar",
#'   palette = "Paired",
#'   top_n = 100
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "scatter"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "lr_contribution",
#'   signaling = "MK"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "gene",
#'   signaling = "MK"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "comparison",
#'   measure = "count",
#'   compare_by = "overall"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "comparison",
#'   measure = "weight",
#'   compare_by = "celltype",
#'   pattern = "all"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "ranknet"
#' )
#'
#' CCCStatPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   idents.use = "Ductal",
#'   plot_type = "role_change"
#' )
#' }
CCCStatPlot <- function(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c(
    "bar",
    "sankey",
    "box",
    "violin",
    "role_scatter",
    "role_network",
    "role_network_marsilea",
    "pathway_summary",
    "comparison",
    "lr_contribution",
    "gene",
    "ranknet",
    "scatter",
    "role_change"
  ),
  display_by = c("aggregation", "interaction"),
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  measure = c("count", "weight"),
  pattern = c("all", "outgoing", "incoming"),
  compare_by = c("overall", "celltype"),
  value = "score",
  stat_type = c("score", "count"),
  top_n = 20,
  x_text_angle = 90,
  min_receiver_flow = 0,
  link_alpha = 0.6,
  facet_by = NULL,
  edge_value = c("sum", "mean", "max", "count"),
  edge_threshold = 0,
  palette = "Chinese",
  palcolor = NULL,
  cell_palette = NULL,
  cell_palcolor = NULL,
  link_palette = NULL,
  link_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  plot_type <- match.arg(plot_type)
  display_by <- match.arg(display_by)
  measure <- match.arg(measure)
  pattern <- match.arg(pattern)
  compare_by <- match.arg(compare_by)
  edge_value <- match.arg(edge_value)
  stat_type <- match.arg(stat_type)

  palette_cfg <- ccc_palettes(
    palette = palette,
    palcolor = palcolor,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    link_palette = link_palette,
    link_palcolor = link_palcolor
  )

  method <- detect_method(srt = srt, method = method)
  dots <- list(...)
  finish_plot <- function(plot) {
    plot
  }
  return.data <- isTRUE(dots[["return.data"]])
  idents.use <- dots[["idents.use"]] %||% NULL
  idents.use <- ccc_alias_arg(dots, "idents_use", idents.use)
  sender.use <- ccc_alias_arg(dots, "sender_use", sender.use)
  receiver.use <- ccc_alias_arg(dots, "receiver_use", receiver.use)
  ligand.use <- ccc_alias_arg(dots, "ligand_use", ligand.use)
  receptor.use <- ccc_alias_arg(dots, "receptor_use", receptor.use)
  interaction.use <- ccc_alias_arg(dots, "interaction_use", interaction.use)
  pairLR.use <- ccc_alias_arg(dots, "pair_lr_use", pairLR.use)
  thresh <- ccc_alias_arg(dots, "pvalue_threshold", thresh)
  signaling.exclude <- dots[["signaling.exclude"]] %||% NULL
  aspect.ratio <- dots[["aspect.ratio"]] %||% NULL
  sample_col <- dots[["sample_col"]] %||% NULL
  sample_col <- ccc_alias_arg(
    dots,
    c("context_col", "condition_col", "dataset_col"),
    sample_col
  )

  if (identical(plot_type, "role_scatter")) {
    plot_type <- "scatter"
  }
  if (plot_type %in% c("role_network", "role_network_marsilea")) {
    return(finish_plot(CCCHeatmap(
      srt = srt,
      method = method,
      condition = condition,
      dataset = dataset,
      comparison = comparison,
      plot_type = "role_heatmap",
      signaling = signaling,
      pattern = pattern,
      top_n = top_n,
      title = title,
      subtitle = subtitle,
      palette = palette,
      palcolor = palcolor,
      cell_palette = cell_palette,
      cell_palcolor = cell_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args,
      verbose = verbose
    )))
  }

  if (plot_type %in% c("ranknet", "role_change")) {
    if (!identical(method, "CellChat")) {
      log_message(
        paste0(
          "{.arg plot_type} = {.val {plot_type}} is currently only supported ",
          "for {.pkg CellChat} results"
        ),
        message_type = "error"
      )
    }
    if (identical(plot_type, "ranknet")) {
      return(finish_plot(ccc_cellchat_ranknet_plot(
        srt = srt,
        condition = condition,
        comparison = comparison,
        signaling = signaling,
        sender.use = sender.use,
        receiver.use = receiver.use,
        thresh = thresh,
        return.data = return.data,
        title = title,
        subtitle = subtitle,
        palette = palette_cfg$cell_palette,
        palcolor = palette_cfg$cell_palcolor,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme_use = theme_use,
        theme_args = theme_args,
        font.size = font.size
      )))
    }

    if (is.null(idents.use)) {
      log_message(
        "{.arg idents.use} must be provided for {.val plot_type = 'role_change'}",
        message_type = "error"
      )
    }
    return(finish_plot(ccc_cellchat_role_change_plot(
      srt = srt,
      condition = condition,
      comparison = comparison,
      idents.use = idents.use,
      signaling = signaling,
      signaling.exclude = signaling.exclude,
      return.data = return.data,
      title = title,
      subtitle = subtitle,
      palette = palette_cfg$cell_palette,
      palcolor = palette_cfg$cell_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio,
      font.size = font.size
    )))
  }

  if (identical(plot_type, "lr_contribution")) {
    if (is.null(signaling)) {
      df <- ccc_long_table_for_method(
        srt = srt,
        method = method,
        condition = condition,
        dataset = dataset,
        slot.name = slot.name,
        signaling = signaling,
        pairLR.use = pairLR.use,
        sources.use = sender.use,
        targets.use = receiver.use,
        thresh = thresh
      )
      df <- standardize_long_df(df)
      df <- filter_long_df(
        df = df,
        sender.use = sender.use,
        receiver.use = receiver.use,
        ligand.use = ligand.use,
        receptor.use = receptor.use,
        interaction.use = interaction.use,
        signaling = signaling,
        pairLR.use = pairLR.use
      )
      df <- ccc_assign_plot_score(df = df, value = value)
      df <- ccc_mark_significance(df, thresh = thresh)
      df <- prepare_plot_df(df)
      return(finish_plot(ccc_stat_bar_plot(
        df = df,
        pair_df = pair_plot_df(df),
        stat_type = if (identical(value, "count")) "count" else "score",
        display_by = "interaction",
        value = value,
        top_n = top_n,
        title = title,
        subtitle = subtitle,
        link_palette = palette_cfg$link_palette,
        link_palcolor = palette_cfg$link_palcolor,
        legend.position = legend.position,
        legend.direction = legend.direction,
        font.size = font.size,
        theme_use = theme_use,
        theme_args = theme_args
      )))
    }
    if (!identical(method, "CellChat")) {
      return(finish_plot(ccc_generic_lr_contribution_plot(
        srt = srt,
        method = method,
        signaling = signaling,
        sender.use = sender.use,
        receiver.use = receiver.use,
        ligand.use = ligand.use,
        receptor.use = receptor.use,
        interaction.use = interaction.use,
        pairLR.use = pairLR.use,
        top_n = top_n,
        value = value,
        return.data = return.data,
        title = title,
        subtitle = subtitle,
        link_palette = palette_cfg$link_palette,
        link_palcolor = palette_cfg$link_palcolor,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme_use = theme_use,
        theme_args = theme_args,
        font.size = font.size,
        combine = combine,
        nrow = nrow,
        ncol = ncol
      )))
    }
    return(finish_plot(ccc_lr_contribution_plot(
      srt = srt,
      condition = condition,
      dataset = dataset,
      signaling = signaling,
      sender.use = sender.use,
      receiver.use = receiver.use,
      thresh = thresh,
      top_n = top_n,
      return.data = return.data,
      title = title,
      subtitle = subtitle,
      link_palette = palette_cfg$link_palette,
      link_palcolor = palette_cfg$link_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size,
      combine = combine,
      nrow = nrow,
      ncol = ncol
    )))
  }

  if (identical(plot_type, "gene")) {
    return(finish_plot(
      ccc_feature_plot(
        srt = srt,
        method = method,
        condition = condition,
        dataset = dataset,
        signaling = signaling,
        sender.use = sender.use,
        receiver.use = receiver.use,
        ligand.use = ligand.use,
        receptor.use = receptor.use,
        interaction.use = interaction.use,
        pairLR.use = pairLR.use,
        slot.name = slot.name,
        thresh = thresh,
        top_n = top_n,
        return.data = return.data,
        title = title,
        subtitle = subtitle,
        palette = palette_cfg$cell_palette,
        palcolor = palette_cfg$cell_palcolor,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme_use = theme_use,
        theme_args = theme_args,
        aspect.ratio = aspect.ratio,
        combine = combine,
        nrow = nrow,
        ncol = ncol,
        dots = dots
      )
    ))
  }

  df <- ccc_long_table_for_method(
    srt = srt,
    method = method,
    condition = condition,
    dataset = dataset,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sender.use,
    targets.use = receiver.use,
    thresh = thresh
  )

  df <- standardize_long_df(df)
  df <- filter_long_df(
    df = df,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    signaling = signaling,
    pairLR.use = pairLR.use
  )

  df <- ccc_assign_plot_score(df = df, value = value)
  df <- ccc_mark_significance(df, thresh = thresh)
  df <- prepare_plot_df(df)
  pair_df <- pair_plot_df(df)
  interaction_df <- interaction_plot_df(df)

  if (identical(plot_type, "pathway_summary")) {
    return(finish_plot(ccc_pathway_summary_plot(
      df = df,
      value = value,
      top_n = top_n,
      title = title,
      subtitle = subtitle,
      palette = palette_cfg$link_palette,
      palcolor = palette_cfg$link_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args,
      grid_major = grid_major,
      grid_major_colour = grid_major_colour,
      grid_major_linetype = grid_major_linetype,
      grid_major_linewidth = grid_major_linewidth
    )))
  }

  if (identical(plot_type, "comparison")) {
    return(finish_plot(ccc_stat_comparison_plot(
      srt = srt,
      method = method,
      condition = condition,
      comparison = comparison,
      sender.use = sender.use,
      receiver.use = receiver.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      interaction.use = interaction.use,
      signaling = signaling,
      pairLR.use = pairLR.use,
      thresh = thresh,
      sample_col = sample_col,
      measure = measure,
      compare_by = compare_by,
      pattern = pattern,
      title = title,
      subtitle = subtitle,
      palette = palette_cfg$cell_palette,
      palcolor = palette_cfg$cell_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    )))
  }

  # --- bar ---
  if (identical(plot_type, "bar")) {
    return(finish_plot(ccc_stat_bar_plot(
      df = df,
      pair_df = pair_df,
      stat_type = stat_type,
      display_by = display_by,
      value = value,
      top_n = top_n,
      title = title,
      subtitle = subtitle,
      link_palette = palette_cfg$link_palette,
      link_palcolor = palette_cfg$link_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args,
      grid_major = grid_major,
      grid_major_colour = grid_major_colour,
      grid_major_linetype = grid_major_linetype,
      grid_major_linewidth = grid_major_linewidth
    )))
  }

  # --- sankey ---
  if (identical(plot_type, "sankey")) {
    return(finish_plot(ccc_sankey_plot(
      pair_df = pair_df,
      interaction_df = interaction_df,
      display_by = display_by,
      top_n = top_n,
      edge_value = edge_value,
      min_receiver_flow = min_receiver_flow,
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      link_palette = palette_cfg$link_palette,
      link_palcolor = palette_cfg$link_palcolor,
      title = title,
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    )))
  }

  # --- box / violin ---
  if (plot_type %in% c("box", "violin")) {
    return(finish_plot(ccc_stat_distribution_plot(
      interaction_df = interaction_df,
      plot_type = plot_type,
      top_n = top_n,
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      facet_by = facet_by,
      x_text_angle = x_text_angle,
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      title = title,
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    )))
  }

  # --- scatter ---
  if (identical(plot_type, "scatter")) {
    return(finish_plot(ccc_scatter_plot(
      srt = srt,
      method = method,
      pair_df = pair_df,
      condition = condition,
      dataset = dataset,
      signaling = signaling,
      title = title,
      subtitle = subtitle,
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      link_palette = palette_cfg$link_palette,
      link_palcolor = palette_cfg$link_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args,
      ...
    )))
  }

  log_message(
    "Unsupported {.arg plot_type}: {.val {plot_type}}",
    message_type = "error"
  )
}

ccc_cellchat_subset_df <- function(
  object,
  signaling = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  pairLR.use = NULL,
  thresh = 0.05,
  dataset = NULL,
  slot.name = "net"
) {
  df <- subset_cc_table(
    object = object,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = resolve_group_index_single(object, sender.use),
    targets.use = resolve_group_index_single(object, receiver.use),
    thresh = thresh,
    dataset = dataset
  )
  df <- standardize_long_df(df)
  if (nrow(df) == 0L) {
    return(df)
  }
  if (!"pathway_name" %in% colnames(df)) {
    df$pathway_name <- if ("signaling" %in% colnames(df)) {
      as.character(df$signaling)
    } else if (!is.null(signaling) && length(signaling) == 1L) {
      as.character(signaling)
    } else {
      NA_character_
    }
  }
  if (!"interaction_label" %in% colnames(df)) {
    df$interaction_label <- if ("interaction_name_2" %in% colnames(df)) {
      as.character(df$interaction_name_2)
    } else if ("interaction_name" %in% colnames(df)) {
      as.character(df$interaction_name)
    } else {
      paste(df$ligand, df$receptor, sep = " - ")
    }
  }
  if (!"prob" %in% colnames(df)) {
    df$prob <- suppressWarnings(as.numeric(df$score))
  }
  df
}

ccc_cellchat_weight_col <- function(df) {
  if ("prob" %in% colnames(df)) {
    return("prob")
  }
  if ("score" %in% colnames(df)) {
    return("score")
  }
  NULL
}

ccc_cellchat_pathway_summary_df <- function(
  object,
  dataset,
  signaling = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  thresh = 0.05
) {
  df <- ccc_cellchat_subset_df(
    object = object,
    signaling = signaling,
    sender.use = sender.use,
    receiver.use = receiver.use,
    thresh = thresh,
    dataset = dataset
  )
  if (nrow(df) == 0L) {
    return(data.frame(
      dataset = character(0),
      pathway = character(0),
      value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  weight_col <- ccc_cellchat_weight_col(df)
  if (is.null(weight_col)) {
    log_message(
      "No communication weight column was found in CellChat results",
      message_type = "error"
    )
  }
  agg <- stats::aggregate(
    df[[weight_col]],
    by = list(pathway = as.character(df$pathway_name)),
    FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
  )
  colnames(agg)[colnames(agg) == "x"] <- "value"
  agg$dataset <- dataset
  agg[, c("dataset", "pathway", "value"), drop = FALSE]
}

ccc_lr_contribution_plot <- function(
  srt,
  condition = NULL,
  dataset = 1,
  signaling,
  sender.use = NULL,
  receiver.use = NULL,
  thresh = 0.05,
  top_n = 20,
  return.data = FALSE,
  title = NULL,
  subtitle = NULL,
  link_palette = "Dark2",
  link_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  font.size = 10,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL
) {
  info <- get_dataset_object(srt, condition = condition, dataset = dataset)
  signaling_use <- unique(as.character(signaling))
  contrib_list <- lapply(signaling_use, function(sig) {
    df <- ccc_cellchat_subset_df(
      object = info$object,
      signaling = sig,
      sender.use = sender.use,
      receiver.use = receiver.use,
      thresh = thresh,
      dataset = info$label
    )
    if (nrow(df) == 0L) {
      return(NULL)
    }
    weight_col <- ccc_cellchat_weight_col(df)
    agg <- stats::aggregate(
      df[[weight_col]],
      by = list(interaction = as.character(df$interaction_label)),
      FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    colnames(agg)[colnames(agg) == "x"] <- "value"
    agg <- agg[is.finite(agg$value) & agg$value > 0, , drop = FALSE]
    if (nrow(agg) == 0L) {
      return(NULL)
    }
    agg <- agg[order(agg$value, decreasing = TRUE), , drop = FALSE]
    if (
      is.numeric(top_n) &&
        length(top_n) == 1L &&
        top_n > 0L &&
        nrow(agg) > top_n
    ) {
      agg <- agg[seq_len(top_n), , drop = FALSE]
    }
    agg$contribution <- agg$value / sum(agg$value, na.rm = TRUE)
    agg$signaling <- sig
    agg$dataset <- info$label
    agg
  })
  contrib_df <- do.call(rbind, Filter(Negate(is.null), contrib_list))
  if (is.null(contrib_df) || nrow(contrib_df) == 0L) {
    log_message(
      "No ligand-receptor contribution data available for the selected signaling pathway/pathways",
      message_type = "error"
    )
  }
  rownames(contrib_df) <- NULL
  if (isTRUE(return.data)) {
    return(contrib_df)
  }

  plot_list <- lapply(
    split(contrib_df, contrib_df$signaling),
    function(df_sig) {
      df_sig <- df_sig[
        order(df_sig$contribution, decreasing = TRUE),
        ,
        drop = FALSE
      ]
      df_sig$interaction <- factor(
        df_sig$interaction,
        levels = rev(unique(df_sig$interaction))
      )
      cols <- palette_colors(
        unique(as.character(df_sig$interaction)),
        palette = link_palette,
        palcolor = link_palcolor
      )
      p <- ggplot2::ggplot(
        df_sig,
        ggplot2::aes(x = interaction, y = contribution, fill = interaction)
      ) +
        ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
        ggplot2::scale_y_continuous(
          labels = scales::label_percent(accuracy = 1),
          expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::labs(x = NULL, y = "Contribution")
      finalize_cc_plot(
        p,
        title = if (length(signaling_use) == 1L) {
          title %||% df_sig$signaling[1]
        } else {
          df_sig$signaling[1]
        },
        subtitle = if (length(signaling_use) == 1L) subtitle else NULL,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme_use = theme_use,
        theme_args = theme_args,
        font.size = font.size
      )
    }
  )
  out <- plot_cc_list(plot_list, combine = combine, nrow = nrow, ncol = ncol)
  if (inherits(out, "patchwork")) {
    out <- patchwork_cc(
      out,
      title = if (length(signaling_use) > 1L) title %||% info$label else NULL,
      subtitle = if (length(signaling_use) > 1L) subtitle else NULL,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size
    )
  }
  out
}

ccc_cellchat_stat_comparison <- function(
  srt,
  condition = NULL,
  comparison = c(1, 2),
  min_n = 1L
) {
  cmp <- .cc_get_cmp(srt = srt, condition = condition)
  comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)
  if (length(comp_idx) < min_n) {
    log_message(
      "At least {min_n} dataset(s) are required for this CellChat comparison plot",
      message_type = "error"
    )
  }
  object_names <- names(cmp$object.list)[comp_idx]
  list(
    cmp = cmp,
    comp_idx = comp_idx,
    object_names = object_names,
    object_list = cmp$object.list[object_names]
  )
}

ccc_cellchat_ranknet_plot <- function(
  srt,
  condition = NULL,
  comparison = c(1, 2),
  signaling = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  thresh = 0.05,
  return.data = FALSE,
  title = NULL,
  subtitle = NULL,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  font.size = 10
) {
  cc_cmp <- ccc_cellchat_stat_comparison(
    srt = srt,
    condition = condition,
    comparison = comparison,
    min_n = 2L
  )
  if (length(cc_cmp$comp_idx) != 2L) {
    log_message(
      "{.val plot_type = 'ranknet'} currently requires exactly two datasets",
      message_type = "error"
    )
  }
  df_list <- lapply(seq_along(cc_cmp$comp_idx), function(i) {
    ccc_cellchat_pathway_summary_df(
      object = cc_cmp$object_list[[i]],
      dataset = cc_cmp$object_names[i],
      signaling = signaling,
      sender.use = sender.use,
      receiver.use = receiver.use,
      thresh = thresh
    )
  })
  pathway_union <- unique(unlist(lapply(df_list, function(x) x$pathway)))
  pathway_union <- pathway_union[!is.na(pathway_union)]
  if (length(pathway_union) == 0L) {
    log_message(
      "No signaling pathways are available for the selected comparison",
      message_type = "error"
    )
  }
  wide <- data.frame(pathway = pathway_union, stringsAsFactors = FALSE)
  for (i in seq_along(df_list)) {
    cur <- df_list[[i]]
    value <- rep(0, length(pathway_union))
    if (nrow(cur) > 0L) {
      idx <- match(cur$pathway, pathway_union)
      value[idx[!is.na(idx)]] <- cur$value[!is.na(idx)]
    }
    wide[[paste0("value_", i)]] <- value
    wide[[paste0("dataset_", i)]] <- cc_cmp$object_names[i]
  }
  wide$diff <- wide$value_2 - wide$value_1
  wide$rank_score <- pmax(wide$value_1, wide$value_2, na.rm = TRUE) +
    abs(wide$diff)
  wide <- wide[order(wide$rank_score, decreasing = TRUE), , drop = FALSE]
  wide$pathway <- factor(wide$pathway, levels = rev(unique(wide$pathway)))

  long_df <- rbind(
    data.frame(
      pathway = wide$pathway,
      dataset = cc_cmp$object_names[1],
      value = wide$value_1,
      stringsAsFactors = FALSE
    ),
    data.frame(
      pathway = wide$pathway,
      dataset = cc_cmp$object_names[2],
      value = wide$value_2,
      stringsAsFactors = FALSE
    )
  )
  if (isTRUE(return.data)) {
    return(data.frame(
      pathway = as.character(wide$pathway),
      dataset_1 = cc_cmp$object_names[1],
      value_1 = wide$value_1,
      dataset_2 = cc_cmp$object_names[2],
      value_2 = wide$value_2,
      diff = wide$diff,
      stringsAsFactors = FALSE
    ))
  }

  cols <- palette_colors(
    cc_cmp$object_names,
    palette = palette,
    palcolor = palcolor
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = wide,
      ggplot2::aes(
        y = pathway,
        yend = pathway,
        x = value_1,
        xend = value_2
      ),
      color = "grey75",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = long_df,
      ggplot2::aes(x = value, y = pathway, color = dataset),
      size = 3
    ) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = "Information flow", y = NULL, color = NULL)
  finalize_cc_plot(
    p,
    title = title %||% cmp_cc_label(cc_cmp$cmp, cc_cmp$comp_idx),
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}

ccc_cellchat_role_components <- function(object, signaling = NULL) {
  outgoing <- ccc_cellchat_role_matrix(
    object = object,
    signaling = signaling,
    pattern = "outgoing",
    scale_rows = FALSE
  )$raw
  incoming <- ccc_cellchat_role_matrix(
    object = object,
    signaling = signaling,
    pattern = "incoming",
    scale_rows = FALSE
  )$raw
  list(outgoing = outgoing, incoming = incoming)
}

ccc_cellchat_role_change_plot <- function(
  srt,
  condition = NULL,
  comparison = c(1, 2),
  idents.use,
  signaling = NULL,
  signaling.exclude = NULL,
  return.data = FALSE,
  title = NULL,
  subtitle = NULL,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  aspect.ratio = NULL,
  font.size = 10
) {
  cc_cmp <- ccc_cellchat_stat_comparison(
    srt = srt,
    condition = condition,
    comparison = comparison,
    min_n = 2L
  )
  if (length(cc_cmp$comp_idx) != 2L) {
    log_message(
      "{.val plot_type = 'role_change'} currently requires exactly two datasets",
      message_type = "error"
    )
  }
  role_list <- lapply(cc_cmp$object_list, function(object) {
    ccc_cellchat_role_components(object, signaling = signaling)
  })
  safe_role_vector <- function(mat, pathways, celltype) {
    out <- stats::setNames(rep(0, length(pathways)), pathways)
    if (
      is.null(mat) ||
        !is.matrix(mat) ||
        !length(pathways) ||
        !celltype %in% colnames(mat)
    ) {
      return(out)
    }
    keep <- intersect(pathways, rownames(mat))
    if (length(keep) == 0L) {
      return(out)
    }
    vals <- mat[keep, celltype, drop = TRUE]
    vals <- as.numeric(vals)
    names(vals) <- keep
    out[keep] <- vals
    out[is.na(out)] <- 0
    out
  }
  pathway_union <- union(
    rownames(role_list[[1]]$outgoing),
    rownames(role_list[[2]]$outgoing)
  )
  if (!is.null(signaling.exclude)) {
    pathway_union <- setdiff(pathway_union, as.character(signaling.exclude))
  }
  if (length(pathway_union) == 0L) {
    log_message(
      "No signaling pathways remain after filtering for role-change plotting",
      message_type = "error"
    )
  }
  idents_use <- unique(as.character(idents.use))
  role_df <- do.call(
    rbind,
    lapply(idents_use, function(celltype) {
      out1 <- safe_role_vector(role_list[[1]]$outgoing, pathway_union, celltype)
      out2 <- safe_role_vector(role_list[[2]]$outgoing, pathway_union, celltype)
      in1 <- safe_role_vector(role_list[[1]]$incoming, pathway_union, celltype)
      in2 <- safe_role_vector(role_list[[2]]$incoming, pathway_union, celltype)
      data.frame(
        celltype = celltype,
        signaling = pathway_union,
        outgoing_1 = as.numeric(out1),
        outgoing_2 = as.numeric(out2),
        incoming_1 = as.numeric(in1),
        incoming_2 = as.numeric(in2),
        stringsAsFactors = FALSE
      )
    })
  )
  role_df$outgoing_diff <- role_df$outgoing_2 - role_df$outgoing_1
  role_df$incoming_diff <- role_df$incoming_2 - role_df$incoming_1
  role_df$magnitude <- sqrt(role_df$outgoing_diff^2 + role_df$incoming_diff^2)
  role_df$direction <- ifelse(
    role_df$outgoing_diff >= 0 & role_df$incoming_diff >= 0,
    paste0(cc_cmp$object_names[2], " up"),
    ifelse(
      role_df$outgoing_diff <= 0 & role_df$incoming_diff <= 0,
      paste0(cc_cmp$object_names[1], " up"),
      "mixed"
    )
  )
  if (isTRUE(return.data)) {
    return(role_df)
  }
  cell_cols <- palette_colors(
    idents_use,
    palette = palette,
    palcolor = palcolor
  )
  p <- ggplot2::ggplot(
    role_df,
    ggplot2::aes(
      x = outgoing_diff,
      y = incoming_diff,
      size = magnitude
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    ggplot2::geom_point(
      ggplot2::aes(fill = celltype),
      shape = 21,
      color = "grey20",
      alpha = 0.85,
      stroke = 0.25
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = signaling),
      size = max(2.8, font.size * 0.28),
      min.segment.length = 0,
      box.padding = 0.3,
      point.padding = 0.15,
      seed = 11,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = cell_cols, drop = FALSE) +
    ggplot2::scale_size_continuous(name = "Change", range = c(2.5, 8)) +
    ggplot2::labs(
      x = paste0(
        "Outgoing change (",
        cc_cmp$object_names[2],
        " - ",
        cc_cmp$object_names[1],
        ")"
      ),
      y = paste0(
        "Incoming change (",
        cc_cmp$object_names[2],
        " - ",
        cc_cmp$object_names[1],
        ")"
      ),
      fill = NULL
    )
  if (length(idents_use) > 1L) {
    p <- p + ggplot2::facet_wrap(~celltype)
  }
  finalize_cc_plot(
    p,
    title = title %||% cmp_cc_label(cc_cmp$cmp, cc_cmp$comp_idx),
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    aspect.ratio = aspect.ratio,
    font.size = font.size
  )
}

ccc_match_seurat_genes <- function(genes, seurat_genes) {
  genes <- unique(stats::na.omit(as.character(genes)))
  if (length(genes) == 0L) {
    return(character(0))
  }
  seurat_genes_upper <- toupper(seurat_genes)
  hit_direct <- genes[toupper(genes) %in% seurat_genes_upper]
  miss <- genes[!(toupper(genes) %in% seurat_genes_upper)]
  miss_trim <- gsub("[-_.].*$", "", miss)
  hit_trim <- miss_trim[toupper(miss_trim) %in% seurat_genes_upper]
  match_idx <- match(toupper(c(hit_direct, hit_trim)), seurat_genes_upper)
  unique(stats::na.omit(seurat_genes[match_idx]))
}

ccc_split_gene_tokens <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(character(0))
  }
  unique(unlist(strsplit(x, "_", fixed = TRUE), use.names = FALSE))
}

ccc_rank_token_table <- function(x, weight = NULL) {
  tokens <- lapply(as.character(x), ccc_split_gene_tokens)
  if (length(tokens) == 0L) {
    return(data.frame(gene = character(0), weight = numeric(0)))
  }
  if (is.null(weight)) {
    weight <- rep(1, length(tokens))
  }
  pieces <- lapply(seq_along(tokens), function(i) {
    tok <- tokens[[i]]
    if (length(tok) == 0L) {
      return(NULL)
    }
    data.frame(
      gene = tok,
      weight = rep(weight[i], length(tok)),
      stringsAsFactors = FALSE
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0L) {
    return(data.frame(gene = character(0), weight = numeric(0)))
  }
  out <- stats::aggregate(weight ~ gene, data = do.call(rbind, pieces), FUN = sum, na.rm = TRUE)
  out[order(out$weight, decreasing = TRUE), , drop = FALSE]
}

ccc_standardize_generic_ligand_target <- function(df) {
  df <- standardize_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
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
    } else {
      df$weight <- 1
    }
  }
  if (!"ligand" %in% colnames(df)) {
    df$ligand <- NA_character_
  }
  if (!"target" %in% colnames(df)) {
    df$target <- NA_character_
  }
  if (!"sender" %in% colnames(df)) {
    df$sender <- NA_character_
  }
  if (!"receiver" %in% colnames(df)) {
    df$receiver <- NA_character_
  }
  df$weight <- suppressWarnings(as.numeric(df$weight))
  df
}

ccc_generic_lr_contribution_plot <- function(
  srt,
  method,
  signaling,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  pairLR.use = NULL,
  top_n = 20,
  value = "score",
  return.data = FALSE,
  title = NULL,
  subtitle = NULL,
  link_palette = "Dark2",
  link_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  font.size = 10,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL
) {
  long_df <- ccc_long_table_for_method(
    srt = srt,
    method = method
  )
  long_df <- standardize_long_df(long_df)
  available_pathways <- unique(as.character(long_df$pathway_name))
  available_pathways <- available_pathways[
    !is.na(available_pathways) & nzchar(available_pathways)
  ]
  if (length(available_pathways) == 0L) {
    log_message(
      paste0(
        "{.val plot_type = 'lr_contribution'} requires pathway annotations ",
        "(for example {.code pathway_name} / {.code classification}) in the ",
        "stored CCC result table"
      ),
      message_type = "error"
    )
  }

  long_df <- filter_long_df(
    df = long_df,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    signaling = signaling,
    pairLR.use = pairLR.use
  )
  if (nrow(long_df) == 0L) {
    log_message(
      "No ligand-receptor contribution data available for the selected signaling pathway/pathways",
      message_type = "error"
    )
  }

  score_col <- value
  if (!score_col %in% colnames(long_df)) {
    score_col <- c("score", "means", "prob")[
      c("score", "means", "prob") %in% colnames(long_df)
    ][1]
  }
  if (is.null(score_col) || is.na(score_col)) {
    log_message(
      "No communication weight column was found in the selected CCC results",
      message_type = "error"
    )
  }

  long_df <- prepare_plot_df(long_df)
  long_df$score_use <- suppressWarnings(as.numeric(long_df[[score_col]]))
  long_df <- long_df[is.finite(long_df$score_use) & long_df$score_use > 0, , drop = FALSE]
  if (nrow(long_df) == 0L) {
    log_message(
      "No positive pathway-specific interaction scores remain after filtering",
      message_type = "error"
    )
  }

  signaling_use <- unique(as.character(signaling))
  contrib_list <- lapply(signaling_use, function(sig) {
    df_sig <- long_df[long_df$pathway_name %in% sig, , drop = FALSE]
    if (nrow(df_sig) == 0L) {
      return(NULL)
    }
    agg <- stats::aggregate(
      df_sig$score_use,
      by = list(interaction = as.character(df_sig$interaction_label)),
      FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
    )
    colnames(agg)[colnames(agg) == "x"] <- "value"
    agg <- agg[is.finite(agg$value) & agg$value > 0, , drop = FALSE]
    if (nrow(agg) == 0L) {
      return(NULL)
    }
    agg <- agg[order(agg$value, decreasing = TRUE), , drop = FALSE]
    if (
      is.numeric(top_n) &&
        length(top_n) == 1L &&
        top_n > 0L &&
        nrow(agg) > top_n
    ) {
      agg <- agg[seq_len(top_n), , drop = FALSE]
    }
    agg$contribution <- agg$value / sum(agg$value, na.rm = TRUE)
    agg$signaling <- sig
    agg
  })
  contrib_df <- do.call(rbind, Filter(Negate(is.null), contrib_list))
  if (is.null(contrib_df) || nrow(contrib_df) == 0L) {
    log_message(
      "No ligand-receptor contribution data available for the selected signaling pathway/pathways",
      message_type = "error"
    )
  }
  rownames(contrib_df) <- NULL
  if (isTRUE(return.data)) {
    return(contrib_df)
  }

  plot_list <- lapply(
    split(contrib_df, contrib_df$signaling),
    function(df_sig) {
      df_sig <- df_sig[
        order(df_sig$contribution, decreasing = TRUE),
        ,
        drop = FALSE
      ]
      df_sig$interaction <- factor(
        df_sig$interaction,
        levels = rev(unique(df_sig$interaction))
      )
      cols <- palette_colors(
        unique(as.character(df_sig$interaction)),
        palette = link_palette,
        palcolor = link_palcolor
      )
      p <- ggplot2::ggplot(
        df_sig,
        ggplot2::aes(x = interaction, y = contribution, fill = interaction)
      ) +
        ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
        ggplot2::scale_y_continuous(
          labels = scales::label_percent(accuracy = 1),
          expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::labs(x = NULL, y = "Contribution")
      finalize_cc_plot(
        p,
        title = if (length(signaling_use) == 1L) {
          title %||% df_sig$signaling[1]
        } else {
          df_sig$signaling[1]
        },
        subtitle = if (length(signaling_use) == 1L) subtitle else NULL,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme_use = theme_use,
        theme_args = theme_args,
        font.size = font.size
      )
    }
  )

  if (!isTRUE(combine) || length(plot_list) <= 1L) {
    return(simplify_cc_plot_list(plot_list))
  }
  patchwork::wrap_plots(plot_list, nrow = nrow, ncol = ncol)
}

ccc_generic_feature_genes <- function(
  srt,
  method,
  condition = NULL,
  dataset = 1,
  signaling = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  pairLR.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  top_n = 20
) {
  long_df <- ccc_long_table_for_method(
    srt = srt,
    method = method,
    condition = condition,
    dataset = dataset,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sender.use,
    targets.use = receiver.use,
    thresh = thresh
  )
  long_df <- standardize_long_df(long_df)
  long_df <- filter_long_df(
    df = long_df,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    signaling = signaling,
    pairLR.use = pairLR.use
  )

  gene_panels <- list()
  if (nrow(long_df) > 0L) {
    ligand_tbl <- ccc_rank_token_table(long_df$ligand)
    receptor_tbl <- ccc_rank_token_table(long_df$receptor)
    if (nrow(ligand_tbl) > 0L) {
      gene_panels[["Ligand"]] <- utils::head(ligand_tbl$gene, top_n)
    }
    if (nrow(receptor_tbl) > 0L) {
      gene_panels[["Receptor"]] <- utils::head(receptor_tbl$gene, top_n)
    }
  }

  target_df <- ccc_standardize_generic_ligand_target(bundle$ligand_target_df %||% data.frame())
  if (nrow(target_df) > 0L) {
    if (!all(is.na(target_df$sender)) && !is.null(sender.use)) {
      target_df <- target_df[target_df$sender %in% sender.use, , drop = FALSE]
    }
    if (!all(is.na(target_df$receiver)) && !is.null(receiver.use)) {
      target_df <- target_df[target_df$receiver %in% receiver.use, , drop = FALSE]
    }
    if (!is.null(ligand.use)) {
      target_df <- target_df[target_df$ligand %in% ligand.use, , drop = FALSE]
    }
    if (nrow(target_df) > 0L) {
      target_tbl <- ccc_rank_token_table(target_df$target, weight = target_df$weight %||% NULL)
      ligand_target_tbl <- ccc_rank_token_table(target_df$ligand, weight = target_df$weight %||% NULL)
      if (nrow(target_tbl) > 0L) {
        gene_panels[["Target"]] <- utils::head(target_tbl$gene, top_n)
      }
      if (identical(method, "MultiNichenetr") && nrow(ligand_target_tbl) > 0L && is.null(gene_panels[["Ligand"]])) {
        gene_panels[["Ligand"]] <- utils::head(ligand_target_tbl$gene, top_n)
      }
    }
  }

  if (length(gene_panels) == 0L) {
    return(data.frame(panel = character(0), gene = character(0)))
  }

  gene_df <- do.call(rbind, lapply(names(gene_panels), function(panel) {
    genes <- ccc_match_seurat_genes(gene_panels[[panel]], rownames(srt))
    if (length(genes) == 0L) {
      return(NULL)
    }
    data.frame(
      panel = panel,
      gene = genes,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(gene_df)) {
    gene_df <- data.frame(panel = character(0), gene = character(0))
  }
  rownames(gene_df) <- NULL
  gene_df
}

ccc_cellchat_signaling_genes <- function(
  object,
  seurat_object,
  signaling,
  sender.use = NULL,
  receiver.use = NULL,
  thresh = 0.05
) {
  split_gene_tokens <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
      return(character(0))
    }
    unique(unlist(strsplit(x, "_", fixed = TRUE), use.names = FALSE))
  }

  df <- ccc_cellchat_subset_df(
    object = object,
    signaling = signaling,
    sender.use = sender.use,
    receiver.use = receiver.use,
    thresh = thresh
  )
  if (nrow(df) == 0L) {
    return(character(0))
  }
  gene_pool <- unique(c(
    split_gene_tokens(df$ligand),
    split_gene_tokens(df$receptor),
    split_gene_tokens(df$interaction_name)
  ))
  ccc_match_seurat_genes(gene_pool, rownames(seurat_object))
}

ccc_feature_plot <- function(
  srt,
  method = "CellChat",
  condition = NULL,
  dataset = 1,
  signaling = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  pairLR.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  top_n = 20,
  return.data = FALSE,
  title = NULL,
  subtitle = NULL,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  aspect.ratio = NULL,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  dots = list()
) {
  if (identical(method, "CellChat")) {
    info <- get_dataset_object(srt, condition = condition, dataset = dataset)
    signaling_use <- unique(as.character(signaling %||% info$object@netP$pathways))
    if (
      is.numeric(top_n) &&
        length(top_n) == 1L &&
        top_n > 0L &&
        length(signaling_use) > top_n
    ) {
      signaling_use <- signaling_use[seq_len(top_n)]
    }
    gene_df_list <- lapply(signaling_use, function(sig) {
      genes <- ccc_cellchat_signaling_genes(
        object = info$object,
        seurat_object = srt,
        signaling = sig,
        sender.use = sender.use,
        receiver.use = receiver.use,
        thresh = thresh
      )
      if (length(genes) == 0L) {
        return(NULL)
      }
      data.frame(
        panel = sig,
        gene = genes,
        stringsAsFactors = FALSE
      )
    })
    gene_df <- do.call(rbind, Filter(Negate(is.null), gene_df_list))
    if (is.null(gene_df) || nrow(gene_df) == 0L) {
      log_message(
        "No signaling genes were matched in the Seurat object for the selected CellChat pathway/pathways",
        message_type = "error"
      )
    }
    rownames(gene_df) <- NULL
  } else {
    gene_df <- ccc_generic_feature_genes(
      srt = srt,
      method = method,
      condition = condition,
      dataset = dataset,
      signaling = signaling,
      sender.use = sender.use,
      receiver.use = receiver.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      slot.name = slot.name,
      thresh = thresh,
      top_n = top_n
    )
    if (nrow(gene_df) == 0L) {
      log_message(
        "No ligand/receptor/target genes were matched in the Seurat object for the selected CCC filters",
        message_type = "error"
      )
    }
    signaling_use <- unique(gene_df$panel)
  }

  if (isTRUE(return.data)) {
    return(gene_df)
  }

  group.by <- get_group_by(srt, method = method)
  if (is.null(group.by) || !group.by %in% colnames(srt@meta.data)) {
    group.by <- ".ccc_gene_group"
    srt[[group.by]] <- as.character(SeuratObject::Idents(srt))
  }
  gene_plot_type <- dots[["gene_plot_type"]] %||% "violin"
  gene_split <- split(gene_df$gene, gene_df$panel)
  plot_list <- lapply(names(gene_split), function(panel_name) {
    genes <- unique(unname(gene_split[[panel_name]]))
    local_dots <- dots
    local_dots[c(
      "return.data",
      "idents.use",
      "signaling.exclude",
      "aspect.ratio",
      "gene_plot_type"
    )] <- NULL
    p <- do.call(
      FeatureStatPlot,
      c(
        list(
          srt = srt,
          stat.by = genes,
          group.by = group.by,
          plot.by = "feature",
          plot_type = gene_plot_type,
          palette = palette,
          palcolor = palcolor,
          stack = isTRUE(local_dots[["stack"]]) ||
            is.null(local_dots[["stack"]]),
          legend.position = legend.position,
          legend.direction = legend.direction,
          theme_use = theme_use,
          theme_args = theme_args,
          aspect.ratio = aspect.ratio,
          title = if (length(signaling_use) == 1L) title %||% panel_name else panel_name,
          subtitle = if (length(signaling_use) == 1L) subtitle else NULL
        ),
        local_dots
      )
    )
    p
  })
  plot_cc_list(plot_list, combine = combine, nrow = nrow, ncol = ncol)
}

ccc_stat_bar_plot <- function(
  df,
  pair_df,
  stat_type = "score",
  display_by = "aggregation",
  value = "sum",
  top_n = 20,
  title = NULL,
  subtitle = NULL,
  link_palette = "Dark2",
  link_palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3
) {
  grid_major_element <- if (isTRUE(grid_major)) {
    ggplot2::element_line(
      colour = grid_major_colour,
      linetype = grid_major_linetype,
      linewidth = grid_major_linewidth
    )
  } else {
    ggplot2::element_blank()
  }
  display_by <- match.arg(display_by, c("aggregation", "interaction"))
  value_use <- if (identical(stat_type, "count")) "count" else value
  if (identical(display_by, "aggregation")) {
    if (is.null(pair_df) || nrow(pair_df) == 0L) {
      log_message(
        "No aggregated sender-receiver communication records remain after filtering",
        message_type = "error"
      )
    }
    if (!"pair" %in% colnames(pair_df)) {
      pair_df$pair <- paste(pair_df$sender, pair_df$receiver, sep = " -> ")
    }
    value_col <- value_use
    if (!value_col %in% colnames(pair_df)) {
      value_col <- c("sum", "score", "mean", "max", "count")[
        c("sum", "score", "mean", "max", "count") %in% colnames(pair_df)
      ][1]
    }
    if (is.null(value_col) || is.na(value_col)) {
      log_message(
        "No communication score column was found in aggregated CCC results",
        message_type = "error"
      )
    }
    plot_group <- as.character(pair_df$pair)
    metric <- suppressWarnings(as.numeric(pair_df[[value_col]]))
    value_use <- value_col
  } else {
    if (is.null(df) || nrow(df) == 0L) {
      log_message(
        "No interaction-level communication records remain after filtering",
        message_type = "error"
      )
    }
    df <- ccc_mark_significance(df)
    df <- prepare_plot_df(df)
    plot_group <- as.character(df$interaction_label)
    metric <- if (identical(value_use, "count")) {
      as.numeric(df$significant)
    } else {
      suppressWarnings(as.numeric(df$score))
    }
  }
  plot_group[is.na(plot_group) | !nzchar(plot_group)] <- "Unclassified"
  metric[!is.finite(metric)] <- 0
  summary <- switch(value_use,
    mean = tapply(metric, plot_group, mean, na.rm = TRUE),
    max = tapply(metric, plot_group, max, na.rm = TRUE),
    count = tapply(metric, plot_group, sum, na.rm = TRUE),
    tapply(metric, plot_group, sum, na.rm = TRUE)
  )
  summary <- sort(summary, decreasing = TRUE)
  summary <- utils::head(summary, top_n)
  if (length(summary) == 0L) {
    log_message(
      "No communication groups remain after summarizing",
      message_type = "error"
    )
  }
  plot_df <- data.frame(
    group = names(summary),
    value = as.numeric(summary),
    stringsAsFactors = FALSE
  )
  plot_df$group <- factor(plot_df$group, levels = rev(plot_df$group))
  cols <- palette_colors(
    as.character(plot_df$group),
    palette = link_palette,
    palcolor = link_palcolor
  )
  value_label <- if (identical(value_use, "count")) {
    "Significant interactions"
  } else if (identical(value_use, "mean")) {
    "Mean communication score"
  } else if (identical(value_use, "max")) {
    "Max communication score"
  } else {
    "Communication score"
  }
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = group, y = value, fill = group)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::labs(x = NULL, y = value_label)
  return(finalize_cc_plot(
    p,
    title = title %||% paste0(
      "Top communication ",
      if (identical(display_by, "aggregation")) "pairs" else "interactions"
    ),
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  ) +
    ggplot2::theme(panel.grid.major = grid_major_element))
}

ccc_stat_distribution_plot <- function(
  interaction_df,
  plot_type = c("box", "violin"),
  top_n = 20,
  interaction.use = NULL,
  pairLR.use = NULL,
  facet_by = NULL,
  x_text_angle = 90,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  ccc_distribution_plot(
    interaction_df = interaction_df,
    plot_type = plot_type,
    top_n = top_n,
    interaction.use = interaction.use,
    pairLR.use = pairLR.use,
    facet_by = facet_by,
    x_text_angle = x_text_angle,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    font.size = font.size,
    theme_use = theme_use,
    theme_args = theme_args
  )
}

ccc_pathway_summary_plot <- function(
  df,
  value = "sum",
  top_n = 20,
  title = NULL,
  subtitle = NULL,
  palette = "Dark2",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3
) {
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No communication records remain after filtering",
      message_type = "error"
    )
  }
  df <- ccc_mark_significance(df)
  pathway <- as.character(df$pathway_name)
  pathway[is.na(pathway) | !nzchar(pathway)] <- "Unclassified"
  score <- suppressWarnings(as.numeric(df$score))
  score[!is.finite(score)] <- 0
  metric <- if (identical(value, "count")) {
    as.numeric(df$significant)
  } else {
    score
  }
  total_strength <- switch(value,
    mean = tapply(metric, pathway, mean, na.rm = TRUE),
    max = tapply(metric, pathway, max, na.rm = TRUE),
    count = tapply(metric, pathway, sum, na.rm = TRUE),
    tapply(metric, pathway, sum, na.rm = TRUE)
  )
  sig_pairs <- tapply(as.numeric(df$significant), pathway, sum, na.rm = TRUE)
  active_pairs <- tapply(score > 0, pathway, sum, na.rm = TRUE)
  summary <- data.frame(
    pathway = names(total_strength),
    total_strength = as.numeric(total_strength),
    n_significant_pairs = as.numeric(sig_pairs[names(total_strength)]),
    n_active_cell_pairs = as.numeric(active_pairs[names(total_strength)]),
    stringsAsFactors = FALSE
  )
  summary$n_significant_pairs[is.na(summary$n_significant_pairs)] <- 0
  summary$n_active_cell_pairs[is.na(summary$n_active_cell_pairs)] <- 0
  summary$is_significant <- summary$n_significant_pairs > 0
  summary <- summary[
    order(summary$is_significant, summary$total_strength, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  summary <- utils::head(summary, top_n)
  if (nrow(summary) == 0L) {
    log_message(
      "No pathway communication records remain after summarizing",
      message_type = "error"
    )
  }
  summary$pathway <- factor(summary$pathway, levels = rev(summary$pathway))
  cols <- palette_colors(
    as.character(summary$pathway),
    palette = palette,
    palcolor = palcolor
  )
  bar_cols <- cols[as.character(summary$pathway)]
  bar_cols[!summary$is_significant] <- "#D9D9D9"
  summary$label <- paste0(
    as.integer(summary$n_significant_pairs),
    "/",
    as.integer(summary$n_active_cell_pairs),
    " sig"
  )
  grid_major_element <- if (isTRUE(grid_major)) {
    ggplot2::element_line(
      colour = grid_major_colour,
      linetype = grid_major_linetype,
      linewidth = grid_major_linewidth
    )
  } else {
    ggplot2::element_blank()
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = pathway, y = total_strength, fill = pathway)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust = -0.08,
      size = max(3, font.size * 0.28)
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_fill_manual(values = bar_cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(x = NULL, y = "Total pathway communication strength")

  finalize_cc_plot(
    p,
    title = title %||% "Significant pathway communication summary",
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  ) +
    ggplot2::theme(panel.grid.major = grid_major_element)
}

ccc_stat_comparison_plot <- function(
  srt,
  method,
  condition = NULL,
  comparison = c(1, 2),
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  thresh = 0.05,
  sample_col = NULL,
  measure = "count",
  compare_by = "overall",
  pattern = "all",
  title = NULL,
  subtitle = NULL,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  if (!identical(method, "CellChat")) {
    return(ccc_generic_stat_comparison_plot(
      srt = srt,
      method = method,
      comparison = comparison,
      sender.use = sender.use,
      receiver.use = receiver.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      interaction.use = interaction.use,
      signaling = signaling,
      pairLR.use = pairLR.use,
      thresh = thresh,
      sample_col = sample_col,
      measure = measure,
      compare_by = compare_by,
      pattern = pattern,
      title = title,
      subtitle = subtitle,
      palette = palette,
      palcolor = palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }
  cc_cmp <- ccc_cellchat_stat_comparison(
    srt = srt,
    condition = condition,
    comparison = comparison,
    min_n = 1L
  )
  cols <- palette_colors(
    cc_cmp$object_names,
    palette = palette,
    palcolor = palcolor
  )

  if (identical(compare_by, "overall")) {
    df <- do.call(
      rbind,
      lapply(seq_along(cc_cmp$object_list), function(i) {
        obj <- cc_cmp$object_list[[i]]
        value <- switch(
          measure,
          count = sum(obj@net$count, na.rm = TRUE),
          weight = sum(obj@net$weight, na.rm = TRUE)
        )
        data.frame(
          dataset = cc_cmp$object_names[i],
          group = cc_cmp$object_names[i],
          value = value,
          stringsAsFactors = FALSE
        )
      })
    )
    ylab <- if (identical(measure, "count")) {
      "Number of inferred interactions"
    } else {
      "Interaction strength"
    }
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = dataset, y = value, fill = group)
    ) +
      ggplot2::geom_col(width = 0.65) +
      ggplot2::geom_text(
        ggplot2::aes(label = signif(value, 3)),
        vjust = -0.25,
        size = 3
      ) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
      ggplot2::labs(x = NULL, y = ylab, fill = NULL)
    return(finalize_cc_plot(
      p,
      title = title %||% cmp_cc_label(cc_cmp$cmp, cc_cmp$comp_idx),
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size
    ))
  }

  df <- do.call(
    rbind,
    lapply(seq_along(cc_cmp$object_list), function(i) {
      obj <- cc_cmp$object_list[[i]]
      mat <- switch(
        measure,
        count = obj@net$count,
        weight = obj@net$weight
      )
      celltypes <- rownames(mat)
      value <- switch(
        pattern,
        outgoing = rowSums(mat, na.rm = TRUE),
        incoming = colSums(mat, na.rm = TRUE),
        all = rowSums(mat, na.rm = TRUE) +
          colSums(mat, na.rm = TRUE) -
          diag(mat)
      )
      data.frame(
        dataset = cc_cmp$object_names[i],
        group = cc_cmp$object_names[i],
        celltype = celltypes,
        value = as.numeric(value),
        stringsAsFactors = FALSE
      )
    })
  )
  ylab <- if (identical(measure, "count")) {
    "Number of inferred interactions"
  } else {
    "Interaction strength"
  }
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = group, y = value, fill = group)
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = signif(value, 3)),
      vjust = -0.25,
      size = 2.6
    ) +
    ggplot2::facet_wrap(~celltype, scales = "free_y") +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = NULL, y = ylab, fill = NULL)
  finalize_cc_plot(
    p,
    title = title %||%
      paste0(cmp_cc_label(cc_cmp$cmp, cc_cmp$comp_idx), ": ", pattern),
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}

ccc_generic_stat_comparison_plot <- function(
  srt,
  method,
  comparison = c(1, 2),
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  thresh = 0.05,
  sample_col = NULL,
  measure = "count",
  compare_by = "overall",
  pattern = "all",
  title = NULL,
  subtitle = NULL,
  palette = "Chinese",
  palcolor = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list()
) {
  plot_data <- ccc_plot_data(
    srt = srt,
    method = method,
    condition = NULL,
    dataset = 1,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sender.use = sender.use,
    receiver.use = receiver.use,
    ligand.use = ligand.use,
    receptor.use = receptor.use,
    interaction.use = interaction.use,
    value = if (identical(measure, "count")) "count" else "score",
    thresh = thresh
  )
  long_df <- plot_data$long_df
  if (nrow(long_df) == 0L) {
    log_message(
      "No CCC records are available for generic comparison plotting",
      message_type = "error"
    )
  }

  cmp_data <- ccc_context_comparison_data(
    long_df = long_df,
    comparison = comparison,
    sample_col = sample_col,
    measure = measure,
    compare_by = compare_by,
    pattern = pattern
  )
  context_names <- cmp_data$context_names
  cols <- palette_colors(
    context_names,
    palette = palette,
    palcolor = palcolor
  )

  if (identical(compare_by, "overall")) {
    df <- cmp_data$data
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = dataset, y = value, fill = group)
    ) +
      ggplot2::geom_col(width = 0.65) +
      ggplot2::geom_text(
        ggplot2::aes(label = signif(value, 3)),
        vjust = -0.25,
        size = 3
      ) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
      ggplot2::labs(x = NULL, y = cmp_data$ylab, fill = NULL)
    return(finalize_cc_plot(
      p,
      title = title %||% paste(context_names, collapse = " vs "),
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size
    ))
  }

  df <- cmp_data$data
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = group, y = value, fill = group)
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = signif(value, 3)),
      vjust = -0.25,
      size = 2.6
    ) +
    ggplot2::facet_wrap(~celltype, scales = "free_y") +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = NULL, y = cmp_data$ylab, fill = NULL)
  finalize_cc_plot(
    p,
    title = title %||% paste0(paste(context_names, collapse = " vs "), ": ", pattern),
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = font.size
  )
}
