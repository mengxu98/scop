#' @title CCC heatmap and dot matrix plot
#'
#' @md
#' @inheritParams CCCStatPlot
#' @param plot_type Plot type. One of `"heatmap"` or `"dot"`.
#' `"bubble"` is a CellChat-specific interaction bubble matrix.
#' `"ligand_target"` is a special heatmap path available only with
#' `NicheNet`/`MultiNicheNet` results. `"role_heatmap"` and `"diff_heatmap"`
#' are CellChat-specific pathway role views.
#' @param top_anno,bottom_anno Column-side annotations for sender groups.
#' Each side accepts `NULL`, `"bar"`, `"box"`, `"point"`, `"line"`,
#' `"histogram"`, `"density"`, `"violin"`, `"cell"`, or a character vector
#' containing multiple values.
#' Defaults are `top_anno = "bar"` and `bottom_anno = "cell"`.
#' @param left_anno,right_anno Row-side annotations for receiver groups.
#' Each side accepts `NULL`, `"bar"`, `"box"`, `"point"`, `"line"`,
#' `"histogram"`, `"density"`, `"violin"`, `"cell"`, or a character vector
#' containing multiple values.
#' Defaults are `left_anno = "bar"` and `right_anno = "cell"`.
#' @param bar_value Aggregation metric shown in the bar annotations.
#' One or more of `"count"`, `"sum"`, `"mean"`, or `"max"`.
#' Multiple values add multiple annotation tracks. Default `"sum"`.
#' @param add_text Logical. Show numeric value labels inside each cell (heatmap mode only).
#' Default `TRUE` for aggregation mode, `FALSE` for interaction mode.
#' @param cluster_rows,cluster_columns Whether to cluster heatmap rows/columns.
#'   Defaults are both `FALSE`.
#' @param border Logical. Whether to draw borders for the heatmap body and all
#'   annotation tracks. Default `TRUE`.
#' @param value_palette Palette used for heatmap value fills.
#' @param value_palcolor Optional custom colors for `value_palette`.
#' @param width,height Optional heatmap body width and height. When only one is
#'   supplied, the other is inferred from the matrix dimensions to keep cells
#'   square. When both are `NULL`, a square-cell size is computed automatically.
#' @param units Units for `width` and `height`. Default `"inch"`.
#'
#' @return A ggplot / patchwork object wrapping the ComplexHeatmap grob.
#' @export
#' 
#' @examples
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
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "dot",
#'   display_by = "aggregation",
#'   top_n = 20
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "dot",
#'   display_by = "interaction",
#'   facet_by = "sender",
#'   top_n = 20
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "dot",
#'   display_by = "interaction",
#'   facet_by = "receiver",
#'   top_n = 20
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "heatmap",
#'   display_by = "aggregation",
#'   top_n = 20
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "heatmap",
#'   display_by = "interaction",
#'   facet_by = "sender",
#'   top_n = 10
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "heatmap",
#'   display_by = "interaction",
#'   facet_by = "sender",
#'   color.by = "pvalue",
#'   top_n = 10
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "bubble",
#'   top_n = 5
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "bubble",
#'   top_n = 5
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "role_heatmap",
#'   #' pattern = "outgoing",
#'   width = 0.6,
#'   height = 2.5
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "role_heatmap",
#'   pattern = "outgoing",
#'   width = 0.6,
#'   height = 2.5
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "role_heatmap",
#'   pattern = "incoming",
#'   palette = "Paired",
#'   width = 0.6,
#'   height = 3.5
#' )
#'
#' CCCHeatmap(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "diff_heatmap",
#'   pattern = "all",
#'   palette = "Paired",
#'   top_n = 20,
#'   width = 0.6,
#'   height = 3.5
#' )
CCCHeatmap <- function(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c(
    "heatmap",
    "dot",
    "bubble",
    "ligand_target",
    "role_heatmap",
    "diff_heatmap"
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
  pattern = c("outgoing", "incoming", "all"),
  value = "sum",
  top_n = 500,
  top_anno = "bar",
  right_anno = "cell",
  left_anno = "bar",
  bottom_anno = "cell",
  bar_value = "sum",
  add_text = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  color.by = c("score", "pvalue"),
  x_text_angle = 90,
  facet_by = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  edge_value = c("sum", "mean", "max", "count"),
  border = TRUE,
  value_palette = "RdBu",
  value_palcolor = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  palette = "Chinese",
  palcolor = NULL,
  width = NULL,
  height = NULL,
  units = "inch",
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
) {
  value_supplied <- !missing(value)
  edge_value_supplied <- !missing(edge_value)
  bar_value_supplied <- !missing(bar_value)
  dots <- list(...)
  dot_names <- names(dots)
  return.data <- isTRUE(dots[["return.data"]])
  remove.isolate <- if (is.null(dots[["remove.isolate"]])) {
    TRUE
  } else {
    isTRUE(dots[["remove.isolate"]])
  }
  bubble_color.by <- dots[["bubble_color.by"]] %||% "prob"
  bubble_size.range <- dots[["bubble_size.range"]] %||% c(1.5, 8)
  angle.x <- dots[["angle.x"]] %||% x_text_angle
  hjust.x <- dots[["hjust.x"]] %||% 1
  vjust.x <- dots[["vjust.x"]] %||% 1
  aspect.ratio <- dots[["aspect.ratio"]] %||% NULL

  if ("add_dot" %in% dot_names) {
    log_message(
      "{.arg add_dot} has been removed from {.fn CCCHeatmap}. Use {.code plot_type = 'dot'} instead.",
      message_type = "error"
    )
  }

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  display_by <- match.arg(display_by)
  plot_type <- match.arg(plot_type)
  measure <- match.arg(measure)
  pattern <- match.arg(pattern)
  edge_value <- match.arg(edge_value)
  color.by <- match.arg(color.by)
  default_bar_source <- if (isTRUE(value_supplied)) {
    value
  } else if (isTRUE(edge_value_supplied)) {
    edge_value
  } else {
    NULL
  }
  if (
    !isTRUE(bar_value_supplied) &&
      length(default_bar_source) == 1L &&
      !is.na(default_bar_source) &&
      default_bar_source %in% c("count", "sum", "mean", "max")
  ) {
    bar_value <- default_bar_source
  }
  bar_value <- match.arg(
    bar_value,
    choices = c("count", "sum", "mean", "max"),
    several.ok = TRUE
  )
  top_anno <- ccc_match_side_anno(top_anno)
  right_anno <- ccc_match_side_anno(right_anno)
  left_anno <- ccc_match_side_anno(left_anno)
  bottom_anno <- ccc_match_side_anno(bottom_anno)

  if ("add_top_bar" %in% dot_names) {
    top_anno <- ccc_update_side_anno_legacy(
      side_anno = top_anno,
      show_bar = isTRUE(dots[["add_top_bar"]])
    )
  }
  if ("add_right_bar" %in% dot_names) {
    left_anno <- ccc_update_side_anno_legacy(
      side_anno = left_anno,
      show_bar = isTRUE(dots[["add_right_bar"]])
    )
  }

  # Allow palette/palcolor as shorthand for cell_palette/cell_palcolor
  if (is.null(cell_palette) || identical(cell_palette, "Chinese")) {
    cell_palette <- palette
    cell_palcolor <- palcolor
  }

  palette_cfg <- ccc_palettes(
    palette = palette,
    palcolor = palcolor,
    value_palette = value_palette,
    value_palcolor = value_palcolor,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor
  )

  method <- detect_method(srt = srt, method = method)

  if (plot_type %in% c("role_heatmap", "diff_heatmap")) {
    if (!identical(method, "CellChat")) {
      log_message(
        paste0(
          "{.arg plot_type} = {.val {plot_type}} is currently only supported ",
          "for {.pkg CellChat} results"
        ),
        message_type = "error"
      )
    }
    if (identical(plot_type, "role_heatmap")) {
      return(ccc_cellchat_role_heatmap_plot(
        srt = srt,
        condition = condition,
        dataset = dataset,
        comparison = comparison,
        signaling = signaling,
        pattern = pattern,
        top_anno = top_anno,
        right_anno = right_anno,
        left_anno = left_anno,
        bottom_anno = bottom_anno,
        bar_value = bar_value,
        add_text = add_text,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        x_text_angle = x_text_angle,
        border = border,
        width = width,
        height = height,
        units = units,
        title = title,
        subtitle = subtitle,
        value_palette = palette_cfg$value_palette,
        value_palcolor = palette_cfg$value_palcolor,
        cell_palette = palette_cfg$cell_palette,
        cell_palcolor = palette_cfg$cell_palcolor,
        legend.position = legend.position,
        legend.direction = legend.direction,
        font.size = font.size,
        theme_use = theme_use,
        theme_args = theme_args
      ))
    }
    return(ccc_cellchat_diff_heatmap_plot(
      srt = srt,
      condition = condition,
      comparison = comparison,
      signaling = signaling,
      pattern = pattern,
      top_anno = top_anno,
      right_anno = right_anno,
      left_anno = left_anno,
      bottom_anno = bottom_anno,
      bar_value = bar_value,
      add_text = add_text,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      x_text_angle = x_text_angle,
      border = border,
      width = width,
      height = height,
      units = units,
      title = title,
      subtitle = subtitle,
      value_palette = palette_cfg$value_palette,
      value_palcolor = palette_cfg$value_palcolor,
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }

  if (identical(plot_type, "bubble")) {
    if (!identical(method, "CellChat")) {
      log_message(
        "{.val plot_type = 'bubble'} is currently only supported for {.pkg CellChat} results",
        message_type = "error"
      )
    }
    if (use_cc_single_condition(srt, condition = condition)) {
      cond1 <- resolve_single_cc_condition(srt, condition = condition)
      obj1 <- get_single_cc_obj(srt, condition = cond1)
      df_bubble <- prepare_cc_bubble_data(
        object = obj1,
        slot.name = slot.name,
        signaling = signaling,
        pairLR.use = pairLR.use,
        sources.use = resolve_group_index_single(obj1, sender.use),
        targets.use = resolve_group_index_single(obj1, receiver.use),
        thresh = thresh,
        dataset = NULL,
        top_n = top_n,
        remove.isolate = remove.isolate
      )
      if (isTRUE(return.data)) {
        return(df_bubble)
      }
      return(custom_cc_bubble_plot(
        df = df_bubble,
        title = title %||% cond1,
        subtitle = subtitle,
        xlab = NULL,
        ylab = NULL,
        color.by = bubble_color.by,
        bubble_size.range = bubble_size.range,
        palette = palette_cfg$value_palette,
        palcolor = palette_cfg$value_palcolor,
        font.size = font.size,
        angle.x = angle.x,
        hjust.x = hjust.x,
        vjust.x = vjust.x,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme_use = theme_use,
        theme_args = theme_args,
        aspect.ratio = aspect.ratio,
        remove.isolate = remove.isolate
      ))
    }

    cmp <- .cc_get_cmp(srt = srt, condition = condition)
    comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)
    obj.list <- cmp$object.list
    df_list <- lapply(names(obj.list)[comp_idx], function(ds) {
      obj <- obj.list[[ds]]
      prepare_cc_bubble_data(
        object = obj,
        slot.name = slot.name,
        signaling = signaling,
        pairLR.use = pairLR.use,
        sources.use = resolve_group_index_single(obj, sender.use),
        targets.use = resolve_group_index_single(obj, receiver.use),
        thresh = thresh,
        dataset = ds,
        top_n = top_n,
        remove.isolate = remove.isolate
      )
    })
    names(df_list) <- names(obj.list)[comp_idx]
    if (isTRUE(return.data)) {
      return(df_list)
    }
    df_list <- Filter(function(x) !is.null(x) && nrow(x) > 0L, df_list)
    if (length(df_list) == 0L) {
      log_message(
        "No communication records found under the current filtering conditions",
        message_type = "error"
      )
    }
    df_bubble <- do.call(rbind, df_list)
    return(custom_cc_bubble_plot(
      df = df_bubble,
      title = title %||% cmp_cc_label(cmp, comp_idx),
      subtitle = subtitle,
      xlab = NULL,
      ylab = NULL,
      color.by = bubble_color.by,
      bubble_size.range = bubble_size.range,
      palette = palette_cfg$value_palette,
      palcolor = palette_cfg$value_palcolor,
      font.size = font.size,
      angle.x = angle.x,
      hjust.x = hjust.x,
      vjust.x = vjust.x,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio,
      remove.isolate = remove.isolate
    ))
  }

  # --- ligand_target: special path (NicheNet/MultiNicheNet only) ---
  if (identical(plot_type, "ligand_target")) {
    ligand_method <- method
    if (
      !identical(method, "Nichenetr") && !identical(method, "MultiNichenetr")
    ) {
      available <- names(srt@tools)
      if ("Nichenetr" %in% available) {
        ligand_method <- "Nichenetr"
      } else if ("MultiNichenetr" %in% available) {
        ligand_method <- "MultiNichenetr"
      } else {
        log_message(
          "Ligand-target plots require {.pkg NicheNet}/{.pkg MultiNicheNet} results. Run {.fn RunNichenetr} or {.fn RunMultiNichenetr} first.",
          message_type = "error"
        )
      }
    }
    bundle <- get_bundle(srt, method = ligand_method)
    ligand_target_df <- bundle$ligand_target_df %||% data.frame()
    return(ccc_ligand_target_heatmap(
      df = ligand_target_df,
      top_n = top_n,
      sender.use = sender.use,
      receiver.use = receiver.use,
      context_df = bundle$long_table %||% data.frame(),
      sender_default = bundle$parameters$sender_celltypes %||%
        bundle$parameters$sender_use %||%
        NULL,
      receiver_default = bundle$parameters$receiver_celltypes %||%
        bundle$parameters$receiver %||%
        NULL,
      top_anno = top_anno,
      right_anno = right_anno,
      left_anno = left_anno,
      bottom_anno = bottom_anno,
      bar_value = bar_value,
      add_text = add_text,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      x_text_angle = x_text_angle,
      border = border,
      width = width,
      height = height,
      units = units,
      title = title,
      subtitle = subtitle,
      value_palette = palette_cfg$value_palette,
      value_palcolor = palette_cfg$value_palcolor,
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }

  if (identical(method, "CellChat")) {
    df <- extract_long_table(
      srt = srt,
      condition = condition,
      dataset = dataset,
      slot.name = slot.name,
      signaling = signaling,
      pairLR.use = pairLR.use,
      sources.use = sender.use,
      targets.use = receiver.use,
      thresh = thresh
    )
  } else {
    bundle <- get_bundle(srt, method = method)
    df <- bundle$long_table %||% data.frame()
  }

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

  score_col <- value
  if (!score_col %in% colnames(df)) {
    score_col <- c("score", "prob", "means")[
      c("score", "prob", "means") %in% colnames(df)
    ][1]
  }
  if (is.null(score_col) || is.na(score_col)) {
    score_col <- "score"
    if (!"score" %in% colnames(df)) df$score <- NA_real_
  }
  if (!"pvalue" %in% colnames(df)) {
    df$pvalue <- if ("pval" %in% colnames(df)) df$pval else NA_real_
  }
  df$score <- df[[score_col]]
  df <- prepare_plot_df(df)

  pair_df <- pair_plot_df(df)
  interaction_df <- interaction_plot_df(df)

  add_dot_eff <- identical(plot_type, "dot")
  value_use <- value
  if (!isTRUE(value_supplied) && isTRUE(edge_value_supplied)) {
    value_use <- edge_value
  }

  ccc_heatmap_full_plot(
    df = df,
    pair_df = pair_df,
    interaction_df = interaction_df,
    display_by = display_by,
    top_n = top_n,
    add_dot = add_dot_eff,
    top_anno = top_anno,
    right_anno = right_anno,
    left_anno = left_anno,
    bottom_anno = bottom_anno,
    value = value_use,
    bar_value = bar_value,
    add_text = add_text,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    edge_value = edge_value,
    color.by = color.by,
    x_text_angle = x_text_angle,
    facet_by = facet_by,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    title = title,
    subtitle = subtitle,
    width = width,
    height = height,
    units = units,
    border = border,
    value_palette = palette_cfg$value_palette,
    value_palcolor = palette_cfg$value_palcolor,
    cell_palette = palette_cfg$cell_palette,
    cell_palcolor = palette_cfg$cell_palcolor,
    legend.position = legend.position,
    legend.direction = legend.direction,
    font.size = font.size,
    theme_use = theme_use,
    theme_args = theme_args
  )
}

# Internal: full ComplexHeatmap rendering
ccc_heatmap_full_plot <- function(
  df,
  pair_df,
  interaction_df,
  display_by = "aggregation",
  top_n = 20,
  add_dot = FALSE,
  top_anno = "bar",
  right_anno = "cell",
  left_anno = "bar",
  bottom_anno = "cell",
  value = "sum",
  bar_value = "sum",
  add_text = NULL,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  edge_value = "sum",
  color.by = "score",
  x_text_angle = 90,
  facet_by = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  title = NULL,
  subtitle = NULL,
  width = NULL,
  height = NULL,
  units = "inch",
  border = TRUE,
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
  check_r("ComplexHeatmap", verbose = FALSE)
  check_r("circlize", verbose = FALSE)
  bar_value <- ccc_match_bar_value(bar_value)
  top_anno <- ccc_match_side_anno(top_anno)
  right_anno <- ccc_match_side_anno(right_anno)
  left_anno <- ccc_match_side_anno(left_anno)
  bottom_anno <- ccc_match_side_anno(bottom_anno)

  if (identical(display_by, "interaction")) {
    plot_value <- ccc_heatmap_value_spec(
      df = interaction_df,
      color.by = color.by,
      value = value
    )
    plot_df <- top_interactions(
      interaction_df,
      top_n = top_n,
      value_col = plot_value$var
    )
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No interaction-level CCC records are available for heatmap plotting",
        message_type = "error"
      )
    }

    # Build matrix: rows = pair or receiver/sender, cols = interaction_label
    if (is.null(facet_by)) {
      row_var <- "pair"
      row_title_str <- "Pair"
    } else if (identical(facet_by, "sender")) {
      row_var <- "receiver"
      row_title_str <- "Receiver"
    } else {
      row_var <- "sender"
      row_title_str <- "Sender"
    }

    score_var <- plot_value$var
    val_label <- plot_value$label

    mat <- ccc_pivot_matrix(
      df = plot_df,
      row_var = row_var,
      col_var = "interaction_label",
      val_var = score_var
    )

    add_text_eff <- if (is.null(add_text)) FALSE else isTRUE(add_text)

    # Bar data: col = interactions, row = pair/receiver/sender
    if (ccc_side_has_bar(top_anno, bottom_anno)) {
      top_bar_vec <- ccc_collect_bar_stats(
        values = plot_df[[score_var]],
        groups = plot_df$interaction_label,
        ordered_levels = colnames(mat),
        metrics = bar_value
      )
    } else {
      top_bar_vec <- NULL
    }
    if (ccc_side_has_distribution(top_anno, bottom_anno)) {
      top_box_vec <- ccc_collect_group_values(
        values = plot_df[[score_var]],
        groups = plot_df$interaction_label,
        ordered_levels = colnames(mat)
      )
    } else {
      top_box_vec <- NULL
    }
    if (ccc_side_has_summary(top_anno, bottom_anno)) {
      top_summary_vec <- ccc_collect_group_summary(
        values = plot_df[[score_var]],
        groups = plot_df$interaction_label,
        ordered_levels = colnames(mat),
        metric = "mean"
      )
    } else {
      top_summary_vec <- NULL
    }

    if (ccc_side_has_bar(left_anno, right_anno)) {
      right_bar_vec <- ccc_collect_bar_stats(
        values = plot_df[[score_var]],
        groups = plot_df[[row_var]],
        ordered_levels = rownames(mat),
        metrics = bar_value
      )
    } else {
      right_bar_vec <- NULL
    }
    if (ccc_side_has_distribution(left_anno, right_anno)) {
      right_box_vec <- ccc_collect_group_values(
        values = plot_df[[score_var]],
        groups = plot_df[[row_var]],
        ordered_levels = rownames(mat)
      )
    } else {
      right_box_vec <- NULL
    }
    if (ccc_side_has_summary(left_anno, right_anno)) {
      right_summary_vec <- ccc_collect_group_summary(
        values = plot_df[[score_var]],
        groups = plot_df[[row_var]],
        ordered_levels = rownames(mat),
        metric = "mean"
      )
    } else {
      right_summary_vec <- NULL
    }

    col_title_str <- "Interaction"
  } else {
    # aggregation: sender × receiver
    plot_value <- ccc_heatmap_value_spec(
      df = pair_df,
      color.by = color.by,
      value = value %||% edge_value
    )
    plot_df <- top_pairs(
      pair_df,
      top_n = top_n,
      value_col = plot_value$var
    )
    if (is.null(plot_df) || nrow(plot_df) == 0L) {
      log_message(
        "No aggregated sender-receiver interactions are available for heatmap plotting",
        message_type = "error"
      )
    }

    score_var <- plot_value$var
    val_label <- plot_value$label

    mat <- ccc_pivot_matrix(
      df = plot_df,
      row_var = "receiver",
      col_var = "sender",
      val_var = score_var
    )

    add_text_eff <- if (is.null(add_text)) TRUE else isTRUE(add_text)

    row_title_str <- "Receiver"
    col_title_str <- "Sender"

    if (ccc_side_has_bar(top_anno, bottom_anno)) {
      top_bar_vec <- ccc_collect_bar_stats(
        values = plot_df[[score_var]],
        groups = plot_df$sender,
        ordered_levels = colnames(mat),
        metrics = bar_value
      )
    } else {
      top_bar_vec <- NULL
    }
    if (ccc_side_has_distribution(top_anno, bottom_anno)) {
      top_box_vec <- ccc_collect_group_values(
        values = plot_df[[score_var]],
        groups = plot_df$sender,
        ordered_levels = colnames(mat)
      )
    } else {
      top_box_vec <- NULL
    }
    if (ccc_side_has_summary(top_anno, bottom_anno)) {
      top_summary_vec <- ccc_collect_group_summary(
        values = plot_df[[score_var]],
        groups = plot_df$sender,
        ordered_levels = colnames(mat),
        metric = "mean"
      )
    } else {
      top_summary_vec <- NULL
    }

    if (ccc_side_has_bar(left_anno, right_anno)) {
      right_bar_vec <- ccc_collect_bar_stats(
        values = plot_df[[score_var]],
        groups = plot_df$receiver,
        ordered_levels = rownames(mat),
        metrics = bar_value
      )
    } else {
      right_bar_vec <- NULL
    }
    if (ccc_side_has_distribution(left_anno, right_anno)) {
      right_box_vec <- ccc_collect_group_values(
        values = plot_df[[score_var]],
        groups = plot_df$receiver,
        ordered_levels = rownames(mat)
      )
    } else {
      right_box_vec <- NULL
    }
    if (ccc_side_has_summary(left_anno, right_anno)) {
      right_summary_vec <- ccc_collect_group_summary(
        values = plot_df[[score_var]],
        groups = plot_df$receiver,
        ordered_levels = rownames(mat),
        metric = "mean"
      )
    } else {
      right_summary_vec <- NULL
    }
  }

  # Color scale
  mat_vals <- mat[is.finite(mat)]
  if (length(mat_vals) == 0L) {
    mat_vals <- c(0, 1)
  }
  val_range <- range(mat_vals, na.rm = TRUE)
  if (val_range[1] == val_range[2]) {
    val_range <- val_range + c(-0.5, 0.5)
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

  # Cell text
  cell_fun <- if (isTRUE(add_text_eff) && !isTRUE(add_dot)) {
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

  # Dot layer_fun (hollow circles scaled by score)
  if (isTRUE(add_dot)) {
    mat_local <- mat
    mat_finite <- mat_local[is.finite(mat_local)]
    dot_max <- if (length(mat_finite) > 0) max(mat_finite) else 1
    if (dot_max == 0) {
      dot_max <- 1
    }
    layer_fun <- function(j, i, x, y, width, height, fill) {
      v <- ComplexHeatmap::pindex(mat_local, i, j)
      sz <- grid::unit(
        ifelse(is.finite(v), pmax(v / dot_max, 0.05), 0) * 8,
        "mm"
      )
      grid::grid.points(
        x,
        y,
        pch = 21,
        size = sz,
        gp = grid::gpar(
          fill = fill,
          col = if (isTRUE(border)) "grey20" else NA,
          lwd = 0.5
        )
      )
    }
    # suppress tile fill when in dot mode
    rect_gp <- grid::gpar(
      col = if (isTRUE(border)) "grey85" else NA,
      fill = NA,
      lwd = 0.6
    )
  } else {
    layer_fun <- NULL
    rect_gp <- grid::gpar(
      col = if (isTRUE(border)) "white" else NA,
      lwd = 1
    )
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
      df = data.frame(Sender = col_nms, row.names = col_nms),
      col = list(Sender = stats::setNames(col_cols, col_nms)),
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
      df = data.frame(Receiver = row_nms, row.names = row_nms),
      col = list(Receiver = stats::setNames(row_cols, row_nms)),
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
    palette_colors(
      top_names,
      palette = cell_palette,
      palcolor = cell_palcolor
    )
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

  ht_args <- list(
    matrix = mat,
    col = col_fun,
    name = val_label,
    cell_fun = cell_fun,
    layer_fun = layer_fun,
    top_annotation = top_combined,
    bottom_annotation = bottom_combined,
    left_annotation = left_combined,
    right_annotation = right_combined,
    row_title = row_title_str,
    column_title = col_title_str,
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
    rect_gp = rect_gp,
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

  ht <- do.call(ComplexHeatmap::Heatmap, ht_args)

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
    p <- thisplot::panel_fix_overall(
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
      do.call(theme_use, theme_args) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = font.size * 1.2, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = font.size, hjust = 0.5)
      )
  }

  p
}

ccc_ligand_target_heatmap <- function(
  df,
  top_n = 20,
  sender.use = NULL,
  receiver.use = NULL,
  context_df = NULL,
  sender_default = NULL,
  receiver_default = NULL,
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
  check_r("ComplexHeatmap", verbose = FALSE)
  check_r("circlize", verbose = FALSE)
  bar_value <- ccc_match_bar_value(bar_value)
  top_anno <- ccc_match_side_anno(top_anno)
  right_anno <- ccc_match_side_anno(right_anno)
  left_anno <- ccc_match_side_anno(left_anno)
  bottom_anno <- ccc_match_side_anno(bottom_anno)

  plot_df <- ccc_standardize_ligand_target_df(
    df = df,
    top_n = top_n,
    sender.use = sender.use,
    receiver.use = receiver.use,
    context_df = context_df,
    sender_default = sender_default,
    receiver_default = receiver_default
  )
  plot_df <- plot_df[is.finite(plot_df$weight), , drop = FALSE]
  if (nrow(plot_df) == 0L) {
    log_message(
      "No finite ligand-target weights are available for heatmap plotting",
      message_type = "error"
    )
  }

  ligand_levels <- levels(plot_df$ligand) %||%
    unique(as.character(plot_df$ligand))
  target_levels <- levels(plot_df$target) %||%
    unique(as.character(plot_df$target))
  mat <- matrix(
    NA_real_,
    nrow = length(ligand_levels),
    ncol = length(target_levels),
    dimnames = list(ligand_levels, target_levels)
  )
  for (k in seq_len(nrow(plot_df))) {
    ligand_k <- as.character(plot_df$ligand[k])
    target_k <- as.character(plot_df$target[k])
    weight_k <- as.numeric(plot_df$weight[k])
    if (is.na(ligand_k) || is.na(target_k) || !is.finite(weight_k)) {
      next
    }
    existing <- mat[ligand_k, target_k]
    mat[ligand_k, target_k] <- if (is.na(existing)) {
      weight_k
    } else {
      existing + weight_k
    }
  }

  add_text_eff <- if (is.null(add_text)) TRUE else isTRUE(add_text)
  top_bar_vec <- if (ccc_side_has_bar(top_anno, bottom_anno)) {
    ccc_collect_bar_stats(
      values = plot_df$weight,
      groups = plot_df$target,
      ordered_levels = colnames(mat),
      metrics = bar_value
    )
  } else {
    NULL
  }
  top_box_vec <- if (ccc_side_has_distribution(top_anno, bottom_anno)) {
    ccc_collect_group_values(
      values = plot_df$weight,
      groups = plot_df$target,
      ordered_levels = colnames(mat)
    )
  } else {
    NULL
  }
  top_summary_vec <- if (ccc_side_has_summary(top_anno, bottom_anno)) {
    ccc_collect_group_summary(
      values = plot_df$weight,
      groups = plot_df$target,
      ordered_levels = colnames(mat),
      metric = "mean"
    )
  } else {
    NULL
  }
  right_bar_vec <- if (ccc_side_has_bar(left_anno, right_anno)) {
    ccc_collect_bar_stats(
      values = plot_df$weight,
      groups = plot_df$ligand,
      ordered_levels = rownames(mat),
      metrics = bar_value
    )
  } else {
    NULL
  }
  right_box_vec <- if (ccc_side_has_distribution(left_anno, right_anno)) {
    ccc_collect_group_values(
      values = plot_df$weight,
      groups = plot_df$ligand,
      ordered_levels = rownames(mat)
    )
  } else {
    NULL
  }
  right_summary_vec <- if (ccc_side_has_summary(left_anno, right_anno)) {
    ccc_collect_group_summary(
      values = plot_df$weight,
      groups = plot_df$ligand,
      ordered_levels = rownames(mat),
      metric = "mean"
    )
  } else {
    NULL
  }

  mat_vals <- mat[is.finite(mat)]
  if (length(mat_vals) == 0L) {
    mat_vals <- c(0, 1)
  }
  val_range <- range(mat_vals, na.rm = TRUE)
  if (val_range[1] == val_range[2]) {
    val_range <- val_range + c(-0.5, 0.5)
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

  cell_fun <- if (isTRUE(add_text_eff)) {
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
      df = data.frame(Target = col_nms, row.names = col_nms),
      col = list(Target = stats::setNames(col_cols, col_nms)),
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
      df = data.frame(Ligand = row_nms, row.names = row_nms),
      col = list(Ligand = stats::setNames(row_cols, row_nms)),
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
    name = "Weight",
    cell_fun = cell_fun,
    top_annotation = top_combined,
    bottom_annotation = bottom_combined,
    left_annotation = left_combined,
    right_annotation = right_combined,
    row_title = "Ligand",
    column_title = "Target",
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
    p <- thisplot::panel_fix_overall(
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
      do.call(theme_use, theme_args) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = font.size * 1.2, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = font.size, hjust = 0.5)
      )
  }

  p
}

ccc_matrix_group_values_from_matrix <- function(mat, margin = c("col", "row")) {
  margin <- match.arg(margin)
  if (is.null(mat) || length(mat) == 0L) {
    return(NULL)
  }
  idx <- if (identical(margin, "col")) {
    seq_len(ncol(mat))
  } else {
    seq_len(nrow(mat))
  }
  values <- lapply(idx, function(i) {
    v <- if (identical(margin, "col")) mat[, i] else mat[i, ]
    v[is.finite(v)]
  })
  names(values) <- if (identical(margin, "col")) {
    colnames(mat)
  } else {
    rownames(mat)
  }
  values
}

ccc_matrix_bar_stats_from_values <- function(values_list, metrics = "sum") {
  if (is.null(values_list) || length(values_list) == 0L) {
    return(NULL)
  }
  metrics <- ccc_match_bar_value(metrics)
  out <- lapply(metrics, function(metric) {
    stats::setNames(
      vapply(
        values_list,
        function(v) {
          if (length(v) == 0L) {
            return(0)
          }
          switch(
            metric,
            count = sum(is.finite(v)),
            sum = sum(v, na.rm = TRUE),
            mean = mean(v, na.rm = TRUE),
            max = max(v, na.rm = TRUE)
          )
        },
        numeric(1)
      ),
      names(values_list)
    )
  })
  names(out) <- metrics
  out
}

ccc_matrix_summary_from_values <- function(values_list, metric = "mean") {
  if (is.null(values_list) || length(values_list) == 0L) {
    return(NULL)
  }
  stats::setNames(
    vapply(
      values_list,
      function(v) {
        if (length(v) == 0L) {
          return(NA_real_)
        }
        switch(
          metric,
          mean = mean(v, na.rm = TRUE),
          max = max(v, na.rm = TRUE),
          sum = sum(v, na.rm = TRUE)
        )
      },
      numeric(1)
    ),
    names(values_list)
  )
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
  check_r("ComplexHeatmap", verbose = FALSE)
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
    p <- thisplot::panel_fix_overall(
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
      do.call(theme_use, theme_args) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = font.size * 1.2, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = font.size, hjust = 0.5)
      )
  }

  p
}

ccc_cellchat_role_matrix <- function(
  object,
  signaling = NULL,
  pattern = c("outgoing", "incoming", "all"),
  scale_rows = TRUE
) {
  pattern <- match.arg(pattern)
  if (length(object@netP$centr) == 0L) {
    log_message(
      "Please run CellChat centrality computation before plotting role heatmaps",
      message_type = "error"
    )
  }
  centr <- object@netP$centr
  groups <- levels(object@idents)
  outgoing <- matrix(0, nrow = length(centr), ncol = length(groups))
  incoming <- matrix(0, nrow = length(centr), ncol = length(groups))
  dimnames(outgoing) <- list(names(centr), groups)
  dimnames(incoming) <- dimnames(outgoing)
  for (i in seq_along(centr)) {
    outgoing[i, ] <- centr[[i]]$outdeg
    incoming[i, ] <- centr[[i]]$indeg
  }
  mat <- switch(
    pattern,
    outgoing = outgoing,
    incoming = incoming,
    all = outgoing + incoming
  )
  if (!is.null(signaling)) {
    signaling <- unique(as.character(signaling))
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(
      0,
      nrow = length(signaling),
      ncol = ncol(mat),
      dimnames = list(signaling, colnames(mat))
    )
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
  }
  raw_mat <- mat
  if (isTRUE(scale_rows)) {
    row_max <- apply(mat, 1, max, na.rm = TRUE)
    row_max[!is.finite(row_max) | row_max == 0] <- 1
    mat <- sweep(mat, 1L, row_max, "/", check.margin = FALSE)
    mat[mat == 0] <- NA_real_
  }
  list(raw = raw_mat, scaled = mat)
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
  build_role_heatmap <- function(mat, plot_title, plot_subtitle = NULL, symmetric = FALSE) {
    ccc_matrix_heatmap_plot(
      mat = mat,
      value_label = if (isTRUE(symmetric)) "Difference" else "Relative strength",
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

ccc_cellchat_heatmap_comparison <- function(srt, condition = NULL, comparison = c(1, 2), min_n = 1L) {
  cmp <- .cc_get_cmp(srt = srt, condition = condition)
  comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)
  if (length(comp_idx) < min_n) {
    log_message(
      paste0(
        "{.arg comparison} must contain at least ", min_n,
        " datasets for CellChat comparison heatmap plotting"
      ),
      message_type = "error"
    )
  }
  object_names <- names(cmp$object.list)[comp_idx]
  object_list <- cmp$object.list[object_names]
  list(
    cmp = cmp,
    comp_idx = comp_idx,
    object_names = object_names,
    object_list = object_list
  )
}

ccc_heatmap_value_spec <- function(df, color.by = "score", value = "sum") {
  value <- value %||% "sum"
  spec <- scale_var(df = df, color.by = color.by, agg_value = value)
  spec$label <- if (identical(spec$var, "specificity")) {
    "-log10(p)"
  } else {
    spec$label
  }
  spec
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

ccc_side_has_box <- function(...) {
  annos <- unlist(list(...), use.names = FALSE)
  "box" %in% annos
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
  anno <- switch(
    type,
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
  anno <- switch(
    type,
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

ccc_pivot_matrix <- function(df, row_var, col_var, val_var) {
  rows <- sort(unique(as.character(df[[row_var]])))
  cols <- sort(unique(as.character(df[[col_var]])))
  mat <- matrix(
    NA_real_,
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(rows, cols)
  )
  for (k in seq_len(nrow(df))) {
    r <- as.character(df[[row_var]][k])
    c <- as.character(df[[col_var]][k])
    v <- df[[val_var]][k]
    if (!is.na(r) && !is.na(c) && r %in% rows && c %in% cols) {
      existing <- mat[r, c]
      mat[r, c] <- if (is.na(existing)) {
        as.numeric(v)
      } else {
        existing + as.numeric(v)
      }
    }
  }
  mat
}
