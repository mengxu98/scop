#' @title CCC network and flow plots
#'
#' @md
#' @inheritParams CCCStatPlot
#' @param plot_type Plot type. One of `"circle"`, `"chord"`, `"pathway"`,
#'   `"individual_lr"`, `"arrow"`, `"sigmoid"`, `"bipartite"`,
#'   `"embedding_network"`, or `"diff_network"`.
#' @param ligand For `plot_type = "bipartite"`: the ligand name to focus on.
#'   If `NULL`, the ligand with the highest total score is used.
#' @param receptor For `plot_type = "bipartite"`: optional receptor names to
#'   restrict to. If `NULL`, all receptors paired with `ligand` are shown.
#' @param reg.by For `plot_type = "bipartite"`: optional metadata column in
#'   `srt` used to color edges by regulation status (e.g. up/down). If `NULL`,
#'   edges are colored by sender cell type.
#' @param reg_palette For `plot_type = "bipartite"`: named character vector or
#'   palette name for regulation categories.
#' @param reg_palcolor For `plot_type = "bipartite"`: custom colors for
#'   regulation palette.
#' @param expr.by For `plot_type = "bipartite"`: optional metadata or score
#'   column used to scale edge line width. If `NULL`, all edges have equal
#'   width.
#' @param group.by For `plot_type = "embedding_network"`: metadata column used
#'   to define cell groups. If `NULL`, the grouping stored in the CCC result is
#'   used when available.
#' @param reduction For `plot_type = "embedding_network"`: dimensional reduction
#'   to use. If `NULL`, the default reduction is used.
#' @param dims For `plot_type = "embedding_network"`: dimensions to plot.
#' @param layout Layout used for graph-based network views. `"chord"` can also
#'   be requested via `plot_type = "circle"` for backward compatibility.
#' @param link_curvature Curvature used for circle-like differential links and
#'   flow edges.
#' @param edge_size Range used for scaling edge widths.
#' @param edge_color Optional edge color override. For differential networks,
#'   this may also be a length-2 vector for negative/positive changes.
#' @param edge_alpha Alpha used for embedding-network edges.
#' @param edge_line Edge geometry for `plot_type = "arrow"`, `"sigmoid"`, and
#'   `"embedding_network"`.
#' @param edge_curvature Curvature used for curved flow/embedding edges.
#' @param directed Whether to draw arrows for directed networks.
#' @param arrow_type Arrow head type passed to `grid::arrow()`.
#' @param arrow_angle Arrow head angle passed to `grid::arrow()`.
#' @param arrow_length Arrow length passed to `grid::arrow()`.
#' @param node_size Base node size.
#' @param node_alpha Node alpha.
#' @param legend.title Legend title.
#' @param ... Additional plot-specific options. For chord plots, `reduce`,
#'   `max.groups`, `small.gap`, `big.gap`, and `lab.cex` can be used to adjust
#'   the CellChat-like chord layout.
#'
#' @return A ggplot, patchwork, or recorded plot object.
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
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "circle",
#'   display_by = "aggregation",
#'   value = "count",
#'   top_n = 20
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "circle",
#'   display_by = "aggregation",
#'   value = "weight",
#'   top_n = 20
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "chord",
#'   display_by = "aggregation",
#'   top_n = 12
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "arrow",
#'   display_by = "interaction",
#'   sender.use = "Ductal",
#'   receiver.use = "Ngn3-low-EP",
#'   edge_line = "straight",
#'   directed = TRUE,
#'   top_n = 3
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "arrow",
#'   display_by = "interaction",
#'   top_n = 20
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "sigmoid",
#'   display_by = "interaction",
#'   top_n = 20
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "bipartite",
#'   display_by = "aggregation",
#'   top_n = 20
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "embedding_network",
#'   group.by = "CellType",
#'   reduction = "UMAP",
#'   top_n = 20,
#'   label = TRUE
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "pathway",
#'   signaling = "MK"
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA",
#'   plot_type = "individual_lr",
#'   signaling = "MK",
#'   pairLR.use = "MDK_SDC1"
#' )
#'
#' CCCNetworkPlot(
#'   pancreas_sub,
#'   method = "CellChat",
#'   condition = "ConditionA_vs_ConditionB",
#'   plot_type = "diff_network",
#'   measure = "count"
#' )
CCCNetworkPlot <- function(
  srt,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c(
    "circle",
    "chord",
    "pathway",
    "individual_lr",
    "arrow",
    "sigmoid",
    "bipartite",
    "embedding_network",
    "diff_network"
  ),
  display_by = c("aggregation", "interaction"),
  sender.use = NULL,
  receiver.use = NULL,
  ligand.use = NULL,
  receptor.use = NULL,
  interaction.use = NULL,
  group.by = NULL,
  reduction = NULL,
  dims = c(1, 2),
  signaling = NULL,
  pairLR.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  measure = c("count", "weight"),
  value = "score",
  top_n = 20,
  ligand = NULL,
  receptor = NULL,
  reg.by = NULL,
  reg_palette = "Set1",
  reg_palcolor = NULL,
  expr.by = NULL,
  layout = c(
    "circle",
    "hierarchy",
    "chord",
    "kk",
    "fr",
    "nicely",
    "lgl",
    "mds",
    "graphopt"
  ),
  link_curvature = 0.2,
  link_alpha = 0.6,
  edge_value = c("sum", "mean", "max", "count"),
  edge_threshold = 0,
  edge_size = c(0.5, 1.8),
  edge_color = NULL,
  edge_alpha = 0.6,
  edge_line = c("curved", "straight"),
  edge_curvature = 0.2,
  directed = FALSE,
  arrow_type = "closed",
  arrow_angle = 20,
  arrow_length = grid::unit(0.02, "npc"),
  node_size = 5,
  node_alpha = 0.9,
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
  legend.title = NULL,
  font.size = 10,
  theme_use = "theme_scop",
  theme_args = list(),
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
  layout <- match.arg(layout)
  edge_value <- match.arg(edge_value)
  edge_line <- match.arg(edge_line)
  if (identical(plot_type, "circle") && identical(layout, "chord")) {
    plot_type <- "chord"
  }
  dots <- list(...)
  label.enable <- isTRUE(dots[["label"]])
  label.size <- dots[["label.size"]] %||% 4
  label.fg <- dots[["label.fg"]] %||% "white"
  label.bg <- dots[["label.bg"]] %||% "black"
  label.bg.r <- dots[["label.bg.r"]] %||% 0.1
  dots[c("label", "label.size", "label.fg", "label.bg", "label.bg.r")] <- NULL
  reduce <- if (is.null(dots[["reduce"]])) TRUE else isTRUE(dots[["reduce"]])
  max.groups <- dots[["max.groups"]] %||% 8
  small.gap <- dots[["small.gap"]] %||% 1
  big.gap <- dots[["big.gap"]] %||% 8
  lab.cex <- dots[["lab.cex"]] %||% 0.6
  palette_cfg <- ccc_palettes(
    palette = palette,
    palcolor = palcolor,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    link_palette = link_palette,
    link_palcolor = link_palcolor
  )

  method <- detect_method(srt = srt, method = method)

  if (plot_type %in% c("pathway", "individual_lr")) {
    if (identical(method, "CellChat")) {
      if (!identical(layout, "circle")) {
        log_message(
          "{.arg layout} is ignored for {.val plot_type = {plot_type}}; scop-native CellChat pathway/LR networks use the circle layout",
          message_type = "warning"
        )
      }
      obj_info <- get_dataset_object(
        srt = srt,
        condition = condition,
        dataset = dataset
      )
      cellchat_object <- obj_info$object
      plot_cellchat_circle <- function(
        sig,
        pairLR = NULL,
        plot_title = NULL,
        plot_subtitle = NULL
      ) {
        ccc_cellchat_circle_network_plot(
          srt = srt,
          condition = condition,
          dataset = dataset,
          signaling = sig,
          pairLR.use = pairLR,
          sender.use = sender.use,
          receiver.use = receiver.use,
          slot.name = slot.name,
          thresh = thresh,
          display_by = display_by,
          value = value,
          top_n = top_n,
          edge_threshold = edge_threshold,
          edge_size = edge_size,
          node_size = node_size,
          node_alpha = node_alpha,
          link_alpha = link_alpha,
          cell_palette = palette_cfg$cell_palette,
          cell_palcolor = palette_cfg$cell_palcolor,
          link_palette = palette_cfg$link_palette,
          link_palcolor = palette_cfg$link_palcolor,
          title = plot_title,
          subtitle = plot_subtitle,
          legend.position = legend.position,
          legend.direction = legend.direction,
          legend.title = legend.title,
          font.size = font.size,
          theme_use = theme_use,
          theme_args = theme_args
        )
      }

      if (identical(plot_type, "pathway")) {
        pathways_to_show <- signaling %||%
          cellchat_object@netP$pathways[seq_len(min(
            top_n,
            length(cellchat_object@netP$pathways)
          ))]
        pathways_to_show <- pathways_to_show[!is.na(pathways_to_show)]
        if (length(pathways_to_show) == 0L) {
          log_message(
            "No signaling pathways available for plotting",
            message_type = "error"
          )
        }
        plots <- lapply(pathways_to_show, function(sig) {
          plot_cellchat_circle(
            sig = sig,
            plot_title = if (length(pathways_to_show) == 1L) {
              title %||% sig
            } else {
              sig
            },
            plot_subtitle = if (length(pathways_to_show) == 1L) subtitle else NULL
          )
        })
        return(simplify_cc_plot_list(plots))
      }

      if (is.null(signaling)) {
        log_message(
          "{.arg signaling} must be provided for {.val plot_type = 'individual_lr'}",
          message_type = "error"
        )
      }
      if (is.null(pairLR.use)) {
        log_message(
          "{.arg pairLR.use} must be provided for {.val plot_type = 'individual_lr'}",
          message_type = "error"
        )
      }
      plots <- lapply(signaling, function(sig) {
        plot_cellchat_circle(
          sig = sig,
          pairLR = pairLR.use,
          plot_title = title %||%
            paste(
              sig,
              paste(as.character(pairLR.use), collapse = ", "),
              sep = ": "
            ),
          plot_subtitle = subtitle
        )
      })
      return(simplify_cc_plot_list(plots))
    }

    if (!identical(layout, "circle")) {
      log_message(
        "{.arg layout} is ignored for {.val plot_type = {plot_type}} when using pathway-aware generic CCC results; the circle layout is used",
        message_type = "warning"
      )
    }

    bundle <- get_bundle(srt, method = method)
    long_df <- standardize_long_df(bundle$long_table %||% data.frame())
    available_pathways <- unique(as.character(long_df$pathway_name))
    available_pathways <- available_pathways[
      !is.na(available_pathways) & nzchar(available_pathways)
    ]
    if (length(available_pathways) == 0L) {
      log_message(
        paste0(
          "{.arg plot_type} = {.val {plot_type}} requires pathway annotations ",
          "(for example {.code pathway_name} / {.code classification}) in the ",
          "stored CCC result table"
        ),
        message_type = "error"
      )
    }

    if (identical(plot_type, "individual_lr")) {
      if (is.null(signaling)) {
        log_message(
          "{.arg signaling} must be provided for {.val plot_type = 'individual_lr'}",
          message_type = "error"
        )
      }
      if (is.null(pairLR.use)) {
        log_message(
          "{.arg pairLR.use} must be provided for {.val plot_type = 'individual_lr'}",
          message_type = "error"
        )
      }
      pathways_to_show <- unique(as.character(signaling))
    } else {
      pathways_to_show <- unique(as.character(signaling %||% available_pathways))
      if (
        is.numeric(top_n) &&
          length(top_n) == 1L &&
          top_n > 0L &&
          length(pathways_to_show) > top_n
      ) {
        pathways_to_show <- utils::head(pathways_to_show, top_n)
      }
    }

    plot_generic_circle <- function(sig) {
      df_sig <- filter_long_df(
        df = long_df,
        sender.use = sender.use,
        receiver.use = receiver.use,
        ligand.use = ligand.use,
        receptor.use = receptor.use,
        interaction.use = interaction.use,
        signaling = sig,
        pairLR.use = if (identical(plot_type, "individual_lr")) pairLR.use else NULL
      )
      if (nrow(df_sig) == 0L) {
        return(NULL)
      }
      df_sig <- ccc_assign_plot_score(df = df_sig, value = value)
      df_sig <- prepare_plot_df(df_sig)
      pair_df_sig <- pair_plot_df(df_sig)
      interaction_df_sig <- interaction_plot_df(df_sig)
      do.call(
        ccc_circle_plot,
        c(
          list(
            pair_df = pair_df_sig,
            interaction_df = interaction_df_sig,
            display_by = display_by,
            top_n = top_n,
            value = value,
            edge_threshold = edge_threshold,
            edge_size = edge_size,
            node_size = node_size,
            node_alpha = node_alpha,
            link_alpha = link_alpha
          ),
          list(
            title = if (identical(plot_type, "individual_lr")) {
              title %||% paste(
                sig,
                paste(as.character(pairLR.use), collapse = ", "),
                sep = ": "
              )
            } else if (length(pathways_to_show) == 1L) {
              title %||% sig
            } else {
              sig
            },
            subtitle = if (length(pathways_to_show) == 1L) subtitle else NULL,
            cell_palette = palette_cfg$cell_palette,
            cell_palcolor = palette_cfg$cell_palcolor,
            link_palette = palette_cfg$link_palette,
            link_palcolor = palette_cfg$link_palcolor,
            legend.position = legend.position,
            legend.direction = legend.direction,
            legend.title = legend.title,
            font.size = font.size,
            theme_use = theme_use,
            theme_args = theme_args,
            label = label.enable,
            label.size = label.size,
            label.fg = label.fg,
            label.bg = label.bg,
            label.bg.r = label.bg.r
          )
        )
      )
    }

    plots <- Filter(Negate(is.null), lapply(pathways_to_show, plot_generic_circle))
    if (length(plots) == 0L) {
      log_message(
        "No pathway-specific communication records remain after filtering",
        message_type = "error"
      )
    }
    return(simplify_cc_plot_list(plots))
  }

  if (identical(plot_type, "diff_network")) {
    if (!identical(method, "CellChat")) {
      log_message(
        "{.val plot_type = 'diff_network'} is currently only supported for {.pkg CellChat}",
        message_type = "error"
      )
    }
    layout_use <- if (layout %in% c("hierarchy", "chord")) "circle" else layout
    return(ccc_diff_network_plot(
      srt = srt,
      condition = condition,
      comparison = comparison,
      measure = measure,
      sender.use = sender.use,
      receiver.use = receiver.use,
      top_n = top_n,
      edge_threshold = edge_threshold,
      edge_size = edge_size,
      edge_color = edge_color,
      link_curvature = link_curvature,
      link_alpha = link_alpha,
      directed = directed,
      arrow_type = arrow_type,
      arrow_angle = arrow_angle,
      arrow_length = arrow_length,
      node_size = node_size,
      node_alpha = node_alpha,
      layout = layout_use,
      title = title,
      subtitle = subtitle,
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

  df <- ccc_assign_plot_score(df = df, value = value)
  df <- prepare_plot_df(df)
  pair_df <- pair_plot_df(df)
  interaction_df <- interaction_plot_df(df)
  network_plot_args <- list(
    title = title,
    subtitle = subtitle,
    cell_palette = palette_cfg$cell_palette,
    cell_palcolor = palette_cfg$cell_palcolor,
    link_palette = palette_cfg$link_palette,
    link_palcolor = palette_cfg$link_palcolor,
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.title = legend.title,
    font.size = font.size,
    theme_use = theme_use,
    theme_args = theme_args
  )
  label_plot_args <- list(
    label = label.enable,
    label.size = label.size,
    label.fg = label.fg,
    label.bg = label.bg,
    label.bg.r = label.bg.r
  )

  if (identical(plot_type, "circle")) {
    return(do.call(
      ccc_circle_plot,
      c(
        list(
          pair_df = pair_df,
          interaction_df = interaction_df,
          display_by = display_by,
          top_n = top_n,
          value = value,
          edge_threshold = edge_threshold,
          edge_size = edge_size,
          node_size = node_size,
          node_alpha = node_alpha,
          link_alpha = link_alpha
        ),
        network_plot_args,
        label_plot_args
      )
    ))
  }

  if (identical(plot_type, "chord")) {
    dots_chord <- dots
    dots_chord[c("reduce", "max.groups", "small.gap", "big.gap", "lab.cex")] <- NULL
    return(do.call(
      ccc_chord_plot,
      c(
        list(
          pair_df = pair_df,
          interaction_df = interaction_df,
          display_by = display_by,
          top_n = top_n,
          edge_value = edge_value,
          edge_threshold = edge_threshold,
          link_alpha = link_alpha,
          reduce = reduce,
          max.groups = max.groups,
          small.gap = small.gap,
          big.gap = big.gap,
          lab.cex = lab.cex
        ),
        network_plot_args,
        dots_chord
      )
    ))
  }

  if (plot_type %in% c("arrow", "sigmoid")) {
    return(do.call(
      ccc_flow_network_plot,
      c(
        list(
          pair_df = pair_df,
          interaction_df = interaction_df,
          plot_type = plot_type,
          display_by = display_by,
          top_n = top_n,
          edge_value = edge_value,
          edge_threshold = edge_threshold,
          edge_size = edge_size,
          edge_color = edge_color,
          link_curvature = link_curvature,
          link_alpha = link_alpha,
          directed = directed,
          arrow_type = arrow_type,
          arrow_angle = arrow_angle,
          arrow_length = arrow_length,
          node_size = node_size,
          node_alpha = node_alpha
        ),
        network_plot_args,
        label_plot_args
      )
    ))
  }

  if (identical(plot_type, "bipartite")) {
    reg_vec <- NULL
    if (!is.null(reg.by)) {
      if (reg.by %in% colnames(srt@meta.data)) {
        reg_vec <- srt@meta.data[[reg.by]]
      } else if (reg.by %in% colnames(df)) {
        reg_vec <- NULL
      }
    }
    return(bipartite_plot(
      df = df,
      ligand = ligand,
      receptor = receptor,
      top_n = top_n,
      reg.by = if (!is.null(reg.by) && reg.by %in% colnames(df)) {
        reg.by
      } else {
        NULL
      },
      reg_palette = reg_palette,
      reg_palcolor = reg_palcolor,
      expr.by = if (!is.null(expr.by) && expr.by %in% colnames(df)) {
        expr.by
      } else {
        NULL
      },
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      node_size = node_size,
      node_alpha = node_alpha,
      link_alpha = link_alpha,
      edge_size = edge_size,
      title = title,
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    ))
  }

  if (identical(plot_type, "embedding_network")) {
    return(do.call(
      ccc_dim_network_plot,
      c(
        list(
          srt = srt,
          pair_df = pair_df,
          method = method,
          group.by = group.by,
          reduction = reduction,
          dims = dims,
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
        ),
        network_plot_args,
        label_plot_args,
        dots
      )
    ))
  }

  log_message(
    "Unsupported {.arg plot_type}: {.val {plot_type}}",
    message_type = "error"
  )
}

ccc_cellchat_circle_network_plot <- function(
  srt,
  condition = NULL,
  dataset = 1,
  signaling = NULL,
  pairLR.use = NULL,
  sender.use = NULL,
  receiver.use = NULL,
  slot.name = "net",
  thresh = 0.05,
  display_by = "aggregation",
  value = "weight",
  top_n = 20,
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
  theme_args = list()
) {
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
  df <- standardize_long_df(df)
  df <- filter_long_df(
    df = df,
    sender.use = sender.use,
    receiver.use = receiver.use,
    signaling = signaling,
    pairLR.use = pairLR.use
  )
  df <- ccc_assign_plot_score(df = df, value = value)
  df <- prepare_plot_df(df)
  pair_df <- pair_plot_df(df)
  interaction_df <- interaction_plot_df(df)
  ccc_circle_plot(
    pair_df = pair_df,
    interaction_df = interaction_df,
    display_by = display_by,
    top_n = top_n,
    value = value,
    edge_threshold = edge_threshold,
    edge_size = edge_size,
    node_size = node_size,
    node_alpha = node_alpha,
    link_alpha = link_alpha,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    link_palette = link_palette,
    link_palcolor = link_palcolor,
    title = title,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    font.size = font.size,
    theme_use = theme_use,
    theme_args = theme_args
  )
}

bipartite_plot <- function(
  df,
  ligand = NULL,
  receptor = NULL,
  top_n = 20,
  reg.by = NULL,
  reg_palette = "Set1",
  reg_palcolor = NULL,
  expr.by = NULL,
  cell_palette = "Chinese",
  cell_palcolor = NULL,
  node_size = 5,
  node_alpha = 0.9,
  link_alpha = 0.6,
  edge_size = c(0.5, 2.5),
  title = NULL,
  subtitle = NULL,
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
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No CCC records are available for {.pkg bipartite} plotting",
      message_type = "error"
    )
  }

  if (!"ligand" %in% colnames(df)) {
    df$ligand <- NA_character_
  }
  if (!"receptor" %in% colnames(df)) {
    df$receptor <- NA_character_
  }

  df <- df[!is.na(df$ligand) & nzchar(df$ligand), , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message(
      "No rows with non-missing ligand remain for {.pkg bipartite} plotting",
      message_type = "error"
    )
  }

  if (is.null(ligand)) {
    lig_scores <- tapply(df$score, df$ligand, sum, na.rm = TRUE)
    ligand <- names(which.max(lig_scores))
  }
  df <- df[as.character(df$ligand) == as.character(ligand), , drop = FALSE]

  if (!is.null(receptor)) {
    df <- df[
      as.character(df$receptor) %in% as.character(receptor),
      ,
      drop = FALSE
    ]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No CCC records remain for ligand {.val {ligand}} after filtering",
      message_type = "error"
    )
  }

  if (is.numeric(top_n) && top_n > 0L) {
    sender_scores <- tapply(df$score, df$sender, sum, na.rm = TRUE)
    receiver_scores <- tapply(df$score, df$receiver, sum, na.rm = TRUE)
    top_senders <- names(sort(sender_scores, decreasing = TRUE))[
      seq_len(min(top_n, length(sender_scores)))
    ]
    top_receivers <- names(sort(receiver_scores, decreasing = TRUE))[
      seq_len(min(top_n, length(receiver_scores)))
    ]
    df <- df[
      df$sender %in% top_senders & df$receiver %in% top_receivers,
      ,
      drop = FALSE
    ]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No CCC records remain after top_n filtering for {.pkg bipartite} plot",
      message_type = "error"
    )
  }

  senders <- unique(as.character(df$sender))
  receptors <- unique(as.character(df$receptor))
  receivers <- unique(as.character(df$receiver))
  pretty_label <- function(x, node_type = c("cell", "ligand", "receptor")) {
    node_type <- match.arg(node_type)
    x <- as.character(x)
    if (identical(node_type, "cell")) {
      return(x)
    }
    x <- gsub("_", "\n", x, fixed = TRUE)
    x
  }

  make_col_nodes <- function(labels, x_pos, col_type) {
    n <- length(labels)
    if (n == 0L) {
      return(data.frame())
    }
    y_pos <- if (n == 1L) {
      0.45
    } else {
      seq(0.75, 0.15, length.out = n)
    }
    data.frame(
      id = paste0(col_type, "::", labels),
      label = labels,
      label_plot = if (
        identical(col_type, "sender") || identical(col_type, "receiver")
      ) {
        pretty_label(labels, "cell")
      } else if (identical(col_type, "ligand")) {
        pretty_label(labels, "ligand")
      } else {
        pretty_label(labels, "receptor")
      },
      col_type = col_type,
      x = x_pos,
      y = y_pos,
      stringsAsFactors = FALSE
    )
  }

  sender_nodes <- make_col_nodes(rev(senders), 0, "sender")
  ligand_nodes <- data.frame(
    id = paste0("ligand::", ligand),
    label = ligand,
    label_plot = pretty_label(ligand, "ligand"),
    col_type = "ligand",
    x = 1,
    y = if (nrow(sender_nodes) > 0L) stats::median(sender_nodes$y) else 0.45,
    stringsAsFactors = FALSE
  )
  receptor_nodes <- make_col_nodes(rev(receptors), 2, "receptor")
  receiver_nodes <- make_col_nodes(rev(receivers), 3, "receiver")

  node_df <- rbind(sender_nodes, ligand_nodes, receptor_nodes, receiver_nodes)

  all_cell_types <- unique(c(senders, receivers))
  cell_cols <- palette_colors(
    all_cell_types,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  node_df$fill <- ifelse(
    node_df$col_type %in% c("sender", "receiver"),
    unname(cell_cols[node_df$label]),
    "white"
  )
  node_df$border <- ifelse(
    node_df$col_type %in% c("ligand", "receptor"),
    "grey50",
    unname(cell_cols[node_df$label])
  )
  node_df$shape <- ifelse(
    node_df$col_type %in% c("ligand", "receptor"),
    22,
    21
  )
  node_df$is_cell <- node_df$col_type %in% c("sender", "receiver")

  node_lookup <- stats::setNames(
    seq_len(nrow(node_df)),
    node_df$id
  )
  get_pos <- function(id_vec) {
    idx <- node_lookup[id_vec]
    list(x = node_df$x[idx], y = node_df$y[idx])
  }

  sender_edge_ids <- paste0("sender::", df$sender)
  ligand_edge_id <- paste0("ligand::", ligand)
  sl_pos_from <- get_pos(sender_edge_ids)
  sl_pos_to <- get_pos(rep(ligand_edge_id, nrow(df)))

  receptor_edge_ids <- paste0("receptor::", df$receptor)
  lr_pos_from <- get_pos(rep(ligand_edge_id, nrow(df)))
  lr_pos_to <- get_pos(receptor_edge_ids)

  receiver_edge_ids <- paste0("receiver::", df$receiver)
  rr_pos_from <- get_pos(receptor_edge_ids)
  rr_pos_to <- get_pos(receiver_edge_ids)

  rr_pos_from <- get_pos(receptor_edge_ids)
  rr_pos_to <- get_pos(receiver_edge_ids)

  weight_col <- expr.by %||% "score"
  if (!weight_col %in% colnames(df)) {
    weight_col <- "score"
  }
  weights <- as.numeric(df[[weight_col]])
  weights[!is.finite(weights)] <- 0

  if (!is.null(reg.by) && reg.by %in% colnames(df)) {
    reg_vals <- as.character(df[[reg.by]])
    reg_levels <- unique(reg_vals)
    reg_cols <- palette_colors(
      reg_levels,
      palette = reg_palette,
      palcolor = reg_palcolor,
      NA_keep = TRUE
    )
    edge_fill <- unname(reg_cols[reg_vals])
  } else {
    edge_fill <- unname(cell_cols[as.character(df$sender)])
  }
  edge_fill[is.na(edge_fill)] <- "grey60"

  w_range <- range(weights, na.rm = TRUE)
  if (diff(w_range) > 0) {
    lwd <- scales::rescale(weights, to = edge_size, from = w_range)
  } else {
    lwd <- rep(mean(edge_size), length(weights))
  }

  edge_df <- data.frame(
    x_from_sl = sl_pos_from$x,
    y_from_sl = sl_pos_from$y,
    x_to_sl = sl_pos_to$x,
    y_to_sl = sl_pos_to$y,
    x_from_lr = lr_pos_from$x,
    y_from_lr = lr_pos_from$y,
    x_to_lr = lr_pos_to$x,
    y_to_lr = lr_pos_to$y,
    x_from_rr = rr_pos_from$x,
    y_from_rr = rr_pos_from$y,
    x_to_rr = rr_pos_to$x,
    y_to_rr = rr_pos_to$y,
    edge_col = edge_fill,
    lwd = lwd,
    weight = weights,
    stringsAsFactors = FALSE
  )

  sender_label_df <- node_df[node_df$col_type == "sender", , drop = FALSE]
  ligand_label_df <- node_df[node_df$col_type == "ligand", , drop = FALSE]
  receptor_label_df <- node_df[node_df$col_type == "receptor", , drop = FALSE]
  receiver_label_df <- node_df[node_df$col_type == "receiver", , drop = FALSE]
  structural_nodes <- node_df[!node_df$is_cell, , drop = FALSE]
  cell_nodes <- node_df[node_df$is_cell, , drop = FALSE]
  sender_label_df$x_label <- sender_label_df$x - 0.08
  receiver_label_df$x_label <- receiver_label_df$x + 0.08
  ligand_label_df$n_lines <- lengths(strsplit(
    ligand_label_df$label_plot,
    "\n",
    fixed = TRUE
  ))
  receptor_label_df$n_lines <- lengths(strsplit(
    receptor_label_df$label_plot,
    "\n",
    fixed = TRUE
  ))
  ligand_label_df$y_label <- ligand_label_df$y +
    0.035 +
    0.018 * pmax(ligand_label_df$n_lines - 1, 0)
  receptor_label_df$y_label <- receptor_label_df$y +
    0.035 +
    0.018 * pmax(receptor_label_df$n_lines - 1, 0)
  y_min <- 0.08
  y_max <- min(
    0.96,
    max(
      c(node_df$y, ligand_label_df$y_label, receptor_label_df$y_label),
      na.rm = TRUE
    ) +
      0.03
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = x_from_sl,
        y = y_from_sl,
        xend = x_to_sl,
        yend = y_to_sl
      ),
      color = edge_df$edge_col,
      linewidth = edge_df$lwd,
      alpha = link_alpha,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = x_from_lr,
        y = y_from_lr,
        xend = x_to_lr,
        yend = y_to_lr
      ),
      color = "grey40",
      linewidth = 0.4,
      linetype = 2,
      alpha = link_alpha,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = x_from_rr,
        y = y_from_rr,
        xend = x_to_rr,
        yend = y_to_rr
      ),
      color = edge_df$edge_col,
      linewidth = edge_df$lwd,
      alpha = link_alpha,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = structural_nodes,
      ggplot2::aes(x = x, y = y),
      fill = structural_nodes$fill,
      color = structural_nodes$border,
      shape = structural_nodes$shape,
      size = node_size,
      alpha = node_alpha,
      stroke = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = cell_nodes,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      fill = cell_nodes$fill,
      color = "grey20",
      stroke = 0.8,
      size = node_size * 0.82,
      alpha = node_alpha,
      show.legend = FALSE
    ) +
    ccc_network_label_layer(
      data = sender_label_df,
      mapping = ggplot2::aes(x = x_label, y = y, label = label_plot),
      label = TRUE,
      label_size = label.size,
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = FALSE,
      hjust = 1
    ) +
    ccc_network_label_layer(
      data = ligand_label_df,
      mapping = ggplot2::aes(x = x, y = y_label, label = label_plot),
      label = TRUE,
      label_size = max(3, label.size * 0.92),
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = FALSE,
      hjust = 0.5,
      vjust = 0,
      lineheight = 0.9,
      fontface = "italic"
    ) +
    ccc_network_label_layer(
      data = receptor_label_df,
      mapping = ggplot2::aes(x = x, y = y_label, label = label_plot),
      label = TRUE,
      label_size = max(3, label.size * 0.92),
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = FALSE,
      hjust = 0.5,
      vjust = 0,
      lineheight = 0.9,
      fontface = "italic"
    ) +
    ccc_network_label_layer(
      data = receiver_label_df,
      mapping = ggplot2::aes(x = x_label, y = y, label = label_plot),
      label = TRUE,
      label_size = label.size,
      label_fg = label.fg,
      label_bg = label.bg,
      label_bg_r = label.bg.r,
      repel = FALSE,
      hjust = 0
    ) +
    ggplot2::geom_point(
      data = data.frame(
        x = NA_real_,
        y = NA_real_,
        label = names(cell_cols),
        fill = unname(cell_cols),
        stringsAsFactors = FALSE
      ),
      ggplot2::aes(x = x, y = y, fill = label),
      shape = 21,
      size = node_size * 0.7,
      show.legend = TRUE,
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = legend.title %||% "Cell type",
      values = cell_cols,
      breaks = names(cell_cols),
      guide = ggplot2::guide_legend(
        override.aes = list(
          shape = 21,
          size = 3,
          color = "grey20",
          stroke = 0.8
        )
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 1, 2, 3),
      labels = c("Sender", "Ligand", "Receptor", "Receiver"),
      expand = ggplot2::expansion(mult = c(0.22, 0.22))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(y_min, y_max),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = subtitle
    )

  if (!is.null(reg.by) && reg.by %in% colnames(df)) {
    reg_df_leg <- data.frame(
      x = NA_real_,
      y = NA_real_,
      reg = reg_levels,
      col = unname(reg_cols[reg_levels]),
      stringsAsFactors = FALSE
    )
    p <- p +
      ggnewscale::new_scale_color() +
      ggplot2::geom_segment(
        data = reg_df_leg,
        ggplot2::aes(
          x = x,
          y = y,
          xend = x,
          yend = y,
          color = reg
        ),
        linewidth = 1,
        show.legend = TRUE,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        name = reg.by,
        values = reg_cols,
        breaks = reg_levels,
        guide = ggplot2::guide_legend(
          override.aes = list(linewidth = 2)
        )
      )
  }

  p <- p +
    do.call(theme_use, theme_args) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction,
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = font.size * 1.15,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(size = font.size),
      plot.margin = ggplot2::margin(5, 40, 5, 40)
    )
  p
}
