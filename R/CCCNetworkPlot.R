#' @title CCC network and flow plots
#'
#' @md
#' @inheritParams CCCStatPlot
#' @inheritParams scop-params
#' @param plot_type Plot type. One of `"circle"`, `"chord"`, `"pathway"`,
#'   `"individual_lr"`, `"arrow"`, `"sigmoid"`, `"bipartite"`,
#'   `"embedding_network"`, `"diff_network"`, or `"spatial"`. The spatial
#'   view overlays stored communication on the selected slice coordinates and
#'   supports `SpatialCellChat`, `SpaTalk`, and `COMMOT` results.
#' @param ligand For `plot_type = "bipartite"`: the ligand name to focus on.
#'   If `NULL`, the ligand with the highest total score is used.
#' @param receptor For `plot_type = "bipartite"`: optional receptor names to
#'   restrict to. If `NULL`, all receptors paired with `ligand` are shown.
#' @param reg.by For `plot_type = "bipartite"`: optional metadata column in
#'   `object` used to color edges by regulation status (e.g. up/down). If `NULL`,
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
#' @param spot_size,spot_alpha Size and alpha of ordinary spatial observations.
#' @param composition_display How composition-mode observations are rendered:
#'   proportional pies, dominant cell type, or neutral points.
#' @param composition_radius Radius of composition pies in display-coordinate
#'   units. `NULL` derives a density-aware default.
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
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
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
  object,
  method = NULL,
  condition = NULL,
  dataset = 1,
  comparison = c(1, 2),
  plot_type = c("circle", "circle_focused", "chord", "lr_chord", "gene_chord", "pathway", "individual_lr", "individual", "individual_outgoing", "individual_incoming", "arrow", "sigmoid", "bipartite", "embedding_network", "diff_network", "spatial", "diffusion"),
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
  measure = c("weight", "count"),
  value = "sum",
  top_n = 20,
  ligand = NULL,
  receptor = NULL,
  reg.by = NULL,
  reg_palette = "Set1",
  reg_palcolor = NULL,
  expr.by = NULL,
  layout = c("circle", "hierarchy", "chord", "kk", "fr", "nicely", "lgl", "mds", "graphopt"),
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
  combine_methods = c("separate", "support", "rank", "legacy"),
  resource = NULL,
  sample = NULL,
  spot_size = 1.2,
  spot_alpha = 0.35,
  composition_display = c("pie", "dominant", "none"),
  composition_radius = NULL,
  ...,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
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
  composition_display <- match.arg(composition_display)
  plot_type_requested <- plot_type
  if (identical(plot_type, "circle") && identical(layout, "chord")) {
    plot_type <- "chord"
  }
  dots <- list(...)
  finish_plot <- function(plot) {
    plot
  }
  finish_base_plot <- function(expr) {
    force(expr)
  }
  label.enable <- isTRUE(dots[["label"]])
  label.spatial <- if (is.null(dots[["label"]])) TRUE else isTRUE(dots[["label"]])
  label.size <- dots[["label.size"]] %||% 4
  label.fg <- dots[["label.fg"]] %||% "white"
  label.bg <- dots[["label.bg"]] %||% "black"
  label.bg.r <- dots[["label.bg.r"]] %||% 0.1
  dots[c("label", "label.size", "label.fg", "label.bg", "label.bg.r")] <- NULL
  ncols <- dots[["ncols"]] %||% dots[["ncol"]] %||% NULL
  nrows <- dots[["nrows"]] %||% dots[["nrow"]] %||% NULL
  combine_panels <- if ("combine" %in% names(dots)) {
    isTRUE(dots[["combine"]])
  } else {
    TRUE
  }
  sender.use <- ccc_alias_arg(dots, "sender_use", sender.use)
  receiver.use <- ccc_alias_arg(dots, "receiver_use", receiver.use)
  ligand.use <- ccc_alias_arg(dots, "ligand_use", ligand.use)
  receptor.use <- ccc_alias_arg(dots, "receptor_use", receptor.use)
  interaction.use <- ccc_alias_arg(dots, "interaction_use", interaction.use)
  pairLR.use <- ccc_alias_arg(dots, "pair_lr_use", pairLR.use)
  thresh <- ccc_alias_arg(dots, "pvalue_threshold", thresh)
  min_interaction_threshold <- ccc_alias_arg(
    dots,
    "min_interaction_threshold",
    NULL
  )
  reduce <- if (is.null(dots[["reduce"]])) TRUE else isTRUE(dots[["reduce"]])
  max.groups <- dots[["max.groups"]] %||% 8
  small.gap <- dots[["small.gap"]] %||% 1
  big.gap <- dots[["big.gap"]] %||% 8
  lab.cex <- dots[["lab.cex"]] %||% 0.6
  if (identical(plot_type, "circle_focused")) {
    ccc_assert_unsupported(
      plot_type = "circle_focused",
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
    if (is.null(signaling)) {
      ccc_assert_unsupported(
        plot_type = "circle_focused",
        sender.use = sender.use,
        receiver.use = receiver.use
      )
    }
    plot_type <- "circle"
    if (!is.null(min_interaction_threshold) && edge_threshold == 0) {
      edge_threshold <- min_interaction_threshold
    }
  }
  if (plot_type %in% c("individual_outgoing", "individual_incoming")) {
    ccc_assert_unsupported(
      plot_type = plot_type,
      signaling = signaling,
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
  }
  if (identical(plot_type, "pathway")) {
    ccc_assert_unsupported(
      plot_type = "pathway",
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
  }
  if (identical(plot_type, "chord")) {
    ccc_assert_unsupported(
      plot_type = "chord",
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
  }
  if (identical(plot_type, "gene_chord")) {
    ccc_assert_unsupported(
      plot_type = "gene_chord",
      interaction.use = interaction.use,
      pairLR.use = pairLR.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
  }
  if (identical(plot_type, "lr_chord")) {
    ccc_assert_unsupported(
      plot_type = "lr_chord",
      signaling = signaling,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
  }
  if (plot_type %in% c("individual", "individual_lr")) {
    ccc_assert_unsupported(
      plot_type = plot_type,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      ligand = ligand,
      receptor = receptor
    )
  }
  if (identical(plot_type, "diffusion")) {
    diffusion_method <- detect_method(srt = srt, method = method)
    if (identical(diffusion_method, "SpatialCellChat")) {
      if (is.null(signaling)) {
        log_message(
          "{.arg signaling} is required for a SpatialCellChat diffusion field",
          message_type = "error"
        )
      }
      check_r("jinworks/SpatialCellChat", verbose = FALSE)
      native <- GetCCCObject(
        object = srt,
        method = "SpatialCellChat",
        result.name = condition,
        sample = sample
      )
      diffusion_fun <- get_namespace_fun(
        "SpatialCellChat",
        "netVisual_CommunField"
      )
      native_formals <- names(formals(diffusion_fun))
      native_dots <- dots[intersect(names(dots), native_formals)]
      native_dots$object <- native
      native_dots$signaling <- signaling
      native_dots$pattern <- native_dots$pattern %||% "outgoing"
      return(do.call(diffusion_fun, native_dots))
    }
    log_message(
      "{.val plot_type = 'diffusion'} requires stored full SpatialCellChat results",
      message_type = "error"
    )
  }
  palette_cfg <- ccc_palettes(
    palette = palette,
    palcolor = palcolor,
    cell_palette = cell_palette,
    cell_palcolor = cell_palcolor,
    link_palette = link_palette,
    link_palcolor = link_palcolor
  )

  method <- detect_method(srt = srt, method = method)
  if (identical(plot_type, "spatial")) {
    stored <- ccc_spatial_stored(
      object = srt,
      method = method,
      condition = condition,
      sample = sample
    )
    pair_df_spatial <- ccc_spatial_network_data(
      stored = stored,
      sender.use = sender.use,
      receiver.use = receiver.use,
      ligand.use = ligand.use,
      receptor.use = receptor.use,
      interaction.use = interaction.use,
      signaling = signaling,
      pairLR.use = pairLR.use,
      thresh = thresh,
      value = value,
      edge_value = edge_value,
      edge_threshold = edge_threshold,
      top_n = top_n
    )
    label_spatial <- label.spatial
    directed_spatial <- if (isTRUE(missing(directed))) TRUE else isTRUE(directed)
    return(finish_plot(ccc_spatial_network_plot(
      stored = stored,
      pair_df = pair_df_spatial,
      cell_palette = palette_cfg$cell_palette,
      cell_palcolor = palette_cfg$cell_palcolor,
      link_palette = palette_cfg$link_palette,
      link_palcolor = palette_cfg$link_palcolor,
      edge_size = edge_size,
      edge_color = edge_color,
      edge_alpha = edge_alpha,
      edge_line = edge_line,
      edge_curvature = edge_curvature,
      directed = directed_spatial,
      arrow_type = arrow_type,
      arrow_angle = arrow_angle,
      arrow_length = arrow_length,
      node_size = node_size,
      node_alpha = node_alpha,
      spot_size = spot_size,
      spot_alpha = spot_alpha,
      composition_display = composition_display,
      composition_radius = composition_radius,
      label = label_spatial,
      label.size = label.size,
      title = title,
      subtitle = subtitle %||% paste(
        stored$method,
        stored$result.name,
        stored$sample %||% "",
        sep = " | "
      ),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      font.size = font.size,
      theme_use = theme_use,
      theme_args = theme_args
    )))
  }
  combine_methods <- match.arg(combine_methods)
  srt <- ccc_prepare_filtered_object(
    srt = srt,
    method = method,
    resource = resource,
    condition = if (identical(method, "CellChat")) NULL else condition,
    sample = sample
  )
  if (identical(method, "CCC") && identical(combine_methods, "separate")) {
    return(ccc_plot_methods_separately(match.call(), srt = srt, env = parent.frame()))
  }
  srt <- ccc_prepare_combined_object(
    srt = srt,
    method = method,
    combine_methods = combine_methods
  )

  if (identical(plot_type, "pathway") && is.null(signaling)) {
    log_message(
      "{.arg signaling} must be provided for {.val plot_type = 'pathway'}",
      message_type = "error"
    )
  }

  if (plot_type %in% c("pathway", "individual", "individual_lr")) {
    if (identical(method, "CellChat")) {
      if (!identical(layout, "circle")) {
        log_message(
          "{.arg layout} is ignored for {.val plot_type = {plot_type}}; CellChat pathway/LR networks use the circle layout",
          message_type = "warning"
        )
      }
      obj_info <- get_dataset_object(
        srt = srt,
        condition = condition,
        dataset = dataset
      )
      cellchat_object <- obj_info$object
      plot_cellchat_circle <- function(sig,
                                       pairLR = NULL,
                                       plot_title = NULL,
                                       plot_subtitle = NULL) {
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
        return(finish_plot(simplify_cc_plot_list(plots)))
      }

      if (is.null(signaling)) {
        log_message(
          "{.arg signaling} must be provided for {.val plot_type = {plot_type_requested}}",
          message_type = "error"
        )
      }
      if (identical(plot_type, "individual_lr") && is.null(pairLR.use)) {
        log_message(
          "{.arg pairLR.use} must be provided for {.val plot_type = 'individual_lr'}",
          message_type = "error"
        )
      }
      plots <- lapply(signaling, function(sig) {
        plot_cellchat_circle(
          sig = sig,
          pairLR = pairLR.use,
          plot_title = title %||% if (is.null(pairLR.use)) {
            sig
          } else {
            paste(
              sig,
              paste(as.character(pairLR.use), collapse = ", "),
              sep = ": "
            )
          },
          plot_subtitle = subtitle
        )
      })
      return(finish_plot(simplify_cc_plot_list(plots)))
    }

    if (!identical(layout, "circle")) {
      log_message(
        "{.arg layout} is ignored for {.val plot_type = {plot_type}} when using pathway-aware generic CCC results; the circle layout is used",
        message_type = "warning"
      )
    }

    long_df <- ccc_long_table_for_method(
      srt = srt,
      method = method,
      condition = condition,
      dataset = dataset,
      slot.name = slot.name,
      thresh = thresh
    )
    long_df <- standardize_long_df(long_df)
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

    if (plot_type %in% c("individual", "individual_lr")) {
      if (is.null(signaling)) {
        log_message(
          "{.arg signaling} must be provided for {.val plot_type = {plot_type_requested}}",
          message_type = "error"
        )
      }
      if (identical(plot_type, "individual_lr") && is.null(pairLR.use)) {
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
    return(finish_plot(simplify_cc_plot_list(plots)))
  }

  if (identical(plot_type, "diff_network")) {
    if (!identical(method, "CellChat")) {
      log_message(
        "{.val plot_type = 'diff_network'} is currently only supported for {.pkg CellChat}",
        message_type = "error"
      )
    }
    layout_use <- if (layout %in% c("hierarchy", "chord")) "circle" else layout
    return(finish_plot(ccc_diff_network_plot(
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
    )))
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

  if (identical(plot_type, "lr_chord")) {
    if (is.null(pairLR.use) && is.null(interaction.use)) {
      log_message(
        "{.arg pairLR.use} or {.arg interaction.use} must be provided for {.val plot_type = 'lr_chord'}",
        message_type = "error"
      )
    }
    dots_chord <- dots
    dots_chord[c("reduce", "max.groups", "small.gap", "big.gap", "lab.cex")] <- NULL
    return(finish_base_plot(do.call(
      ccc_chord_plot,
      c(
        list(
          pair_df = pair_df,
          interaction_df = interaction_df,
          display_by = "interaction",
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
    )))
  }

  if (identical(plot_type, "gene_chord")) {
    return(finish_base_plot(do.call(
      ccc_gene_chord_plot,
      c(
        list(
          df = df,
          top_n = top_n,
          edge_threshold = edge_threshold,
          link_alpha = link_alpha,
          small.gap = small.gap,
          big.gap = big.gap,
          lab.cex = lab.cex
        ),
        network_plot_args
      )
    )))
  }

  if (plot_type %in% c("individual_outgoing", "individual_incoming")) {
    split_var <- if (identical(plot_type, "individual_outgoing")) {
      "sender"
    } else {
      "receiver"
    }
    group_levels <- unique(as.character(df[[split_var]]))
    group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
    if (length(group_levels) == 0L) {
      log_message(
        "No cell groups are available for individual network plotting",
        message_type = "error"
      )
    }
    grid <- ccc_panel_grid(
      n_panels = length(group_levels),
      ncols = ncols %||% min(4L, length(group_levels)),
      nrows = nrows,
      context = paste0(plot_type, " network plot")
    )
    plots <- lapply(group_levels, function(group_i) {
      df_i <- df[df[[split_var]] == group_i, , drop = FALSE]
      pair_df_i <- pair_plot_df(df_i)
      interaction_df_i <- interaction_plot_df(df_i)
      plot_args_i <- network_plot_args
      plot_args_i$title <- title %||% group_i
      do.call(
        ccc_circle_plot,
        c(
          list(
            pair_df = pair_df_i,
            interaction_df = interaction_df_i,
            display_by = display_by,
            top_n = top_n,
            value = value,
            edge_threshold = edge_threshold,
            edge_size = edge_size,
            node_size = node_size,
            node_alpha = node_alpha,
            link_alpha = link_alpha
          ),
          plot_args_i,
          label_plot_args
        )
      )
    })
    return(finish_plot(plot_cc_list(
      plots,
      combine = combine_panels,
      ncol = grid$ncols,
      nrow = grid$nrows
    )))
  }

  if (identical(plot_type, "circle")) {
    return(finish_base_plot(do.call(
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
    )))
  }

  if (identical(plot_type, "chord")) {
    dots_chord <- dots
    dots_chord[c("reduce", "max.groups", "small.gap", "big.gap", "lab.cex")] <- NULL
    return(finish_base_plot(do.call(
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
    )))
  }

  if (plot_type %in% c("arrow", "sigmoid")) {
    return(finish_plot(do.call(
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
    )))
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
    return(finish_plot(bipartite_plot(
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
    )))
  }

  if (identical(plot_type, "embedding_network")) {
    return(finish_plot(do.call(
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
    )))
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
  df <- ccc_mark_significance(df, thresh = thresh)
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
  sender_chr <- as.character(df$sender)
  receptor_chr <- as.character(df$receptor)
  receiver_chr <- as.character(df$receiver)
  df <- df[
    !is.na(sender_chr) &
      nzchar(sender_chr) &
      !is.na(receptor_chr) &
      nzchar(receptor_chr) &
      !is.na(receiver_chr) &
      nzchar(receiver_chr), ,
    drop = FALSE
  ]

  if (!is.null(receptor)) {
    df <- df[
      as.character(df$receptor) %in% as.character(receptor), ,
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
    receptor_scores <- tapply(df$score, df$receptor, sum, na.rm = TRUE)
    receiver_scores <- tapply(df$score, df$receiver, sum, na.rm = TRUE)
    top_senders <- names(sort(sender_scores, decreasing = TRUE))[
      seq_len(min(top_n, length(sender_scores)))
    ]
    top_receptors <- names(sort(receptor_scores, decreasing = TRUE))[
      seq_len(min(top_n, length(receptor_scores)))
    ]
    top_receivers <- names(sort(receiver_scores, decreasing = TRUE))[
      seq_len(min(top_n, length(receiver_scores)))
    ]
    df <- df[
      df$sender %in% top_senders &
        df$receptor %in% top_receptors &
        df$receiver %in% top_receivers, ,
      drop = FALSE
    ]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No CCC records remain after top_n filtering for {.pkg bipartite} plot",
      message_type = "error"
    )
  }

  rank_stage <- function(column) {
    scores <- tapply(df$score, df[[column]], sum, na.rm = TRUE)
    names(sort(scores, decreasing = TRUE))
  }
  senders <- rank_stage("sender")
  receptors <- rank_stage("receptor")
  receivers <- rank_stage("receiver")
  pretty_label <- function(x, node_type = c("cell", "ligand", "receptor")) {
    node_type <- match.arg(node_type)
    x <- as.character(x)
    if (identical(node_type, "cell")) {
      return(x)
    }
    x <- ccc_display_gene(x)
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
    apply_plot_theme(theme_use, theme_args) +
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

ccc_group_levels <- function(x) {
  if (is.factor(x)) {
    lvls <- levels(x)
  } else {
    lvls <- unique(as.character(x))
  }
  lvls[!is.na(lvls) & nzchar(lvls)]
}

ccc_align_named_palcolor <- function(palcolor, levels) {
  if (is.null(palcolor) || is.null(names(palcolor))) {
    return(palcolor)
  }
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]
  pal_names <- names(palcolor)
  if (!all(levels %in% pal_names)) {
    return(palcolor)
  }
  unname(palcolor[levels])
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
  cmp <- cc_get_cmp(srt = srt, condition = condition)
  comp_idx <- cc_resolve_dataset_index(cmp, comparison = comparison)
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

  cell_palcolor <- ccc_align_named_palcolor(cell_palcolor, node_levels)
  link_palcolor <- ccc_align_named_palcolor(link_palcolor, node_levels)
  cell_cols <- palette_colors(
    node_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  link_cols <- palette_colors(
    node_levels,
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
    link_cols[igraph::V(g)$name[edge_start[, 1]]],
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

ccc_reduce_chord_pairs <- function(pair_plot, strength_df, max.groups = 8) {
  keep_cells <- utils::head(strength_df$cell, max.groups)
  reduced <- pair_plot[
    pair_plot$sender %in% keep_cells &
      pair_plot$receiver %in% keep_cells, ,
    drop = FALSE
  ]
  if (nrow(reduced) > 0L) {
    return(list(
      pair_plot = reduced,
      strength_df = strength_df[strength_df$cell %in% keep_cells, , drop = FALSE]
    ))
  }

  keep_cells <- unique(c(
    as.character(pair_plot$sender),
    as.character(pair_plot$receiver)
  ))
  keep_cells <- keep_cells[!is.na(keep_cells) & nzchar(keep_cells)]
  list(
    pair_plot = pair_plot,
    strength_df = strength_df[strength_df$cell %in% keep_cells, , drop = FALSE]
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
        paste(plot_df$interaction_label, plot_df$receiver, sep = "\r") %in% receiver_keep, ,
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
      order(sender_lookup$interaction_label, sender_lookup$receiver, -sender_lookup$weight), ,
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
    reduced <- ccc_reduce_chord_pairs(
      pair_plot = pair_plot,
      strength_df = strength_df,
      max.groups = max.groups
    )
    pair_plot <- reduced$pair_plot
    strength_df <- reduced$strength_df
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
      nzchar(df$receiver), ,
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
      nzchar(edge_df$target), ,
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
    is.finite(edge_df$score) & edge_df$score > 0, ,
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
    order(strength_df$strength, decreasing = TRUE, na.last = TRUE), ,
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
  label_insitu <- dots[["label_insitu"]] %||% TRUE
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

  group_levels <- ccc_group_levels(srt@meta.data[[group.by]])
  cell_palcolor <- ccc_align_named_palcolor(cell_palcolor, group_levels)
  link_palcolor <- ccc_align_named_palcolor(link_palcolor, group_levels)
  base_plot <- do.call(
    CellDimPlot,
    c(
      list(
        object = srt,
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
    levels = group_levels,
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
  label_insitu = TRUE,
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
    object = srt,
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
  cell_palcolor <- ccc_align_named_palcolor(cell_palcolor, group_levels)
  link_palcolor <- ccc_align_named_palcolor(link_palcolor, group_levels)
  colors <- palette_colors(
    group_levels,
    palette = cell_palette,
    palcolor = cell_palcolor,
    NA_keep = TRUE
  )
  edge_cols <- palette_colors(
    group_levels,
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
