#' @title Plot top regulon specificity scores from SCENIC results
#'
#' @description
#' Calculate regulon specificity score (RSS) from SCENIC regulon activity and
#' plot the top regulons for each group.
#'
#' @md
#' @param srt A Seurat object containing SCENIC results from
#' [RunSCENIC()].
#' @param group.by Metadata column used as the cell group annotation.
#' @param tool_name Name of the `srt@tools` entry storing SCENIC results.
#' @param assay Assay used as a fallback source of regulon activity.
#' @param layer Assay layer used as a fallback source of regulon activity.
#' @param plot_type Plot type. `"rss_rank"` keeps the original regulon RSS rank
#' plot. Other options summarize RSS, regulon activity, regulon sizes, or
#' TF-target subnetworks.
#' @param features Optional TF/regulon names used by activity, network, and
#' target plots. Values can match either `"Sox9"` or `"Sox9(+)"`. Explicit
#' values are resolved in input order; duplicated regulons are drawn once.
#' @param reduction Dimensional reduction used when `plot_type =
#' "activity_dim"`. If `NULL`, a UMAP/tSNE/PCA-like reduction is selected when
#' available.
#' @param dims Two reduction dimensions used when `plot_type =
#' "activity_dim"`.
#' @param top_n Number of top regulons labeled for each group.
#' @param activity_scale Whether to z-score each regulon across groups in
#' `plot_type = "activity_heatmap"`. The default is `FALSE` so that the
#' heatmap shows mean regulon activity and does not collapse constant regulons
#' to zero.
#' @param rss_scale Whether to z-score each regulon across groups in
#' `plot_type = "rss_heatmap"`. Use `rss_scale = TRUE`,
#' `activity_scale = TRUE`, and the same `heatmap_limits` value when RSS and
#' activity heatmaps should use a comparable row-wise relative scale.
#' @param heatmap_show_row_names,heatmap_show_column_names Whether to show row
#' and column names in `plot_type = "rss_heatmap"` and `plot_type =
#' "activity_heatmap"`.
#' @param heatmap_cluster_rows,heatmap_cluster_columns Whether to cluster rows
#' and columns in SCENIC heatmaps.
#' @param heatmap_order Row ordering strategy for `plot_type = "rss_heatmap"`
#' and `plot_type = "activity_heatmap"`. `"cluster"` keeps the existing
#' dendrogram-based order, `"group"` groups regulons by the group where each
#' regulon reaches its maximum heatmap value, and `"input"` keeps the resolved
#' feature order. `"group"` and `"input"` disable row clustering so the chosen
#' order is preserved.
#' @param heatmap_row_names_side,heatmap_column_names_side Sides used for row
#' and column names in SCENIC heatmaps.
#' @param heatmap_row_names_rot,heatmap_column_names_rot Rotation angles for
#' row and column names in SCENIC heatmaps.
#' @param heatmap_border Whether to draw heatmap borders in SCENIC heatmaps.
#' @param heatmap_palette,heatmap_palcolor Palette passed to [GroupHeatmap()]
#' or [FeatureHeatmap()] for SCENIC heatmaps. If `heatmap_palette = NULL`, a
#' sensible default is selected for each heatmap type.
#' @param heatmap_group_palette,heatmap_group_palcolor Group annotation palette
#' passed to [GroupHeatmap()] or [FeatureHeatmap()] for SCENIC heatmaps.
#' @param heatmap_limits Optional two-length numeric vector used as the color
#' scale limits for `plot_type = "rss_heatmap"` and `plot_type =
#' "activity_heatmap"`. For example, `c(-2, 2)` fixes both z-score heatmaps to
#' the same legend range.
#' @param heatmap_args Additional arguments passed to [GroupHeatmap()] for
#' `plot_type = "activity_heatmap"` or [FeatureHeatmap()] for `plot_type =
#' "rss_heatmap"`.
#' @param ... Additional arguments passed directly to the underlying
#' [GroupHeatmap()] or [FeatureHeatmap()] call when `plot_type` is
#' `"activity_heatmap"` or `"rss_heatmap"`. For example, `width` and `height`
#' can be supplied directly.
#' @param max_targets Maximum number of target genes shown per TF/regulon in
#' network-style plots.
#' @param max_edges Maximum number of TF-target edges shown in global network
#' plots. Edges are ranked by absolute weight when a weight column is present.
#' @param network_layout Graph layout used by network plots. `"fr"` matches the
#' force-directed Fruchterman-Reingold layout used in Pando examples.
#' @param network_tf Optional TF names used when `plot_type = "network"`. If
#' `NULL`, `features`, `highlight_tf`, or the top RSS regulons are used.
#' @param label_nodes Which nodes to label in network plots.
#' @param network_label_top_n Maximum number of high-degree TF nodes labeled in
#' `plot_type = "network_graph"` when `label_nodes = "tfs"`.
#' @param combine Whether to combine group plots with [patchwork::wrap_plots()].
#' @param ncol Number of columns used when `combine = TRUE`.
#' @param return_data Whether to return RSS matrices and ranking tables together
#' with plots. If `FALSE`, only the plot object or plot list is returned.
#' @param title Optional title added to the combined plot.
#' @param point_color Color for all regulon rank points.
#' @param top_color Color for top regulon rank points.
#' @param point_size Point size.
#' @param point_alpha Alpha for all regulon rank points.
#' @param highlight_tf Optional TF or regulon names to highlight in every group
#' plot. Values can match either `TF` or `regulon`, for example `"Sox9"` or
#' `"Sox9(+)"`.
#' @param highlight_color Color for highlighted TF or regulon points and rank
#' lines.
#' @param highlight_point_size Point size for highlighted TFs or regulons.
#' @param highlight_linewidth Line width for highlighted TF or regulon rank
#' lines.
#' @param label_size Text size for top regulon labels.
#' @param verbose Whether to print messages.
#'
#' @return A list containing `rss_matrix`, `rank_table`, `top_table`, `plots`,
#' and `plot` when `return_data = TRUE`; otherwise a plot object or list of
#' plots.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSCENIC(
#'   pancreas_sub,
#'   species = "Mus_musculus"
#' )
#'
#' scenic_rss <- SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "rss_rank"
#' )
#' scenic_rss$plot
#' example_regulons <- unique(scenic_rss$top_table$regulon)[1:2]
#' example_tfs <- unique(scenic_rss$top_table$TF)[1:2]
#'
#' SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "rss_heatmap")
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "rss_heatmap",
#'   width = 2,
#'   height = 3
#' )
#' SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "rss_dotplot")
#' SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "activity_heatmap")
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "rss_heatmap",
#'   heatmap_order = "group",
#'   heatmap_cluster_columns = FALSE
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "rss_heatmap",
#'   rss_scale = TRUE,
#'   heatmap_order = "group",
#'   heatmap_limits = c(-2, 2)
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "activity_heatmap",
#'   activity_scale = TRUE,
#'   heatmap_limits = c(-2, 2)
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "activity_violin",
#'   features = example_regulons
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "activity_dim",
#'   features = example_regulons
#' )
#' SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "regulon_size")
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "network_graph",
#'   max_targets = 10,
#'   max_edges = 500,
#'   label_nodes = "tfs"
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "network",
#'   network_tf = example_tfs,
#'   max_targets = 30
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "target_bar",
#'   features = example_regulons,
#'   max_targets = 20
#' )
#' }
SCENICPlot <- function(
  srt,
  group.by,
  tool_name = "SCENIC",
  assay = "scenic",
  layer = "data",
  plot_type = c(
    "rss_rank",
    "rss_heatmap",
    "rss_dotplot",
    "activity_heatmap",
    "activity_violin",
    "activity_dim",
    "regulon_size",
    "network_graph",
    "network",
    "target_bar"
  ),
  features = NULL,
  reduction = NULL,
  dims = c(1, 2),
  top_n = 12,
  activity_scale = FALSE,
  rss_scale = FALSE,
  heatmap_show_row_names = FALSE,
  heatmap_show_column_names = FALSE,
  heatmap_cluster_rows = TRUE,
  heatmap_cluster_columns = FALSE,
  heatmap_order = c("cluster", "group", "input"),
  heatmap_row_names_side = "right",
  heatmap_column_names_side = "top",
  heatmap_row_names_rot = 0,
  heatmap_column_names_rot = 45,
  heatmap_border = TRUE,
  heatmap_palette = NULL,
  heatmap_palcolor = NULL,
  heatmap_group_palette = "Chinese",
  heatmap_group_palcolor = NULL,
  heatmap_limits = NULL,
  heatmap_args = list(),
  max_targets = 20,
  max_edges = Inf,
  network_layout = c("fr", "nicely", "kk", "lgl", "drl"),
  network_tf = NULL,
  label_nodes = c("tfs", "all", "none"),
  network_label_top_n = 60,
  combine = TRUE,
  ncol = 3,
  return_data = TRUE,
  title = NULL,
  point_color = "#1F77B4",
  top_color = "#DC050C",
  point_size = 2,
  point_alpha = 0.5,
  highlight_tf = NULL,
  highlight_color = "#7A0177",
  highlight_point_size = 2,
  highlight_linewidth = 0.5,
  label_size = 3,
  verbose = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)
  network_layout <- match.arg(network_layout)
  label_nodes <- match.arg(label_nodes)
  heatmap_order <- match.arg(heatmap_order)
  dot_args <- list(...)
  if (length(dot_args) > 0) {
    if (is.null(names(dot_args)) || any(!nzchar(names(dot_args)))) {
      log_message(
        "Arguments passed through {.arg ...} must be named.",
        message_type = "error"
      )
    }
    if (!plot_type %in% c("rss_heatmap", "activity_heatmap")) {
      log_message(
        "{.arg ...} is only passed to SCENIC heatmap plot types: {.val rss_heatmap} and {.val activity_heatmap}.",
        message_type = "error"
      )
    }
    heatmap_args[names(dot_args)] <- dot_args
  }
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(group.by) != 1 || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  if (length(top_n) != 1 || is.na(top_n) || top_n < 1) {
    log_message(
      "{.arg top_n} must be a positive integer",
      message_type = "error"
    )
  }
  top_n <- as.integer(top_n)
  max_targets <- max(1L, as.integer(max_targets))
  max_edges <- suppressWarnings(as.numeric(max_edges))
  if (length(max_edges) != 1 || is.na(max_edges) || max_edges <= 0) {
    max_edges <- Inf
  }
  network_label_top_n <- max(0L, as.integer(network_label_top_n))
  if (length(dims) != 2 || any(is.na(as.integer(dims))) || any(as.integer(dims) < 1L)) {
    log_message(
      "{.arg dims} must contain two positive dimensions",
      message_type = "error"
    )
  }
  dims <- as.integer(dims)
  if (!is.null(heatmap_limits)) {
    heatmap_limits <- suppressWarnings(as.numeric(heatmap_limits))
    if (length(heatmap_limits) != 2L ||
      any(!is.finite(heatmap_limits)) ||
      heatmap_limits[[1]] >= heatmap_limits[[2]]) {
      log_message(
        "{.arg heatmap_limits} must be an increasing two-length numeric vector",
        message_type = "error"
      )
    }
  }
  if (!is.null(highlight_tf)) {
    highlight_tf <- unique(as.character(highlight_tf))
    highlight_tf <- highlight_tf[!is.na(highlight_tf) & nzchar(highlight_tf)]
    if (length(highlight_tf) == 0) {
      highlight_tf <- NULL
    }
  }

  auc_mat <- scenic_get_rss_auc_matrix(
    srt = srt,
    tool_name = tool_name,
    assay = assay,
    layer = layer
  )
  group_annotation <- srt@meta.data[[group.by]]
  names(group_annotation) <- rownames(srt@meta.data)

  common_cells <- intersect(colnames(auc_mat), names(group_annotation))
  if (length(common_cells) == 0) {
    log_message(
      "No shared cells between SCENIC regulon activity and {.arg srt} metadata",
      message_type = "error"
    )
  }
  auc_mat <- auc_mat[, common_cells, drop = FALSE]
  group_annotation <- group_annotation[common_cells]

  keep_cells <- !is.na(group_annotation)
  if (!any(keep_cells)) {
    log_message(
      "All cells have missing values in {.arg group.by}",
      message_type = "error"
    )
  }
  auc_mat <- auc_mat[, keep_cells, drop = FALSE]
  group_annotation <- group_annotation[keep_cells]
  group_names <- if (is.factor(group_annotation)) {
    levels(droplevels(group_annotation))
  } else {
    unique(as.character(group_annotation))
  }
  group_annotation <- as.character(group_annotation)
  names(group_annotation) <- colnames(auc_mat)

  log_message(
    "Calculating SCENIC RSS for {.val {nrow(auc_mat)}} regulons across {.val {length(group_names)}} group{?s}",
    verbose = verbose
  )
  rss_matrix <- scenic_calc_rss_matrix(
    auc_mat = auc_mat,
    cell_annotation = group_annotation,
    cell_types = group_names
  )
  rss_matrix <- rss_matrix[stats::complete.cases(rss_matrix), , drop = FALSE]
  if (nrow(rss_matrix) == 0) {
    log_message(
      "No valid RSS values remain after removing missing scores",
      message_type = "error"
    )
  }

  rank_table <- do.call(
    rbind,
    lapply(colnames(rss_matrix), function(one_group) {
      specificity_score <- rss_matrix[, one_group]
      keep_regulons <- !is.na(specificity_score)
      regulons <- rownames(rss_matrix)[keep_regulons]
      specificity_score <- specificity_score[keep_regulons]
      regulon_order <- order(specificity_score, decreasing = TRUE)
      regulons <- regulons[regulon_order]
      specificity_score <- specificity_score[regulon_order]
      data.frame(
        group = one_group,
        regulon = regulons,
        TF = sub("\\(\\+\\)$", "", regulons),
        specificity_score = as.numeric(specificity_score),
        rank = seq_along(regulons),
        is_top = seq_along(regulons) <= min(top_n, length(regulons)),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(rank_table) <- NULL
  rank_table[["is_highlight"]] <- FALSE
  if (!is.null(highlight_tf)) {
    rank_table[["is_highlight"]] <- rank_table[["regulon"]] %in%
      highlight_tf |
      rank_table[["TF"]] %in% highlight_tf
  }
  top_table <- rank_table[rank_table[["is_top"]], , drop = FALSE]

  plot_result <- switch(
    plot_type,
    rss_rank = scenic_plot_rss_rank(
      rss_matrix = rss_matrix,
      rank_table = rank_table,
      top_table = top_table,
      highlight_tf = highlight_tf,
      combine = combine,
      ncol = ncol,
      title = title,
      point_color = point_color,
      top_color = top_color,
      point_size = point_size,
      point_alpha = point_alpha,
      highlight_color = highlight_color,
      highlight_point_size = highlight_point_size,
      highlight_linewidth = highlight_linewidth,
      label_size = label_size
    ),
    rss_heatmap = scenic_plot_rss_heatmap(
      rss_matrix = rss_matrix,
      top_table = top_table,
      features = features,
      scale = rss_scale,
      title = title,
      show_row_names = heatmap_show_row_names,
      show_column_names = heatmap_show_column_names,
      cluster_rows = heatmap_cluster_rows,
      cluster_columns = heatmap_cluster_columns,
      heatmap_order = heatmap_order,
      row_names_side = heatmap_row_names_side,
      column_names_side = heatmap_column_names_side,
      row_names_rot = heatmap_row_names_rot,
      column_names_rot = heatmap_column_names_rot,
      border = heatmap_border,
      heatmap_palette = heatmap_palette,
      heatmap_palcolor = heatmap_palcolor,
      group_palette = heatmap_group_palette,
      group_palcolor = heatmap_group_palcolor,
      heatmap_limits = heatmap_limits,
      heatmap_args = heatmap_args
    ),
    rss_dotplot = scenic_plot_rss_dotplot(
      rss_matrix = rss_matrix,
      top_table = top_table,
      features = features,
      title = title
    ),
    activity_heatmap = scenic_plot_activity_heatmap(
      srt = srt,
      auc_mat = auc_mat,
      group_annotation = group_annotation,
      group_names = group_names,
      group.by = group.by,
      top_table = top_table,
      features = features,
      assay = assay,
      layer = layer,
      scale = activity_scale,
      title = title,
      show_row_names = heatmap_show_row_names,
      show_column_names = heatmap_show_column_names,
      cluster_rows = heatmap_cluster_rows,
      cluster_columns = heatmap_cluster_columns,
      heatmap_order = heatmap_order,
      row_names_side = heatmap_row_names_side,
      column_names_side = heatmap_column_names_side,
      row_names_rot = heatmap_row_names_rot,
      column_names_rot = heatmap_column_names_rot,
      border = heatmap_border,
      heatmap_palette = heatmap_palette,
      heatmap_palcolor = heatmap_palcolor,
      group_palette = heatmap_group_palette,
      group_palcolor = heatmap_group_palcolor,
      heatmap_limits = heatmap_limits,
      heatmap_args = heatmap_args
    ),
    activity_violin = scenic_plot_activity_violin(
      srt = srt,
      auc_mat = auc_mat,
      group_annotation = group_annotation,
      group.by = group.by,
      top_table = top_table,
      features = features,
      assay = assay,
      layer = layer,
      combine = combine,
      ncol = ncol,
      title = title
    ),
    activity_dim = scenic_plot_activity_dim(
      srt = srt,
      auc_mat = auc_mat,
      top_table = top_table,
      features = features,
      assay = assay,
      layer = layer,
      reduction = reduction,
      dims = dims,
      combine = combine,
      ncol = ncol,
      title = title,
      point_size = point_size,
      point_alpha = point_alpha
    ),
    regulon_size = scenic_plot_regulon_size(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      top_n = top_n,
      title = title
    ),
    network_graph = scenic_plot_network_graph(
      srt = srt,
      tool_name = tool_name,
      features = features,
      highlight_tf = highlight_tf,
      max_targets = max_targets,
      max_edges = max_edges,
      network_layout = network_layout,
      label_nodes = label_nodes,
      network_label_top_n = network_label_top_n,
      title = title
    ),
    network = scenic_plot_network(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      network_tf = network_tf,
      highlight_tf = highlight_tf,
      max_targets = max_targets,
      label_nodes = label_nodes,
      title = title
    ),
    target_bar = scenic_plot_target_bar(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      max_targets = max_targets,
      title = title
    )
  )

  plots <- plot_result[["plots"]] %||% list(plot_result[["plot"]])
  plot <- plot_result[["plot"]]

  if (isFALSE(return_data)) {
    return(plot)
  }

  list(
    rss_matrix = rss_matrix,
    rank_table = rank_table,
    top_table = top_table,
    plots = plots,
    plot = plot,
    plot_type = plot_type,
    plot_data = plot_result[["data"]]
  )
}

scenic_plot_rss_rank <- function(
  rss_matrix,
  rank_table,
  top_table,
  highlight_tf = NULL,
  combine = TRUE,
  ncol = 3,
  title = NULL,
  point_color = "#1F77B4",
  top_color = "#DC050C",
  point_size = 2,
  point_alpha = 0.5,
  highlight_color = "#7A0177",
  highlight_point_size = 2,
  highlight_linewidth = 0.5,
  label_size = 3
) {
  plots <- lapply(colnames(rss_matrix), function(one_group) {
    data_rank_plot <- rank_table[
      rank_table[["group"]] == one_group,
      ,
      drop = FALSE
    ]
    top_df <- top_table[top_table[["group"]] == one_group, , drop = FALSE]
    highlight_df <- data_rank_plot[
      data_rank_plot[["is_highlight"]],
      ,
      drop = FALSE
    ]
    label_df <- unique(rbind(top_df, highlight_df))

    plot_title <- one_group
    if (!is.null(highlight_tf)) {
      highlight_title <- if (nrow(highlight_df) > 0) {
        paste0(highlight_df[["TF"]], " rank = ", highlight_df[["rank"]])
      } else {
        paste0(highlight_tf, " not found")
      }
      plot_title <- paste0(
        one_group,
        "\n",
        paste(highlight_title, collapse = "; ")
      )
    }

    p <- ggplot2::ggplot(
      data_rank_plot,
      ggplot2::aes(x = .data[["rank"]], y = .data[["specificity_score"]])
    ) +
      ggplot2::geom_point(
        size = point_size,
        shape = 16,
        color = point_color,
        alpha = point_alpha
      ) +
      ggplot2::geom_point(
        data = top_df,
        size = point_size,
        color = top_color
      ) +
      scenic_plot_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = "Regulon rank",
        y = "Specificity Score",
        title = plot_title
      )

    if (nrow(highlight_df) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = highlight_df,
          size = highlight_point_size,
          color = highlight_color
        ) +
        ggplot2::geom_vline(
          data = highlight_df,
          ggplot2::aes(xintercept = .data[["rank"]]),
          linetype = "dashed",
          color = highlight_color,
          linewidth = highlight_linewidth
        )
    }

    p +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(label = .data[["TF"]]),
        color = "black",
        size = label_size,
        fontface = "italic",
        arrow = grid::arrow(ends = "first", length = grid::unit(0.01, "npc")),
        box.padding = 0.2,
        point.padding = 0.3,
        segment.color = "black",
        segment.size = 0.3,
        force = 1,
        max.iter = 3000
      )
  })
  names(plots) <- colnames(rss_matrix)

  plot <- scenic_combine_plots(plots, combine = combine, ncol = ncol, title = title)
  list(plot = plot, plots = plots, data = rank_table)
}

scenic_plot_rss_heatmap <- function(
  rss_matrix,
  top_table,
  features = NULL,
  scale = FALSE,
  title = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  heatmap_order = c("cluster", "group", "input"),
  row_names_side = "right",
  column_names_side = "top",
  row_names_rot = 0,
  column_names_rot = 45,
  border = TRUE,
  heatmap_palette = "viridis",
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  heatmap_limits = NULL,
  heatmap_args = list()
) {
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = rownames(rss_matrix),
    top_table = top_table
  )
  heatmap_order <- match.arg(heatmap_order)
  rss_subset <- rss_matrix[regulons, , drop = FALSE]
  value_name <- "RSS"
  legend_title <- "RSS"
  heatmap_palette <- heatmap_palette %||% if (isTRUE(scale)) "RdBu" else "viridis"
  if (isTRUE(scale)) {
    rss_subset <- scenic_scale_rows(rss_subset)
    constant_rows <- attr(rss_subset, "constant_rows")
    if (length(constant_rows) == nrow(rss_subset)) {
      log_message(
        "All selected regulons have identical group RSS values; z-score heatmap values are zero. Use {.code rss_scale = FALSE} to show raw RSS.",
        message_type = "warning"
      )
    }
    value_name <- "RSS_z"
    legend_title <- "RSS z"
  }
  regulons <- scenic_order_heatmap_features(rss_subset, regulons, heatmap_order)
  rss_subset <- rss_subset[regulons, , drop = FALSE]
  if (!identical(heatmap_order, "cluster")) {
    cluster_rows <- FALSE
  }
  plot_data <- scenic_matrix_to_long(
    rss_subset,
    row_name = "regulon",
    col_name = "group",
    value_name = value_name
  )
  plot_data[["TF"]] <- scenic_tf_from_regulon(plot_data[["regulon"]])
  plot_data[["regulon"]] <- factor(plot_data[["regulon"]], levels = rev(regulons))
  plot_data[["group"]] <- factor(plot_data[["group"]], levels = colnames(rss_subset))

  plot <- scenic_plot_feature_heatmap_from_matrix(
    mat = rss_subset,
    group_names = colnames(rss_subset),
    features = regulons,
    legend_title = legend_title,
    title = title %||% "SCENIC regulon specificity",
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_names_side = row_names_side,
    column_names_side = column_names_side,
    row_names_rot = row_names_rot,
    column_names_rot = column_names_rot,
    border = border,
    heatmap_palette = heatmap_palette,
    heatmap_palcolor = heatmap_palcolor,
    group_palette = group_palette,
    group_palcolor = group_palcolor,
    heatmap_limits = heatmap_limits,
    heatmap_args = heatmap_args
  )

  list(plot = plot, plots = list(plot), data = plot_data)
}

scenic_plot_rss_dotplot <- function(rss_matrix, top_table, features = NULL, title = NULL) {
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = rownames(rss_matrix),
    top_table = top_table
  )
  rss_subset <- rss_matrix[regulons, , drop = FALSE]
  plot_data <- scenic_matrix_to_long(
    rss_subset,
    row_name = "regulon",
    col_name = "group",
    value_name = "RSS"
  )
  plot_data[["regulon"]] <- factor(plot_data[["regulon"]], levels = rev(regulons))
  plot_data[["group"]] <- factor(plot_data[["group"]], levels = colnames(rss_subset))

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[["group"]],
      y = .data[["regulon"]],
      size = .data[["RSS"]],
      color = .data[["RSS"]]
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradientn(colors = c("#2C7BB6", "#FFFFBF", "#D7191C")) +
    ggplot2::scale_size(range = c(1, 6)) +
    scenic_plot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, size = "RSS", color = "RSS", title = title %||% "SCENIC RSS dot plot")

  list(plot = plot, plots = list(plot), data = plot_data)
}

scenic_plot_activity_heatmap <- function(
  srt,
  auc_mat,
  group_annotation,
  group_names,
  group.by,
  top_table,
  features = NULL,
  assay = "scenic",
  layer = "data",
  scale = FALSE,
  title = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  heatmap_order = c("cluster", "group", "input"),
  row_names_side = "right",
  column_names_side = "top",
  row_names_rot = 0,
  column_names_rot = 45,
  border = TRUE,
  heatmap_palette = NULL,
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  heatmap_limits = NULL,
  heatmap_args = list()
) {
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = rownames(auc_mat),
    top_table = top_table
  )
  heatmap_order <- match.arg(heatmap_order)
  avg_mat <- scenic_group_average_matrix(auc_mat[regulons, , drop = FALSE], group_annotation, group_names)
  feature_split <- scenic_align_heatmap_feature_split(
    feature_split = heatmap_args[["feature_split"]],
    features = features,
    regulons = regulons,
    available = rownames(auc_mat)
  )
  value_name <- "activity"
  if (isTRUE(scale)) {
    avg_mat <- scenic_scale_rows(avg_mat)
    constant_rows <- attr(avg_mat, "constant_rows")
    if (length(constant_rows) == nrow(avg_mat)) {
      log_message(
        "All selected regulons have identical group-average activity; z-score heatmap values are zero. Use {.code activity_scale = FALSE} to show mean activity.",
        message_type = "warning"
      )
    }
    value_name <- "activity_z"
  }
  regulons <- scenic_order_heatmap_features(avg_mat, regulons, heatmap_order)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[regulons]
    heatmap_args[["feature_split"]] <- feature_split
  }
  avg_mat <- avg_mat[regulons, , drop = FALSE]
  if (!identical(heatmap_order, "cluster")) {
    cluster_rows <- FALSE
  }
  plot_data <- scenic_matrix_to_long(
    avg_mat,
    row_name = "regulon",
    col_name = "group",
    value_name = value_name
  )
  plot_data[["regulon"]] <- factor(plot_data[["regulon"]], levels = rev(regulons))
  plot_data[["group"]] <- factor(plot_data[["group"]], levels = colnames(avg_mat))

  srt_use <- scenic_attach_auc_assay(srt = srt, auc_mat = auc_mat, assay = assay)
  heatmap_palette <- heatmap_palette %||% if (isTRUE(scale)) "RdBu" else "viridis"
  plot <- scenic_call_with_args(
    GroupHeatmap,
    args = list(
      srt = srt_use,
      features = regulons,
      group.by = group.by,
      aggregate_fun = base::mean,
      border = border,
      layer = layer,
      assay = assay,
      exp_method = if (isTRUE(scale)) "zscore" else "raw",
      exp_legend_title = if (isTRUE(scale)) "Activity z" else "Mean activity",
      limits = heatmap_limits,
      lib_normalize = FALSE,
      max_cells = Inf,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      column_title = title %||% "SCENIC regulon activity",
      heatmap_palette = heatmap_palette,
      heatmap_palcolor = heatmap_palcolor,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      verbose = FALSE
    ),
    extra_args = heatmap_args
  )

  list(plot = plot, plots = list(plot), data = plot_data)
}

scenic_plot_activity_violin <- function(
  srt,
  auc_mat,
  group_annotation,
  group.by,
  top_table,
  features = NULL,
  assay = "scenic",
  layer = "data",
  combine = TRUE,
  ncol = 3,
  title = NULL
) {
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = rownames(auc_mat),
    top_table = top_table,
    max_features = if (is.null(features)) 6 else NULL
  )
  plot_data <- do.call(
    rbind,
    lapply(regulons, function(regulon) {
      data.frame(
        regulon = regulon,
        cell = colnames(auc_mat),
        group = group_annotation[colnames(auc_mat)],
        activity = as.numeric(auc_mat[regulon, ]),
        stringsAsFactors = FALSE
      )
    })
  )
  plot_data[["group"]] <- factor(plot_data[["group"]], levels = unique(group_annotation))

  srt_use <- scenic_attach_auc_assay(srt = srt, auc_mat = auc_mat, assay = assay)
  plot <- FeatureStatPlot(
    srt = srt_use,
    stat.by = regulons,
    group.by = group.by,
    assay = assay,
    layer = layer,
    plot_type = "violin",
    add_box = TRUE,
    ylab = "Regulon activity",
    title = title,
    combine = combine,
    ncol = ncol,
    force = TRUE
  )
  plots <- scenic_as_plot_list(plot, regulons)
  list(plot = plot, plots = plots, data = plot_data)
}

scenic_plot_activity_dim <- function(
  srt,
  auc_mat,
  top_table,
  features = NULL,
  assay = "scenic",
  layer = "data",
  reduction = NULL,
  dims = c(1, 2),
  combine = TRUE,
  ncol = 3,
  title = NULL,
  point_size = 2,
  point_alpha = 0.8
) {
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = rownames(auc_mat),
    top_table = top_table,
    max_features = if (is.null(features)) 6 else NULL
  )
  reduction <- scenic_select_reduction(srt, reduction)
  emb <- Seurat::Embeddings(srt, reduction = reduction)
  if (max(dims) > ncol(emb)) {
    log_message(
      "{.arg dims} requests dimension {.val {max(dims)}}, but reduction {.val {reduction}} has only {.val {ncol(emb)}} dimension{?s}",
      message_type = "error"
    )
  }
  srt_use <- scenic_attach_auc_assay(srt = srt, auc_mat = auc_mat, assay = assay)
  plots <- FeatureDimPlot(
    srt = srt_use,
    features = regulons,
    assay = assay,
    layer = layer,
    reduction = reduction,
    dims = dims,
    pt.size = point_size,
    pt.alpha = point_alpha,
    legend.title = "Activity",
    combine = FALSE,
    force = TRUE
  )
  plot <- scenic_combine_plots(plots, combine = combine, ncol = ncol, title = title)
  list(plot = plot, plots = plots, data = data.frame(regulon = regulons))
}

scenic_plot_regulon_size <- function(srt, tool_name, top_table, features = NULL, top_n = 12, title = NULL) {
  regulon_list <- scenic_get_regulon_list(srt, tool_name)
  size_data <- data.frame(
    regulon = names(regulon_list),
    TF = scenic_tf_from_regulon(names(regulon_list)),
    target_count = lengths(regulon_list),
    stringsAsFactors = FALSE
  )
  if (!is.null(features)) {
    regulons <- scenic_resolve_regulon_features(features, size_data[["regulon"]], top_table = NULL)
    size_data <- size_data[size_data[["regulon"]] %in% regulons, , drop = FALSE]
  } else {
    size_data <- size_data[order(size_data[["target_count"]], decreasing = TRUE), , drop = FALSE]
    size_data <- utils::head(size_data, top_n)
  }
  size_data[["regulon"]] <- factor(size_data[["regulon"]], levels = rev(size_data[["regulon"]]))

  plot <- ggplot2::ggplot(
    size_data,
    ggplot2::aes(x = .data[["target_count"]], y = .data[["regulon"]])
  ) +
    ggplot2::geom_col(fill = "#E69F00", width = 0.75) +
    scenic_plot_theme() +
    ggplot2::labs(x = "Target genes", y = NULL, title = title %||% "SCENIC regulon size")

  list(plot = plot, plots = list(plot), data = size_data)
}

scenic_plot_network_graph <- function(
  srt,
  tool_name,
  features = NULL,
  highlight_tf = NULL,
  max_targets = 30,
  max_edges = Inf,
  network_layout = c("fr", "nicely", "kk", "lgl", "drl"),
  label_nodes = c("tfs", "all", "none"),
  network_label_top_n = 60,
  title = NULL
) {
  network_layout <- match.arg(network_layout)
  label_nodes <- match.arg(label_nodes)
  adjacency <- scenic_get_adjacency(srt, tool_name)
  cols <- scenic_adjacency_columns(adjacency)
  if (!is.null(features)) {
    features <- unique(c(as.character(features), scenic_tf_from_regulon(features)))
    adjacency <- adjacency[
      adjacency[[cols[["tf"]]]] %in% features |
        adjacency[[cols[["target"]]]] %in% features,
      ,
      drop = FALSE
    ]
  }
  if (nrow(adjacency) == 0) {
    log_message(
      "No SCENIC adjacency edges found for the requested global network",
      message_type = "error"
    )
  }
  adjacency <- scenic_top_edges(adjacency, tf_col = cols[["tf"]], weight_col = cols[["weight"]], max_targets = max_targets)
  adjacency <- scenic_limit_edges(adjacency, weight_col = cols[["weight"]], max_edges = max_edges)
  network_data <- scenic_network_plot_data(adjacency = adjacency, cols = cols, layout = network_layout)
  node_data <- network_data[["nodes"]]
  edge_plot <- network_data[["edge_plot"]]
  edge_data <- network_data[["edges"]]
  label_data <- scenic_network_label_data(
    node_data = node_data,
    label_nodes = label_nodes,
    highlight_tf = highlight_tf,
    top_n = network_label_top_n
  )

  plot <- scenic_network_ggplot(
    node_data = node_data,
    edge_plot = edge_plot,
    label_data = label_data,
    title = title %||% "SCENIC TF-target network",
    edge_width_range = c(0.12, 0.8),
    node_size_range = c(1.2, 8.5)
  )

  list(plot = plot, plots = list(plot), data = list(edges = edge_data, nodes = node_data))
}

scenic_plot_network <- function(
  srt,
  tool_name,
  top_table,
  features = NULL,
  network_tf = NULL,
  highlight_tf = NULL,
  max_targets = 30,
  label_nodes = c("tfs", "all", "none"),
  title = NULL
) {
  label_nodes <- match.arg(label_nodes)
  adjacency <- scenic_get_adjacency(srt, tool_name)
  cols <- scenic_adjacency_columns(adjacency)
  tf_candidates <- network_tf %||% features %||% highlight_tf
  if (is.null(tf_candidates)) {
    tf_candidates <- unique(top_table[["TF"]])
  }
  tf_candidates <- unique(scenic_tf_from_regulon(as.character(tf_candidates)))
  adjacency <- adjacency[adjacency[[cols[["tf"]]]] %in% tf_candidates, , drop = FALSE]
  if (nrow(adjacency) == 0) {
    log_message(
      "No SCENIC adjacency edges found for requested TFs",
      message_type = "error"
    )
  }
  adjacency <- scenic_top_edges(adjacency, tf_col = cols[["tf"]], weight_col = cols[["weight"]], max_targets = max_targets)
  network_data <- scenic_network_plot_data(adjacency = adjacency, cols = cols, layout = "fr")
  edge_data <- network_data[["edges"]]
  node_data <- network_data[["nodes"]]
  edge_plot <- network_data[["edge_plot"]]

  label_data <- switch(
    label_nodes,
    all = node_data,
    tfs = node_data[node_data[["node_type"]] == "TF", , drop = FALSE],
    none = node_data[FALSE, , drop = FALSE]
  )
  root_data <- node_data[node_data[["name"]] %in% unique(edge_data[["from"]]), , drop = FALSE]
  label_data <- unique(rbind(label_data, root_data))

  plot <- scenic_network_ggplot(
    node_data = node_data,
    edge_plot = edge_plot,
    label_data = label_data,
    title = title %||% "SCENIC TF-target network",
    edge_width_range = c(0.15, 1.2),
    node_size_range = c(2, 6)
  )

  list(plot = plot, plots = list(plot), data = list(edges = edge_data, nodes = node_data))
}

scenic_plot_target_bar <- function(srt, tool_name, top_table, features = NULL, max_targets = 30, title = NULL) {
  regulon_list <- scenic_get_regulon_list(srt, tool_name)
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = names(regulon_list),
    top_table = top_table,
    max_features = 4
  )
  adjacency <- tryCatch(scenic_get_adjacency(srt, tool_name), error = function(...) NULL)
  cols <- if (!is.null(adjacency)) scenic_adjacency_columns(adjacency) else NULL

  plot_data <- do.call(
    rbind,
    lapply(regulons, function(regulon) {
      tf <- scenic_tf_from_regulon(regulon)
      targets <- regulon_list[[regulon]]
      df <- data.frame(regulon = regulon, TF = tf, target = targets, importance = NA_real_, stringsAsFactors = FALSE)
      if (!is.null(adjacency)) {
        adj_tf <- adjacency[adjacency[[cols[["tf"]]]] == tf & adjacency[[cols[["target"]]]] %in% targets, , drop = FALSE]
        if (nrow(adj_tf) > 0 && !is.null(cols[["weight"]])) {
          df[["importance"]] <- adj_tf[[cols[["weight"]]]][match(df[["target"]], adj_tf[[cols[["target"]]]])]
        }
      }
      df <- df[order(df[["importance"]], decreasing = TRUE, na.last = TRUE), , drop = FALSE]
      utils::head(df, max_targets)
    })
  )
  plot_data[["target"]] <- factor(plot_data[["target"]], levels = rev(unique(plot_data[["target"]])))
  value_col <- if (all(is.na(plot_data[["importance"]]))) "rank_score" else "importance"
  if (identical(value_col, "rank_score")) {
    plot_data[["rank_score"]] <- rev(seq_len(nrow(plot_data)))
  }

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[value_col]], y = .data[["target"]])
  ) +
    ggplot2::geom_col(fill = "#E69F00", width = 0.75) +
    scenic_plot_theme() +
    ggplot2::facet_wrap(ggplot2::vars(.data[["regulon"]]), scales = "free_y") +
    ggplot2::labs(
      x = if (identical(value_col, "importance")) "GRNBoost2 importance" else "Target rank",
      y = NULL,
      title = title %||% "SCENIC regulon targets"
    )

  list(plot = plot, plots = list(plot), data = plot_data)
}

scenic_network_plot_data <- function(adjacency, cols, layout = "fr") {
  edge_data <- data.frame(
    from = adjacency[[cols[["tf"]]]],
    to = adjacency[[cols[["target"]]]],
    weight = if (!is.null(cols[["weight"]])) adjacency[[cols[["weight"]]]] else 1,
    stringsAsFactors = FALSE
  )
  edge_data <- edge_data[edge_data[["from"]] != edge_data[["to"]], , drop = FALSE]
  edge_data <- edge_data[stats::complete.cases(edge_data[, c("from", "to"), drop = FALSE]), , drop = FALSE]
  if (nrow(edge_data) == 0) {
    log_message(
      "No valid TF-target edges remain for SCENIC network plotting",
      message_type = "error"
    )
  }
  graph <- igraph::graph_from_data_frame(edge_data, directed = TRUE)
  xy <- scenic_network_layout(graph, layout = layout)
  node_data <- data.frame(
    name = igraph::V(graph)$name,
    x = xy[, 1],
    y = xy[, 2],
    degree = igraph::degree(graph, mode = "all"),
    out_degree = igraph::degree(graph, mode = "out"),
    in_degree = igraph::degree(graph, mode = "in"),
    stringsAsFactors = FALSE
  )
  node_data[["node_type"]] <- ifelse(node_data[["out_degree"]] > 0, "TF", "target")
  edge_plot <- merge(edge_data, node_data[, c("name", "x", "y")], by.x = "from", by.y = "name")
  edge_plot <- merge(edge_plot, node_data[, c("name", "x", "y")], by.x = "to", by.y = "name", suffixes = c("", "_end"))
  edge_plot[["edge_sign"]] <- ifelse(edge_plot[["weight"]] < 0, "negative", "positive")
  list(edges = edge_data, nodes = node_data, edge_plot = edge_plot, graph = graph)
}

scenic_network_layout <- function(graph, layout = "fr") {
  layout_weights <- NULL
  if ("weight" %in% igraph::edge_attr_names(graph)) {
    layout_weights <- abs(igraph::E(graph)$weight)
    layout_weights[!is.finite(layout_weights) | layout_weights <= 0] <- 1
  }
  xy <- switch(
    layout,
    fr = igraph::layout_with_fr(graph, weights = layout_weights),
    nicely = igraph::layout_nicely(graph),
    kk = igraph::layout_with_kk(graph, weights = layout_weights),
    lgl = igraph::layout_with_lgl(graph),
    drl = igraph::layout_with_drl(graph)
  )
  if (ncol(xy) < 2) {
    xy <- cbind(xy[, 1], 0)
  }
  xy
}

scenic_network_label_data <- function(
  node_data,
  label_nodes = c("tfs", "all", "none"),
  highlight_tf = NULL,
  top_n = 60
) {
  label_nodes <- match.arg(label_nodes)
  label_data <- switch(
    label_nodes,
    all = node_data,
    tfs = {
      tf_nodes <- node_data[node_data[["node_type"]] == "TF", , drop = FALSE]
      tf_nodes <- tf_nodes[order(tf_nodes[["degree"]], decreasing = TRUE), , drop = FALSE]
      utils::head(tf_nodes, top_n)
    },
    none = node_data[FALSE, , drop = FALSE]
  )
  if (!is.null(highlight_tf)) {
    highlight_tf <- scenic_tf_from_regulon(highlight_tf)
    highlight_data <- node_data[node_data[["name"]] %in% highlight_tf, , drop = FALSE]
    label_data <- unique(rbind(label_data, highlight_data))
  }
  label_data
}

scenic_network_ggplot <- function(
  node_data,
  edge_plot,
  label_data,
  title,
  edge_width_range = c(0.12, 0.8),
  node_size_range = c(1.2, 8.5)
) {
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_plot,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        xend = .data[["x_end"]],
        yend = .data[["y_end"]],
        color = .data[["edge_sign"]],
        linewidth = abs(.data[["weight"]])
      ),
      alpha = 0.75,
      arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = node_data,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["node_type"]],
        size = .data[["degree"]]
      ),
      shape = 21,
      color = "white",
      stroke = 0.35
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["name"]]),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      segment.color = "grey60",
      segment.size = 0.2
    ) +
    ggplot2::scale_color_manual(
      values = c(negative = "darkgrey", positive = "#E69F00"),
      breaks = c("positive", "negative"),
      labels = c(positive = "positive", negative = "negative"),
      drop = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = c(target = "#2C7BB6", TF = "#E69F00"),
      breaks = c("target", "TF"),
      labels = c(target = "target", TF = "TF"),
      name = NULL
    ) +
    ggplot2::scale_size_continuous(
      range = node_size_range,
      breaks = scenic_degree_breaks(node_data[["degree"]]),
      name = "Degree"
    ) +
    ggplot2::scale_linewidth(range = edge_width_range, guide = "none") +
    ggplot2::coord_equal() +
    scenic_plot_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 11)
    ) +
    ggplot2::labs(color = "Edge", title = title)
}

scenic_degree_breaks <- function(degree) {
  degree <- degree[is.finite(degree)]
  if (length(degree) == 0) {
    return(NULL)
  }
  max_degree <- max(degree)
  breaks <- pretty(c(1, max_degree), n = 3)
  breaks <- breaks[breaks > 0 & breaks <= max_degree]
  unique(breaks)
}

scenic_plot_theme <- function() {
  theme_scop() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 12),
      axis.text = ggplot2::element_text(colour = "black", size = 10),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

scenic_plot_feature_heatmap_from_matrix <- function(
  mat,
  group_names,
  features,
  legend_title,
  title = NULL,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  column_names_side = "top",
  row_names_rot = 0,
  column_names_rot = 45,
  border = TRUE,
  heatmap_palette = "viridis",
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  heatmap_limits = NULL,
  heatmap_args = list()
) {
  colnames(mat) <- group_names
  mat_sparse <- Matrix::Matrix(mat, sparse = TRUE)
  srt_heatmap <- Seurat::CreateSeuratObject(
    counts = mat_sparse,
    assay = "SCENICHeatmap"
  )
  assay_object <- SeuratObject::SetAssayData(
    object = srt_heatmap[["SCENICHeatmap"]],
    layer = "data",
    new.data = mat_sparse
  )
  srt_heatmap[["SCENICHeatmap"]] <- assay_object
  srt_heatmap[["SCENIC_group"]] <- factor(group_names, levels = group_names)
  scenic_call_with_args(
    FeatureHeatmap,
    args = list(
      srt = srt_heatmap,
      features = features,
      cells = group_names,
      group.by = "SCENIC_group",
      max_cells = Inf,
      cell_order = group_names,
      border = border,
      layer = "data",
      assay = "SCENICHeatmap",
      exp_method = "raw",
      exp_legend_title = legend_title,
      limits = heatmap_limits,
      lib_normalize = FALSE,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      column_title = title,
      heatmap_palette = heatmap_palette,
      heatmap_palcolor = heatmap_palcolor,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      verbose = FALSE
    ),
    extra_args = heatmap_args
  )
}

scenic_call_with_args <- function(fun, args, extra_args = list()) {
  extra_args <- extra_args %||% list()
  if (!is.list(extra_args)) {
    log_message(
      "{.arg heatmap_args} must be a list",
      message_type = "error"
    )
  }
  args[names(extra_args)] <- extra_args
  do.call(fun, args)
}

scenic_scale_rows <- function(mat) {
  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
  variable_rows <- is.finite(row_sd) & row_sd > sqrt(.Machine$double.eps)
  out <- mat
  if (any(variable_rows)) {
    out[variable_rows, ] <- sweep(
      sweep(mat[variable_rows, , drop = FALSE], 1, row_mean[variable_rows], "-"),
      1,
      row_sd[variable_rows],
      "/"
    )
  }
  if (any(!variable_rows)) {
    out[!variable_rows, ] <- 0
  }
  attr(out, "constant_rows") <- rownames(mat)[!variable_rows]
  out
}

scenic_attach_auc_assay <- function(srt, auc_mat, assay = "scenic") {
  cells <- intersect(colnames(srt), colnames(auc_mat))
  if (length(cells) == 0) {
    log_message(
      "No shared cells between SCENIC regulon activity and {.arg srt}",
      message_type = "error"
    )
  }
  if (!identical(cells, colnames(srt))) {
    srt <- srt[, cells, drop = FALSE]
  }
  auc_mat <- auc_mat[, colnames(srt), drop = FALSE]
  auc_mat <- Matrix::Matrix(auc_mat, sparse = TRUE)
  assay_object <- Seurat::CreateAssayObject(
    counts = auc_mat,
    check.matrix = FALSE
  )
  assay_object <- SeuratObject::SetAssayData(
    object = assay_object,
    layer = "data",
    new.data = auc_mat
  )
  srt[[assay]] <- assay_object
  srt
}

scenic_as_plot_list <- function(plot, names_use = NULL) {
  if (is.list(plot) && !inherits(plot, c("ggplot", "patchwork"))) {
    return(plot)
  }
  plots <- list(plot)
  if (!is.null(names_use) && length(names_use) == length(plots)) {
    names(plots) <- names_use
  }
  plots
}

scenic_combine_plots <- function(plots, combine = TRUE, ncol = 3, title = NULL) {
  plot <- plots
  if (isTRUE(combine)) {
    check_r("patchwork", verbose = FALSE)
    plot <- if (length(plots) == 1) {
      plots[[1]]
    } else {
      patchwork::wrap_plots(plotlist = plots, ncol = ncol)
    }
    if (!is.null(title)) {
      plot <- plot + patchwork::plot_annotation(title = title)
    }
  }
  plot
}

scenic_matrix_to_long <- function(mat, row_name, col_name, value_name) {
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(out) <- c(row_name, col_name, value_name)
  out
}

scenic_tf_from_regulon <- function(regulon) {
  sub("\\(\\+\\)$", "", as.character(regulon))
}

scenic_resolve_regulon_features <- function(
  features,
  available,
  top_table = NULL,
  max_features = NULL
) {
  if (is.null(features)) {
    if (!is.null(top_table) && nrow(top_table) > 0) {
      features <- unique(top_table[["regulon"]])
    } else {
      features <- available
    }
  }
  features <- unique(as.character(features))
  regulons <- scenic_match_regulon_features(features, available)
  if (length(regulons) == 0) {
    log_message(
      "None of {.arg features} matched available SCENIC regulons",
      message_type = "error"
    )
  }
  if (!is.null(max_features)) {
    regulons <- utils::head(regulons, max_features)
  }
  regulons
}

scenic_match_regulon_features <- function(features, available) {
  matches <- vapply(
    as.character(features),
    function(feature) {
      candidates <- unique(c(feature, paste0(feature, "(+)")))
      hit <- candidates[candidates %in% available]
      if (length(hit) == 0) {
        return(NA_character_)
      }
      hit[[1]]
    },
    character(1)
  )
  unique(matches[!is.na(matches)])
}

scenic_align_heatmap_feature_split <- function(
  feature_split,
  features,
  regulons,
  available
) {
  if (is.null(feature_split)) {
    return(NULL)
  }
  if (!is.null(names(feature_split)) && all(regulons %in% names(feature_split))) {
    return(feature_split[regulons])
  }
  if (length(feature_split) == length(regulons) && is.null(features)) {
    names(feature_split) <- regulons
    return(feature_split)
  }
  if (is.null(features) || length(feature_split) != length(features)) {
    log_message(
      "{.arg feature_split} must have the same length as the explicit {.arg features} used by {.val activity_heatmap}, or be named by regulon.",
      message_type = "error"
    )
  }

  split_levels <- if (is.factor(feature_split)) {
    levels(feature_split)
  } else {
    unique(as.character(feature_split))
  }
  split_values <- as.character(feature_split)
  matched_all <- vapply(
    as.character(features),
    function(feature) {
      candidates <- unique(c(feature, paste0(feature, "(+)")))
      hit <- candidates[candidates %in% available]
      if (length(hit) == 0) {
        return(NA_character_)
      }
      hit[[1]]
    },
    character(1)
  )

  out <- stats::setNames(rep(NA_character_, length(regulons)), regulons)
  for (idx in seq_along(matched_all)) {
    regulon <- matched_all[[idx]]
    if (!is.na(regulon) && regulon %in% names(out) && is.na(out[[regulon]])) {
      out[[regulon]] <- split_values[[idx]]
    }
  }
  if (any(is.na(out[regulons]))) {
    missing_regulons <- regulons[is.na(out[regulons])]
    log_message(
      "{.arg feature_split} could not be matched to all displayed regulons: {.val {missing_regulons}}",
      message_type = "error"
    )
  }
  stats::setNames(factor(out[regulons], levels = split_levels), regulons)
}

scenic_order_heatmap_features <- function(
  mat,
  features,
  heatmap_order = c("cluster", "group", "input")
) {
  heatmap_order <- match.arg(heatmap_order)
  if (!identical(heatmap_order, "group")) {
    return(features)
  }

  mat <- as.matrix(mat[features, , drop = FALSE])
  mat[!is.finite(mat)] <- NA_real_
  max_group_idx <- apply(mat, 1L, function(x) {
    if (all(is.na(x))) {
      return(NA_integer_)
    }
    which.max(replace(x, is.na(x), -Inf))
  })
  max_value <- apply(mat, 1L, function(x) {
    out <- suppressWarnings(max(x, na.rm = TRUE))
    if (is.finite(out)) out else NA_real_
  })
  max_group <- colnames(mat)[max_group_idx]
  features[
    order(
      factor(max_group, levels = colnames(mat)),
      -max_value,
      features,
      na.last = TRUE
    )
  ]
}

scenic_group_average_matrix <- function(auc_mat, group_annotation, group_names) {
  out <- vapply(
    group_names,
    function(group) {
      cells <- names(group_annotation)[group_annotation == group]
      Matrix::rowMeans(auc_mat[, cells, drop = FALSE])
    },
    numeric(nrow(auc_mat))
  )
  rownames(out) <- rownames(auc_mat)
  colnames(out) <- group_names
  out
}

scenic_select_reduction <- function(srt, reduction = NULL) {
  reductions <- SeuratObject::Reductions(srt)
  if (!is.null(reduction)) {
    if (!reduction %in% reductions) {
      log_message(
        "{.arg reduction} {.val {reduction}} is not present in {.arg srt}",
        message_type = "error"
      )
    }
    return(reduction)
  }
  preferred <- c("StandardUMAP2D", "umap", "UMAP", "tsne", "TSNE", "pca", "PCA")
  reduction <- scenic_first(intersect(preferred, reductions), scenic_first(reductions))
  if (is.null(reduction)) {
    log_message(
      "No dimensional reduction found in {.arg srt}",
      message_type = "error"
    )
  }
  reduction
}

scenic_get_regulon_list <- function(srt, tool_name) {
  regulon_list <- srt@tools[[tool_name]][["regulon_list"]]
  if (!is.null(regulon_list)) {
    return(regulon_list)
  }
  regulons <- srt@tools[[tool_name]][["regulons"]]
  if (is.null(regulons) || !all(c("regulon", "target") %in% colnames(regulons))) {
    log_message(
      "Cannot find SCENIC regulon target lists in tools slot {.val {tool_name}}",
      message_type = "error"
    )
  }
  scenic_regulon_list <- lapply(seq_len(nrow(regulons)), function(idx) {
    unique(unlist(strsplit(regulons[["target"]][[idx]], ",", fixed = TRUE), use.names = FALSE))
  })
  names(scenic_regulon_list) <- regulons[["regulon"]]
  scenic_regulon_list
}

scenic_get_adjacency <- function(srt, tool_name) {
  adjacency <- srt@tools[[tool_name]][["adjacency"]]
  if (is.null(adjacency) || nrow(adjacency) == 0) {
    log_message(
      "Cannot find SCENIC adjacency table in tools slot {.val {tool_name}}",
      message_type = "error"
    )
  }
  as.data.frame(adjacency, stringsAsFactors = FALSE)
}

scenic_adjacency_columns <- function(adjacency) {
  cols <- colnames(adjacency)
  tf_col <- scenic_first(intersect(c("TF", "tf", "regulator", "source", "from"), cols), cols[[1]])
  target_col <- scenic_first(intersect(c("target", "Target", "gene", "to"), cols), cols[[2]])
  weight_col <- scenic_first(
    intersect(c("importance", "weight", "estimate", "score"), cols),
    if (length(cols) >= 3) cols[[3]] else NULL
  )
  list(tf = tf_col, target = target_col, weight = weight_col)
}

scenic_top_edges <- function(adjacency, tf_col, weight_col, max_targets = 30) {
  edge_list <- split(adjacency, adjacency[[tf_col]])
  edge_list <- lapply(edge_list, function(df) {
    if (!is.null(weight_col) && weight_col %in% colnames(df)) {
      df <- df[order(abs(df[[weight_col]]), decreasing = TRUE), , drop = FALSE]
    }
    utils::head(df, max_targets)
  })
  do.call(rbind, edge_list)
}

scenic_limit_edges <- function(adjacency, weight_col = NULL, max_edges = Inf) {
  if (!is.finite(max_edges) || nrow(adjacency) <= max_edges) {
    return(adjacency)
  }
  if (!is.null(weight_col) && weight_col %in% colnames(adjacency)) {
    adjacency <- adjacency[order(abs(adjacency[[weight_col]]), decreasing = TRUE), , drop = FALSE]
  }
  utils::head(adjacency, max_edges)
}

scenic_first <- function(x, default = NULL) {
  if (length(x) == 0 || is.null(x)) {
    return(default)
  }
  x[[1]]
}

scenic_get_rss_auc_matrix <- function(
  srt,
  tool_name = "SCENIC",
  assay = "scenic",
  layer = "data"
) {
  if (!is.null(srt@tools[[tool_name]][["scores_cells_by_regulon"]])) {
    scores_cells_by_regulon <- srt@tools[[tool_name]][[
      "scores_cells_by_regulon"
    ]]
    auc_mat <- as.matrix(scores_cells_by_regulon)
    if (is.null(rownames(auc_mat)) || is.null(colnames(auc_mat))) {
      log_message(
        "{.arg scores_cells_by_regulon} must have cell row names and regulon column names",
        message_type = "error"
      )
    }
    return(t(auc_mat))
  }

  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "Cannot find SCENIC results in tools slot {.val {tool_name}} or assay {.val {assay}}",
      message_type = "error"
    )
  }
  auc_mat <- GetAssayData5(srt, assay = assay, layer = layer)
  auc_mat <- as.matrix(auc_mat)
  if (is.null(rownames(auc_mat)) || is.null(colnames(auc_mat))) {
    log_message(
      "{.arg assay} regulon activity matrix must have regulon row names and cell column names",
      message_type = "error"
    )
  }
  auc_mat
}

scenic_calc_rss_matrix <- function(
  auc_mat,
  cell_annotation,
  cell_types = NULL
) {
  if (any(is.na(cell_annotation))) {
    log_message(
      "{.arg cell_annotation} contains missing values",
      message_type = "error"
    )
  }
  if (is.null(cell_types)) {
    cell_types <- unique(cell_annotation)
  }

  row_sums <- rowSums(auc_mat)
  norm_auc <- auc_mat / row_sums
  rss_list <- lapply(cell_types, function(this_type) {
    p_cell_type <- as.numeric(cell_annotation == this_type)
    p_cell_type <- p_cell_type / sum(p_cell_type)
    vapply(
      seq_len(nrow(norm_auc)),
      function(regulon_idx) {
        scenic_calc_one_rss(norm_auc[regulon_idx, ], p_cell_type)
      },
      numeric(1)
    )
  })
  rss <- do.call(cbind, rss_list)
  rownames(rss) <- rownames(auc_mat)
  colnames(rss) <- cell_types
  rss
}

scenic_calc_one_rss <- function(p_regulon, p_cell_type) {
  jsd <- scenic_calc_jsd(p_regulon, p_cell_type)
  1 - sqrt(jsd)
}

scenic_calc_jsd <- function(p_regulon, p_cell_type) {
  scenic_entropy((p_regulon + p_cell_type) / 2) -
    ((scenic_entropy(p_regulon) + scenic_entropy(p_cell_type)) / 2)
}

scenic_entropy <- function(p_vector) {
  p_vector <- p_vector[p_vector > 0]
  -sum(p_vector * log2(p_vector))
}
