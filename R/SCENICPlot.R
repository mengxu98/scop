#' @title Plot SCENIC or SCENIC+ regulon activity and specificity
#'
#' @description
#' Visualize regulon specificity scores (RSS) and activity stored by
#' [RunSCENIC()] or [RunSCENICPlus()]. Use [SCENICPlusPlot()] for SCENIC+
#' defaults (`tool_name = "SCENICPlus"`, `assay = "scenicplus"`).
#' `"heatmap_dotplot"` matches the official SCENIC+ signature figure (color =
#' TF expression, size = RSS). `"eregulon_dim"` compares TF expression with
#' gene-based AUC, and region-based AUC when it is stored. `"coverage"` draws
#' accessibility and region–gene link tracks for one target locus.
#' `"network"` draws per-TF regulon hubs after iRegulon/SCENIC+ concentrical
#' figures (TF at the center, targets on a ring; SCENIC+ adds an inner region
#' ring when triplets are stored). `"egrn"` is the SCENIC+ enhancer-driven GRN
#' (Nature Methods 2023 Fig. 2e: TF, unlabeled region diamonds, and gene
#' nodes). `"overlap"` is the eRegulon target-gene Jaccard heatmap.
#'
#' @md
#' @inheritParams CellDimPlot
#' @param srt A Seurat object containing SCENIC results from
#' [RunSCENIC()] or SCENIC+ results from [RunSCENICPlus()].
#' @param group.by Metadata column used as the cell group annotation.
#' @param tool_name Name of the `srt@tools` entry storing SCENIC results.
#' @param assay Assay used as a fallback source of regulon activity.
#' @param layer Assay layer used as a fallback source of regulon activity.
#' @param plot_type Plot type. `"rss_rank"` keeps the original regulon RSS rank
#' plot. `"heatmap_dotplot"` colors dots by TF expression and sizes them by
#' RSS. `"eregulon_dim"` compares TF expression and regulon AUC on embeddings.
#' `"coverage"` shows ATAC accessibility with region–gene links for a target
#' locus. `"network"` draws regulon hub networks. `"network_graph"` draws a
#' global TF–target graph. `"egrn"` draws the SCENIC+ TF–region–gene GRN.
#' `"overlap"` shows pairwise eRegulon target overlap. Other options summarize
#' RSS, regulon activity, or regulon sizes.
#' @param features Optional TF/regulon names used by activity, network, and
#' target plots. Values can match either `"Sox9"` or `"Sox9(+)"`. Explicit
#' values are resolved in input order; duplicated regulons are drawn once.
#' @param reduction Dimensional reduction used when `plot_type =
#' "activity_dim"`. If `NULL`, a UMAP/tSNE/PCA-like reduction is selected when
#' available.
#' @param dims Two reduction dimensions used when `plot_type =
#' "activity_dim"`.
#' @param cor.features Two metadata columns or gene names used when
#' `plot_type = "activity_cor_dumbbell"`. SCENIC regulon activity is
#' correlated with each feature across matched cells.
#' @param cor.feature.labels Optional labels for `cor.features` in the plot
#' legend. If `NULL`, `cor.features` are used.
#' @param cor_method Correlation method used by
#' `plot_type = "activity_cor_dumbbell"`.
#' @param p_cutoff P-value cutoff used to mark significant correlations in
#' `plot_type = "activity_cor_dumbbell"`.
#' @param cor_sort.by Ordering used for the dumbbell rows. `"difference"`
#' sorts by the absolute distance between the two correlations, `"max_abs"`
#' sorts by the strongest absolute correlation, and `"input"` keeps the
#' resolved regulon order.
#' @param cor_cols Two colors used for the two correlation targets.
#' @param cor_xlim Optional x-axis limits for the dumbbell plot.
#' @param cor_label Whether to label dumbbell points with correlation
#' coefficients.
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
#' regulon reaches its maximum heatmap value. For `activity_heatmap` without
#' explicit `features`, `"group"` uses the RSS group that first selected each
#' displayed regulon, so `top_n` controls both the displayed feature set and
#' source-group row blocks. `"input"` keeps the resolved feature order.
#' `"group"` and `"input"` disable row clustering so the chosen order is
#' preserved.
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
#' @param palette,palcolor Palette passed to `palette_colors()` for
#' `"rss_dotplot"`, `"heatmap_dotplot"`, `"regulon_size"`, `"target_bar"`,
#' and network TF/edge colors. Network plots follow published SCENIC+/iRegulon
#' figures: TF-colored edges, light-gray gene nodes, unlabeled region diamonds,
#' and concentric or Kamada–Kawai layouts. Multi-TF `"egrn"` and `"network_graph"`
#' plots add a legend of TF colors labeled with the RSS top cell type
#' (`Jun (Ductal)`). Single-TF hubs put that cell type in the panel title.
#' When `palette = "RdYlBu"` (the `SCENICPlot()` default), networks use the
#' scop `"Chinese"` discrete palette.
#' @param compare_expression Whether `"activity_dim"` also plots TF gene
#' expression beside regulon activity. `"eregulon_dim"` always compares
#' expression with activity.
#' @param expression_assay,expression_layer Assay and layer used for TF
#' expression in `"heatmap_dotplot"`, `"activity_dim"`, and `"eregulon_dim"`.
#' If `expression_assay` is `NULL`, `"RNA"` is used when present.
#' @param expression_scale Whether `"heatmap_dotplot"` z-scores TF expression
#' across groups.
#' @param dim_args Additional arguments passed to [FeatureDimPlot()] and
#' [CellDimPlot()] when `plot_type` is `"activity_dim"` or `"eregulon_dim"`.
#' @param atac_assay Chromatin assay used by `plot_type = "coverage"`. If
#' `NULL`, `"peaks"` or `"ATAC"` is used when present.
#' @param extend.upstream,extend.downstream Bases added around the selected
#' region–gene locus in `"coverage"` plots.
#' @param coverage_args Additional arguments passed to [CoverageTrackPlot()]
#' when fragment files are available for `"coverage"` plots.
#' @param ... Additional arguments passed directly to the underlying
#' [GroupHeatmap()] or [FeatureHeatmap()] call when `plot_type` is
#' `"activity_heatmap"`, `"rss_heatmap"`, or `"overlap"`. For example, `width`
#' and `height` can be supplied directly.
#' @param max_targets Maximum number of target genes shown per TF/regulon in
#' network-style plots.
#' @param max_edges Maximum number of TF-target edges shown in global network
#' plots. Edges are ranked by absolute weight when a weight column is present.
#' @param network_layout Graph layout used by network plots. `"auto"` selects
#' `"star"` for `"network"` (and single-TF `"egrn"`), `"kk"` for
#' `"network_graph"` and multi-TF `"egrn"`. `"star"` is the SCENIC+
#' concentrical / iRegulon hub (TF at the origin; inner ring = regions when
#' triplets are present). `"hub"` places TFs on a circle with unique targets
#' around each TF. `"tripartite"` is a column layout (TF | region | gene).
#' `"fr"` and `"kk"` are force-directed layouts initialized from a circle.
#' @param network_tf Optional TF names used when `plot_type` is `"network"` or
#' `"egrn"`. If `NULL`, `features`, `highlight_tf`, or the top RSS regulons
#' are used.
#' @param network_include_regions Whether SCENIC+ `"network"` / `"network_graph"`
#' plots use TF–region–gene triplets when they are stored. Default is `TRUE`.
#' @param label_nodes Which nodes to label in network plots. `"auto"` labels
#' TFs plus target genes in star/eGRN plots (region coordinates stay unlabeled)
#' and high-degree TFs in dense global graphs.
#' @param network_label_top_n Maximum number of high-degree TF nodes labeled in
#' `plot_type = "network_graph"` when `label_nodes = "tfs"`.
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
#' @param label_max_overlaps Maximum number of overlapping labels allowed by
#' [ggrepel::geom_text_repel()] before dropping a label. The default `Inf`
#' keeps all requested top or highlighted TF labels.
#' @param verbose Whether to print messages.
#'
#' @return A list containing `rss_matrix`, `rank_table`, `top_table`, `plots`,
#' and `plot` when `return_data = TRUE`; heatmap plot types also include the
#' full `heatmap` result returned by [FeatureHeatmap()] or [GroupHeatmap()].
#' Otherwise, a plot object or list of plots.
#'
#' @seealso [RunSCENIC], [RunSCENICPlus], [SCENICPlusPlot], [CoverageTrackPlot],
#' [GroupHeatmap], [FeatureHeatmap], [FeatureDimPlot], [FeatureStatPlot]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunSCENIC(
#'   pancreas_sub,
#'   species = "Mus_musculus",
#'   backend = "cpp",
#'   work_dir = "test/scenic"
#' )
#'
#' scenic_rss <- SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "rss_rank"
#' )
#' scenic_rss$plot
#' example_regulons <- unique(scenic_rss$top_table$regulon)[1:3]
#' example_tfs <- unique(scenic_rss$top_table$TF)[1:3]
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = example_regulons,
#'   assay = "scenic",
#'   reduction = "StandardpcaUMAP2D"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = example_regulons,
#'   group.by = "CellType",
#'   assay = "scenic"
#' )
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
#' SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "heatmap_dotplot")
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
#'   features = example_regulons,
#'   compare_expression = TRUE
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "eregulon_dim",
#'   features = example_regulons
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "activity_cor_dumbbell",
#'   features = example_regulons,
#'   cor.features = c("nFeature_RNA", "nCount_RNA")
#' )
#' SCENICPlot(pancreas_sub, group.by = "CellType", plot_type = "regulon_size")
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "network_graph",
#'   features = example_tfs,
#'   max_targets = 12,
#'   max_edges = 200,
#'   label_nodes = "all"
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "network",
#'   features = example_tfs[[1]],
#'   max_targets = 20
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "target_bar",
#'   features = example_regulons,
#'   max_targets = 20
#' )
#' SCENICPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "overlap",
#'   features = example_regulons
#' )
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "egrn",
#'   features = example_tfs,
#'   max_targets = 10
#' )
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "network_graph",
#'   features = example_tfs,
#'   max_targets = 8,
#'   label_nodes = "all"
#' )
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "coverage",
#'   features = example_tfs
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
    "heatmap_dotplot",
    "activity_heatmap",
    "activity_violin",
    "activity_dim",
    "eregulon_dim",
    "activity_cor_dumbbell",
    "regulon_size",
    "network_graph",
    "network",
    "egrn",
    "overlap",
    "target_bar",
    "coverage"
  ),
  features = NULL,
  reduction = NULL,
  dims = c(1, 2),
  cor.features = NULL,
  cor.feature.labels = NULL,
  cor_method = c("spearman", "pearson", "kendall"),
  p_cutoff = 0.05,
  cor_sort.by = c("difference", "max_abs", "input"),
  cor_cols = c("#56B4E9", "#E83F6F"),
  cor_xlim = NULL,
  cor_label = TRUE,
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
  palette = "RdYlBu",
  palcolor = NULL,
  compare_expression = FALSE,
  expression_assay = NULL,
  expression_layer = "data",
  expression_scale = TRUE,
  dim_args = list(),
  atac_assay = NULL,
  extend.upstream = 50000,
  extend.downstream = 50000,
  coverage_args = list(),
  max_targets = 20,
  max_edges = Inf,
  network_layout = c(
    "auto",
    "star",
    "hub",
    "tripartite",
    "bipartite",
    "circle",
    "fr",
    "nicely",
    "kk",
    "lgl",
    "drl"
  ),
  network_tf = NULL,
  network_include_regions = TRUE,
  label_nodes = c("auto", "tfs", "all", "none"),
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
  label_max_overlaps = Inf,
  verbose = TRUE,
  ...
) {
  plot_type <- match.arg(plot_type)
  network_layout <- match.arg(network_layout)
  label_nodes <- match.arg(label_nodes)
  heatmap_order <- match.arg(heatmap_order)
  cor_method <- match.arg(cor_method)
  cor_sort.by <- match.arg(cor_sort.by)
  dot_args <- list(...)
  if (length(dot_args) > 0) {
    if (is.null(names(dot_args)) || any(!nzchar(names(dot_args)))) {
      log_message(
        "Arguments passed through {.arg ...} must be named.",
        message_type = "error"
      )
    }
    if (!plot_type %in% c("rss_heatmap", "activity_heatmap", "overlap")) {
      log_message(
        "{.arg ...} is only passed to SCENIC heatmap plot types: {.val rss_heatmap}, {.val activity_heatmap}, and {.val overlap}.",
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
  if (plot_type == "activity_cor_dumbbell") {
    if (is.null(cor.features) || length(cor.features) != 2L) {
      log_message(
        "{.arg cor.features} must contain exactly two metadata columns or gene names",
        message_type = "error"
      )
    }
    if (length(p_cutoff) != 1L || !is.finite(p_cutoff) || p_cutoff < 0) {
      log_message(
        "{.arg p_cutoff} must be one non-negative finite number",
        message_type = "error"
      )
    }
    if (length(cor_cols) != 2L) {
      log_message(
        "{.arg cor_cols} must contain two colors",
        message_type = "error"
      )
    }
    if (!is.null(cor_xlim)) {
      cor_xlim <- suppressWarnings(as.numeric(cor_xlim))
      if (length(cor_xlim) != 2L ||
        any(!is.finite(cor_xlim)) ||
        cor_xlim[[1]] >= cor_xlim[[2]]) {
        log_message(
          "{.arg cor_xlim} must be an increasing two-length numeric vector",
          message_type = "error"
        )
      }
    }
    if (!is.null(cor.feature.labels) && length(cor.feature.labels) != 2L) {
      log_message(
        "{.arg cor.feature.labels} must contain two labels when provided",
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
        TF = scenic_tf_from_regulon(regulons),
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

  plot_result <- switch(plot_type,
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
      label_size = label_size,
      label_max_overlaps = label_max_overlaps
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
      palette = palette,
      palcolor = palcolor,
      title = title
    ),
    heatmap_dotplot = scenic_plot_heatmap_dotplot(
      srt = srt,
      rss_matrix = rss_matrix,
      top_table = top_table,
      features = features,
      group_annotation = group_annotation,
      group_names = group_names,
      expression_assay = expression_assay,
      expression_layer = expression_layer,
      expression_scale = expression_scale,
      activity_assay = assay,
      palette = palette,
      palcolor = palcolor,
      title = title,
      verbose = verbose
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
      group.by = group.by,
      top_table = top_table,
      features = features,
      assay = assay,
      layer = layer,
      reduction = reduction,
      dims = dims,
      compare_expression = compare_expression,
      expression_assay = expression_assay,
      expression_layer = expression_layer,
      include_region_auc = FALSE,
      tool_name = tool_name,
      palette = palette,
      palcolor = palcolor,
      group_palette = heatmap_group_palette,
      group_palcolor = heatmap_group_palcolor,
      dim_args = dim_args,
      combine = combine,
      ncol = ncol,
      title = title,
      point_size = point_size,
      point_alpha = point_alpha,
      verbose = verbose
    ),
    eregulon_dim = scenic_plot_activity_dim(
      srt = srt,
      auc_mat = auc_mat,
      group.by = group.by,
      top_table = top_table,
      features = features,
      assay = assay,
      layer = layer,
      reduction = reduction,
      dims = dims,
      compare_expression = TRUE,
      expression_assay = expression_assay,
      expression_layer = expression_layer,
      include_region_auc = TRUE,
      tool_name = tool_name,
      palette = palette,
      palcolor = palcolor,
      group_palette = heatmap_group_palette,
      group_palcolor = heatmap_group_palcolor,
      dim_args = dim_args,
      combine = combine,
      ncol = ncol,
      title = title,
      point_size = point_size,
      point_alpha = point_alpha,
      verbose = verbose
    ),
    activity_cor_dumbbell = scenic_plot_activity_cor_dumbbell(
      srt = srt,
      auc_mat = auc_mat,
      top_table = top_table,
      features = features,
      cor.features = cor.features,
      cor.feature.labels = cor.feature.labels,
      cor_method = cor_method,
      p_cutoff = p_cutoff,
      cor_sort.by = cor_sort.by,
      cor_cols = cor_cols,
      cor_xlim = cor_xlim,
      cor_label = cor_label,
      assay = assay,
      layer = layer,
      title = title,
      point_size = point_size,
      label_size = label_size,
      verbose = verbose
    ),
    regulon_size = scenic_plot_regulon_size(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      top_n = top_n,
      palette = palette,
      palcolor = palcolor,
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
      network_include_regions = network_include_regions,
      palette = palette,
      palcolor = palcolor,
      title = title,
      rank_table = rank_table
    ),
    network = scenic_plot_network(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      network_tf = network_tf,
      highlight_tf = highlight_tf,
      max_targets = max_targets,
      network_layout = network_layout,
      label_nodes = label_nodes,
      network_include_regions = network_include_regions,
      palette = palette,
      palcolor = palcolor,
      combine = combine,
      ncol = ncol,
      title = title,
      rank_table = rank_table
    ),
    egrn = scenic_plot_egrn(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      network_tf = network_tf,
      highlight_tf = highlight_tf,
      max_targets = max_targets,
      network_layout = network_layout,
      label_nodes = label_nodes,
      palette = palette,
      palcolor = palcolor,
      title = title,
      rank_table = rank_table
    ),
    overlap = scenic_plot_overlap(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      top_n = top_n,
      heatmap_palette = heatmap_palette,
      heatmap_palcolor = heatmap_palcolor,
      heatmap_group_palette = heatmap_group_palette,
      heatmap_group_palcolor = heatmap_group_palcolor,
      heatmap_limits = heatmap_limits,
      heatmap_args = heatmap_args,
      title = title
    ),
    target_bar = scenic_plot_target_bar(
      srt = srt,
      tool_name = tool_name,
      top_table = top_table,
      features = features,
      max_targets = max_targets,
      palette = palette,
      palcolor = palcolor,
      title = title
    ),
    coverage = scenic_plot_coverage(
      srt = srt,
      tool_name = tool_name,
      group.by = group.by,
      group_annotation = group_annotation,
      group_names = group_names,
      features = features,
      top_table = top_table,
      atac_assay = atac_assay,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      palette = palette,
      palcolor = palcolor,
      group_palette = heatmap_group_palette,
      group_palcolor = heatmap_group_palcolor,
      coverage_args = coverage_args,
      title = title,
      verbose = verbose
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
    plot_data = plot_result[["data"]],
    heatmap = plot_result[["heatmap"]]
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
  label_size = 3,
  label_max_overlaps = Inf
) {
  plots <- lapply(colnames(rss_matrix), function(one_group) {
    data_rank_plot <- rank_table[
      rank_table[["group"]] == one_group, ,
      drop = FALSE
    ]
    top_df <- top_table[top_table[["group"]] == one_group, , drop = FALSE]
    highlight_df <- data_rank_plot[
      data_rank_plot[["is_highlight"]], ,
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
      theme_scop() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(colour = "black", size = 12),
        axis.text = ggplot2::element_text(colour = "black", size = 10),
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
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
        max.iter = 3000,
        max.overlaps = label_max_overlaps
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

  heatmap_result <- scenic_plot_feature_heatmap_from_matrix(
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
  plot <- heatmap_result[["plot"]]

  list(plot = plot, plots = list(plot), data = plot_data, heatmap = heatmap_result)
}

scenic_plot_rss_dotplot <- function(
  rss_matrix,
  top_table,
  features = NULL,
  palette = "RdYlBu",
  palcolor = NULL,
  title = NULL
) {
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
  fill_colors <- scenic_continuous_colors(palette = palette, palcolor = palcolor)

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
    scenic_color_scale(fill_colors, "RSS") +
    ggplot2::scale_size(range = c(1, 6)) +
    theme_scop() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 12),
      axis.text = ggplot2::element_text(colour = "black", size = 10),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, size = "RSS", title = title %||% "SCENIC RSS dot plot")

  list(plot = plot, plots = list(plot), data = plot_data)
}

scenic_plot_heatmap_dotplot <- function(
  srt,
  rss_matrix,
  top_table,
  features = NULL,
  group_annotation,
  group_names,
  expression_assay = NULL,
  expression_layer = "data",
  expression_scale = TRUE,
  activity_assay = "scenic",
  palette = "RdYlBu",
  palcolor = NULL,
  title = NULL,
  verbose = TRUE
) {
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
  plot_data[["TF"]] <- scenic_tf_from_regulon(plot_data[["regulon"]])
  expr_info <- scenic_group_tf_expression(
    srt = srt,
    tfs = unique(plot_data[["TF"]]),
    group_annotation = group_annotation,
    group_names = group_names,
    expression_assay = expression_assay,
    expression_layer = expression_layer,
    activity_assay = activity_assay,
    scale = expression_scale,
    verbose = verbose
  )
  color_name <- "RSS"
  if (!is.null(expr_info)) {
    expr_long <- scenic_matrix_to_long(
      expr_info[["matrix"]],
      row_name = "TF",
      col_name = "group",
      value_name = "TF_expr"
    )
    plot_data <- merge(
      plot_data,
      expr_long,
      by = c("TF", "group"),
      all.x = TRUE
    )
    color_name <- if (isTRUE(expression_scale)) "TF expr z" else "TF expression"
  } else {
    plot_data[["TF_expr"]] <- plot_data[["RSS"]]
    log_message(
      "No matching TF expression found; {.val heatmap_dotplot} colors dots by RSS",
      message_type = "warning",
      verbose = verbose
    )
  }
  plot_data[["regulon"]] <- factor(plot_data[["regulon"]], levels = rev(regulons))
  plot_data[["group"]] <- factor(plot_data[["group"]], levels = colnames(rss_subset))
  fill_colors <- scenic_continuous_colors(palette = palette, palcolor = palcolor)
  plot_title <- title %||% if (!is.null(expr_info)) {
    "TF expression and regulon RSS"
  } else {
    "SCENIC RSS dot plot"
  }

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[["group"]],
      y = .data[["regulon"]],
      size = .data[["RSS"]],
      color = .data[["TF_expr"]]
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    scenic_color_scale(fill_colors, color_name) +
    ggplot2::scale_size(range = c(1, 8), name = "RSS") +
    theme_scop() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 12),
      axis.text = ggplot2::element_text(colour = "black", size = 10),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = plot_title)

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
  rss_source <- scenic_top_regulon_source(top_table, regulons)
  avg_mat <- scenic_group_average_matrix(auc_mat[regulons, , drop = FALSE], group_annotation, group_names)
  feature_split <- scenic_align_heatmap_feature_split(
    feature_split = heatmap_args[["feature_split"]],
    features = features,
    regulons = regulons,
    available = rownames(auc_mat)
  )
  if (
    is.null(feature_split) &&
      is.null(features) &&
      identical(heatmap_order, "group") &&
      !is.null(rss_source)
  ) {
    feature_split <- rss_source[["group"]]
  }
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
  if (identical(heatmap_order, "group") && is.null(features) && !is.null(rss_source)) {
    regulons <- rss_source[["regulon"]]
  } else {
    regulons <- scenic_order_heatmap_features(avg_mat, regulons, heatmap_order)
  }
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
  if (!is.null(rss_source)) {
    source_idx <- match(as.character(plot_data[["regulon"]]), rss_source[["regulon"]])
    plot_data[["rss_group"]] <- rss_source[["group"]][source_idx]
    plot_data[["rss_rank"]] <- rss_source[["rank"]][source_idx]
  }

  srt_use <- scenic_attach_auc_assay(srt = srt, auc_mat = auc_mat, assay = assay)
  heatmap_palette <- heatmap_palette %||% if (isTRUE(scale)) "RdBu" else "viridis"
  heatmap_result <- scenic_call_with_args(
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
  plot <- heatmap_result[["plot"]]

  list(plot = plot, plots = list(plot), data = plot_data, heatmap = heatmap_result)
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
  group.by,
  top_table,
  features = NULL,
  assay = "scenic",
  layer = "data",
  reduction = NULL,
  dims = c(1, 2),
  compare_expression = FALSE,
  expression_assay = NULL,
  expression_layer = "data",
  include_region_auc = FALSE,
  tool_name = "SCENIC",
  palette = "RdYlBu",
  palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  dim_args = list(),
  combine = TRUE,
  ncol = 3,
  title = NULL,
  point_size = 2,
  point_alpha = 0.8,
  verbose = TRUE
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
  assay_features <- tryCatch(rownames(srt_use[[assay]]), error = function(...) character())
  region_auc <- NULL
  if (isTRUE(include_region_auc)) {
    region_auc <- scenic_get_region_auc_matrix(srt, tool_name)
    if (!is.null(region_auc)) {
      region_auc <- region_auc[
        intersect(rownames(region_auc), regulons),
        intersect(colnames(region_auc), colnames(srt_use)),
        drop = FALSE
      ]
      if (nrow(region_auc) == 0L) {
        region_auc <- NULL
      }
    }
  }
  if (!isTRUE(compare_expression) && is.null(region_auc)) {
    plots <- FeatureDimPlot(
      srt = srt_use,
      features = scenic_assay_features(regulons, assay_features),
      assay = assay,
      layer = layer,
      reduction = reduction,
      dims = dims,
      palette = palette,
      palcolor = palcolor,
      pt.size = point_size,
      pt.alpha = point_alpha,
      legend.title = "Activity",
      combine = FALSE,
      force = TRUE
    )
    plot <- scenic_combine_plots(plots, combine = combine, ncol = ncol, title = title)
    return(list(plot = plot, plots = plots, data = data.frame(regulon = regulons)))
  }

  expression_assay <- scenic_expression_assay(
    srt = srt_use,
    expression_assay = expression_assay,
    activity_assay = assay
  )
  expression_layer <- scenic_resolve_expression_layer(
    srt = srt_use,
    assay = expression_assay,
    layer = expression_layer
  )
  expr_names <- tryCatch(
    rownames(srt_use[[expression_assay]]),
    error = function(...) character()
  )
  cell_formals <- names(formals(CellDimPlot))
  feat_formals <- names(formals(FeatureDimPlot))
  cell_args <- list(
    srt = srt_use,
    group.by = group.by,
    reduction = reduction,
    dims = dims,
    palette = group_palette,
    palcolor = group_palcolor,
    label = TRUE,
    title = "Cell types",
    pt.size = point_size,
    pt.alpha = point_alpha
  )
  extra_cell <- dim_args[intersect(names(dim_args), cell_formals)]
  extra_cell <- extra_cell[setdiff(names(extra_cell), c("srt", "group.by"))]
  cell_args[names(extra_cell)] <- extra_cell
  p_cell <- do.call(CellDimPlot, cell_args)

  tf_plots <- lapply(regulons, function(regulon) {
    assay_feature <- scenic_assay_feature(regulon, assay_features)
    if (is.na(assay_feature)) {
      log_message(
        "Regulon {.val {regulon}} is not present in assay {.val {assay}}",
        message_type = "error"
      )
    }
    act_args <- list(
      srt = srt_use,
      features = assay_feature,
      assay = assay,
      layer = layer,
      reduction = reduction,
      dims = dims,
      palette = palette,
      palcolor = palcolor,
      bg_cutoff = -Inf,
      title = paste(regulon, "activity"),
      legend.title = "Activity",
      pt.size = point_size,
      pt.alpha = point_alpha,
      combine = TRUE,
      force = TRUE
    )
    extra_feat <- dim_args[intersect(names(dim_args), feat_formals)]
    extra_feat <- extra_feat[setdiff(
      names(extra_feat),
      c("srt", "features", "assay", "layer", "title", "legend.title")
    )]
    act_args[names(extra_feat)] <- extra_feat
    panels <- list(do.call(FeatureDimPlot, act_args))

    tf <- scenic_tf_from_regulon(regulon)
    expr_feature <- scenic_match_feature(tf, expr_names)
    if (isTRUE(compare_expression) && !is.na(expr_feature)) {
      exp_args <- act_args
      exp_args$features <- expr_feature
      exp_args$assay <- expression_assay
      exp_args$layer <- expression_layer
      exp_args$palette <- dim_args$expression_palette %||% "Spectral"
      exp_args$palcolor <- dim_args$expression_palcolor
      exp_args$bg_cutoff <- dim_args$expression_bg_cutoff %||% 0
      exp_args$title <- paste(expr_feature, "expression")
      exp_args$legend.title <- "Expression"
      panels <- c(panels, list(do.call(FeatureDimPlot, exp_args)))
    } else if (isTRUE(compare_expression) && is.na(expr_feature)) {
      log_message(
        "No matching expression feature for TF {.val {tf}} in assay {.val {expression_assay}}",
        message_type = "warning",
        verbose = verbose
      )
    }

    if (!is.null(region_auc) && regulon %in% rownames(region_auc)) {
      region_assay <- paste0(assay, "_regions")
      srt_region <- scenic_attach_auc_assay(
        srt = srt_use,
        auc_mat = region_auc,
        assay = region_assay
      )
      region_features <- tryCatch(
        rownames(srt_region[[region_assay]]),
        error = function(...) character()
      )
      region_args <- act_args
      region_args$srt <- srt_region
      region_args$features <- scenic_assay_feature(regulon, region_features)
      region_args$assay <- region_assay
      region_args$title <- paste(regulon, "region AUC")
      region_args$legend.title <- "Region AUC"
      panels <- c(panels, list(do.call(FeatureDimPlot, region_args)))
    }
    if (length(panels) == 1L) {
      return(panels[[1L]])
    }
    patchwork::wrap_plots(panels, ncol = length(panels))
  })

  plot <- patchwork::wrap_plots(c(list(p_cell), tf_plots), ncol = 1)
  if (!is.null(title)) {
    plot <- plot + patchwork::plot_annotation(title = title)
  }
  list(
    plot = plot,
    plots = c(list(p_cell), tf_plots),
    data = data.frame(regulon = regulons, TF = scenic_tf_from_regulon(regulons))
  )
}

scenic_plot_activity_cor_dumbbell <- function(
  srt,
  auc_mat,
  top_table,
  features = NULL,
  cor.features,
  cor.feature.labels = NULL,
  cor_method = "spearman",
  p_cutoff = 0.05,
  cor_sort.by = "difference",
  cor_cols = c("#56B4E9", "#E83F6F"),
  cor_xlim = NULL,
  cor_label = TRUE,
  assay = "scenic",
  layer = "data",
  title = NULL,
  point_size = 2,
  label_size = 3,
  verbose = TRUE
) {
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = rownames(auc_mat),
    top_table = top_table
  )
  target_mat <- scenic_resolve_cor_targets(
    srt = srt,
    features = cor.features,
    labels = cor.feature.labels,
    cells = colnames(auc_mat),
    assay = assay,
    layer = layer,
    verbose = verbose
  )
  target_levels <- colnames(target_mat)
  target_features <- attr(target_mat, "features")
  target_sources <- attr(target_mat, "sources")

  plot_data <- do.call(
    rbind,
    lapply(regulons, function(regulon) {
      do.call(
        rbind,
        lapply(seq_along(target_levels), function(target_idx) {
          x <- as.numeric(auc_mat[regulon, ])
          y <- as.numeric(target_mat[, target_idx])
          keep <- is.finite(x) & is.finite(y)
          cor_test <- NULL
          if (sum(keep) >= 3L &&
            stats::var(x[keep]) > 0 &&
            stats::var(y[keep]) > 0) {
            cor_args <- list(
              x = x[keep],
              y = y[keep],
              method = cor_method
            )
            if (cor_method %in% c("spearman", "kendall")) {
              cor_args$exact <- FALSE
            }
            cor_test <- tryCatch(
              do.call(stats::cor.test, cor_args),
              error = function(e) NULL
            )
          }
          cor_value <- if (is.null(cor_test)) {
            NA_real_
          } else {
            as.numeric(cor_test$estimate)
          }
          p_val <- if (is.null(cor_test)) {
            NA_real_
          } else {
            as.numeric(cor_test$p.value)
          }
          data.frame(
            regulon = regulon,
            TF = scenic_tf_from_regulon(regulon),
            target = target_levels[[target_idx]],
            target_feature = target_features[[target_idx]],
            target_source = target_sources[[target_idx]],
            cor = cor_value,
            p_val = p_val,
            n = sum(keep),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  rownames(plot_data) <- NULL
  plot_data[["significant"]] <- is.finite(plot_data[["p_val"]]) &
    plot_data[["p_val"]] < p_cutoff
  if (!any(is.finite(plot_data[["cor"]]))) {
    log_message(
      "No finite correlations remain for {.arg plot_type = 'activity_cor_dumbbell'}",
      message_type = "error"
    )
  }

  segment_data <- do.call(
    rbind,
    lapply(seq_along(regulons), function(regulon_idx) {
      regulon <- regulons[[regulon_idx]]
      one_df <- plot_data[plot_data[["regulon"]] == regulon, , drop = FALSE]
      data.frame(
        regulon = regulon,
        TF = scenic_tf_from_regulon(regulon),
        input_rank = regulon_idx,
        cor_1 = one_df[one_df[["target"]] == target_levels[[1]], "cor"],
        cor_2 = one_df[one_df[["target"]] == target_levels[[2]], "cor"],
        stringsAsFactors = FALSE
      )
    })
  )
  segment_data[["difference"]] <- abs(segment_data[["cor_1"]] - segment_data[["cor_2"]])
  segment_data[["max_abs"]] <- pmax(
    abs(segment_data[["cor_1"]]),
    abs(segment_data[["cor_2"]]),
    na.rm = TRUE
  )
  segment_data[["difference"]][!is.finite(segment_data[["difference"]])] <- -Inf
  segment_data[["max_abs"]][!is.finite(segment_data[["max_abs"]])] <- -Inf
  regulon_order <- switch(cor_sort.by,
    difference = segment_data[order(-segment_data[["difference"]], segment_data[["input_rank"]]), "regulon"],
    max_abs = segment_data[order(-segment_data[["max_abs"]], segment_data[["input_rank"]]), "regulon"],
    input = segment_data[order(segment_data[["input_rank"]]), "regulon"]
  )
  plot_data[["regulon"]] <- factor(plot_data[["regulon"]], levels = rev(regulon_order))
  segment_data[["regulon"]] <- factor(segment_data[["regulon"]], levels = rev(regulon_order))
  plot_data[["target"]] <- factor(plot_data[["target"]], levels = target_levels)

  sig_label <- paste0("p < ", signif(p_cutoff, 3))
  ns_label <- paste0("p >= ", signif(p_cutoff, 3))
  plot_data[["significance"]] <- ifelse(
    plot_data[["significant"]],
    sig_label,
    ns_label
  )
  plot_data[["significance"]] <- factor(
    plot_data[["significance"]],
    levels = c(ns_label, sig_label)
  )
  plot_data[["cor_label"]] <- ifelse(
    is.finite(plot_data[["cor"]]),
    sprintf("%.2f", plot_data[["cor"]]),
    ""
  )
  cor_cols <- stats::setNames(cor_cols, target_levels)
  shape_values <- stats::setNames(c(1, 21), c(ns_label, sig_label))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["cor"]], y = .data[["regulon"]])) +
    ggplot2::geom_segment(
      data = segment_data,
      ggplot2::aes(
        x = .data[["cor_1"]],
        xend = .data[["cor_2"]],
        y = .data[["regulon"]],
        yend = .data[["regulon"]]
      ),
      inherit.aes = FALSE,
      color = "grey80",
      linewidth = 0.45,
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = .data[["target"]],
        fill = .data[["target"]],
        shape = .data[["significance"]]
      ),
      size = max(point_size * 2.2, 4.5),
      stroke = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = cor_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cor_cols, drop = FALSE, guide = "none") +
    ggplot2::scale_shape_manual(values = shape_values, drop = FALSE) +
    theme_scop() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 12),
      axis.text = ggplot2::element_text(colour = "black", size = 10),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey86", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = "Correlation",
        override.aes = list(shape = 21, size = 4)
      ),
      shape = ggplot2::guide_legend(
        title = "Significance",
        override.aes = list(fill = c("white", "black"), color = "black")
      )
    ) +
    ggplot2::labs(
      x = "Correlation coefficient",
      y = NULL,
      title = title %||% "SCENIC regulon activity correlation"
    )
  if (!is.null(cor_xlim)) {
    p <- p + ggplot2::scale_x_continuous(limits = cor_xlim, n.breaks = 5)
  } else {
    p <- p + ggplot2::scale_x_continuous(n.breaks = 5)
  }
  if (isTRUE(cor_label)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[["cor_label"]]),
        size = max(label_size * 0.6, 2),
        color = "black",
        na.rm = TRUE,
        show.legend = FALSE
      )
  }

  list(plot = p, plots = list(p), data = plot_data)
}

scenic_resolve_cor_targets <- function(
  srt,
  features,
  labels = NULL,
  cells,
  assay = "scenic",
  layer = "data",
  verbose = TRUE
) {
  features <- as.character(features)
  labels <- labels %||% features
  labels <- as.character(labels)
  if (anyDuplicated(labels)) {
    log_message(
      "{.arg cor.feature.labels} must be unique",
      message_type = "error"
    )
  }
  target_mat <- matrix(
    NA_real_,
    nrow = length(cells),
    ncol = length(features),
    dimnames = list(cells, labels)
  )
  target_sources <- character(length(features))

  for (feature_idx in seq_along(features)) {
    feature <- features[[feature_idx]]
    if (feature %in% colnames(srt@meta.data)) {
      raw_value <- srt@meta.data[cells, feature, drop = TRUE]
      target_sources[[feature_idx]] <- "metadata"
    } else {
      raw_value <- NULL
      assay_candidates <- unique(c(SeuratObject::DefaultAssay(srt), assay))
      assay_candidates <- assay_candidates[
        assay_candidates %in% SeuratObject::Assays(srt)
      ]
      for (assay_use in assay_candidates) {
        if (!feature %in% rownames(srt@assays[[assay_use]])) {
          next
        }
        assay_mat <- tryCatch(
          GetAssayData5(srt, assay = assay_use, layer = layer),
          error = function(e) NULL
        )
        if (is.null(assay_mat)) {
          next
        }
        raw_value <- assay_mat[feature, cells, drop = TRUE]
        target_sources[[feature_idx]] <- paste0("assay:", assay_use)
        break
      }
      if (is.null(raw_value)) {
        log_message(
          "Cannot find {.val {feature}} in {.arg srt} metadata or selected assay features",
          message_type = "error"
        )
      }
    }
    if (is.factor(raw_value)) {
      raw_value <- as.character(raw_value)
    }
    value <- suppressWarnings(as.numeric(raw_value))
    if (!any(is.finite(value))) {
      log_message(
        "{.val {feature}} must contain finite numeric values for correlation",
        message_type = "error"
      )
    }
    target_mat[, feature_idx] <- value
  }
  attr(target_mat, "features") <- features
  attr(target_mat, "sources") <- target_sources
  log_message(
    "Resolved correlation targets: {.val {features}}",
    verbose = verbose
  )
  target_mat
}

scenic_plot_regulon_size <- function(
  srt,
  tool_name,
  top_table,
  features = NULL,
  top_n = 12,
  palette = "RdYlBu",
  palcolor = NULL,
  title = NULL
) {
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
  fill_colors <- scenic_continuous_colors(palette = palette, palcolor = palcolor)

  plot <- ggplot2::ggplot(
    size_data,
    ggplot2::aes(x = .data[["target_count"]], y = .data[["regulon"]], fill = .data[["target_count"]])
  ) +
    ggplot2::geom_col(width = 0.75) +
    scenic_fill_scale(fill_colors, "Targets") +
    theme_scop() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 12),
      axis.text = ggplot2::element_text(colour = "black", size = 10),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Target genes", y = NULL, title = title %||% "SCENIC regulon size")

  list(plot = plot, plots = list(plot), data = size_data)
}

scenic_plot_target_bar <- function(
  srt,
  tool_name,
  top_table,
  features = NULL,
  max_targets = 30,
  palette = "RdYlBu",
  palcolor = NULL,
  title = NULL
) {
  regulon_list <- scenic_get_regulon_list(srt, tool_name)
  regulons <- scenic_resolve_regulon_features(
    features = features,
    available = names(regulon_list),
    top_table = top_table,
    max_features = if (is.null(features)) 4 else NULL
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
  value_col <- if (all(is.na(plot_data[["importance"]]))) "rank_score" else "importance"
  if (identical(value_col, "rank_score")) {
    plot_data[["rank_score"]] <- ave(
      seq_len(nrow(plot_data)),
      plot_data[["regulon"]],
      FUN = function(idx) rev(seq_along(idx))
    )
  }
  fill_colors <- scenic_continuous_colors(palette = palette, palcolor = palcolor)
  plot_data[["regulon"]] <- factor(plot_data[["regulon"]], levels = unique(plot_data[["regulon"]]))
  plots <- lapply(split(plot_data, plot_data[["regulon"]], drop = TRUE), function(df) {
    df[["target"]] <- factor(df[["target"]], levels = rev(unique(df[["target"]])))
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[value_col]], y = .data[["target"]], fill = .data[[value_col]])
    ) +
      ggplot2::geom_col(width = 0.75) +
      scenic_fill_scale(
        fill_colors,
        if (identical(value_col, "importance")) "Importance" else "Rank"
      ) +
      theme_scop() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(colour = "black", size = 12),
        axis.text = ggplot2::element_text(colour = "black", size = 10),
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = if (identical(value_col, "importance")) "GRN importance" else "Target rank",
        y = NULL,
        title = unique(as.character(df[["regulon"]]))
      )
  })
  plot <- scenic_combine_plots(
    plots,
    combine = TRUE,
    ncol = min(3L, length(plots)),
    title = title %||% "SCENIC regulon targets"
  )

  list(plot = plot, plots = plots, data = plot_data)
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
  mat_sparse <- methods::as(as.matrix(mat), "dgCMatrix")
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
  check_r("matrixStats", verbose = FALSE)
  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- matrixStats::rowSds(mat, na.rm = TRUE)
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
  out <- as.character(regulon)
  out <- sub("\\([+-]\\)$", "", out)
  out <- sub("_direct_.*$", "", out)
  out <- sub("_[+-]_[+-]$", "", out)
  out <- sub("_[+-]$", "", out)
  out
}

scenic_regulon_feature_candidates <- function(feature) {
  feature <- as.character(feature)
  base <- scenic_tf_from_regulon(feature)
  unique(c(feature, base, paste0(base, "(+)"), paste0(base, "(-)")))
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
  available <- as.character(available)
  available_tf <- scenic_tf_from_regulon(available)
  matches <- lapply(
    as.character(features),
    function(feature) {
      candidates <- scenic_regulon_feature_candidates(feature)
      hit <- available[available %in% candidates]
      if (length(hit) == 0L) {
        feature_tf <- scenic_tf_from_regulon(feature)
        hit <- available[available_tf == feature | available_tf == feature_tf]
      }
      unique(hit)
    }
  )
  unique(unlist(matches, use.names = FALSE))
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
  matched_all <- lapply(
    as.character(features),
    function(feature) {
      candidates <- scenic_regulon_feature_candidates(feature)
      hit <- candidates[candidates %in% available]
      if (length(hit) == 0) {
        return(character(0))
      }
      hit
    }
  )

  out <- stats::setNames(rep(NA_character_, length(regulons)), regulons)
  for (idx in seq_along(matched_all)) {
    for (regulon in matched_all[[idx]]) {
      if (regulon %in% names(out) && is.na(out[[regulon]])) {
        out[[regulon]] <- split_values[[idx]]
      }
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

scenic_top_regulon_source <- function(top_table, regulons) {
  if (is.null(top_table) || nrow(top_table) == 0L || length(regulons) == 0L) {
    return(NULL)
  }
  top_table <- top_table[top_table[["regulon"]] %in% regulons, , drop = FALSE]
  if (nrow(top_table) == 0L) {
    return(NULL)
  }
  top_table <- top_table[!duplicated(top_table[["regulon"]]), , drop = FALSE]
  top_table <- top_table[match(regulons, top_table[["regulon"]]), , drop = FALSE]
  top_table <- top_table[!is.na(top_table[["regulon"]]), , drop = FALSE]
  if (nrow(top_table) == 0L) {
    return(NULL)
  }
  source_group <- factor(
    top_table[["group"]],
    levels = unique(as.character(top_table[["group"]]))
  )
  list(
    regulon = top_table[["regulon"]],
    group = stats::setNames(source_group, top_table[["regulon"]]),
    rank = stats::setNames(top_table[["rank"]], top_table[["regulon"]])
  )
}

scenic_group_row_maxima <- function(mat) {
  check_r("matrixStats", verbose = FALSE)
  has_value <- rowSums(!is.na(mat)) > 0L
  mat_for_max <- mat
  mat_for_max[is.na(mat_for_max)] <- -Inf
  max_group_idx <- max.col(mat_for_max, ties.method = "first")
  max_group_idx[!has_value] <- NA_integer_
  max_value <- matrixStats::rowMaxs(mat_for_max)
  max_value[!has_value] <- NA_real_
  names(max_group_idx) <- rownames(mat)
  names(max_value) <- rownames(mat)
  list(group_idx = max_group_idx, value = max_value)
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
  maxima <- scenic_group_row_maxima(mat)
  max_group_idx <- maxima[["group_idx"]]
  max_value <- maxima[["value"]]
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
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = nrow(auc_mat), dimnames = list(rownames(auc_mat), group_names))
  }
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
  result <- srt@tools[[tool_name]]
  adjacency <- result[["adjacency"]]
  if (is.null(adjacency) || nrow(adjacency) == 0) {
    adjacency <- result[["tf_gene"]]
  }
  if ((is.null(adjacency) || nrow(adjacency) == 0) && !is.null(result[["triplets"]])) {
    triplets <- as.data.frame(result[["triplets"]], stringsAsFactors = FALSE)
    if (nrow(triplets) > 0 && all(c("TF", "gene") %in% colnames(triplets))) {
      score_col <- scenic_first(
        intersect(c("score", "importance", "weight"), colnames(triplets)),
        NULL
      )
      adjacency <- unique(triplets[, c("TF", "gene", score_col), drop = FALSE])
      colnames(adjacency)[colnames(adjacency) == "gene"] <- "target"
    }
  }
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
  result <- srt@tools[[tool_name]]
  if (!is.null(result[["scores_cells_by_regulon"]])) {
    scores_cells_by_regulon <- result[["scores_cells_by_regulon"]]
    auc_mat <- as.matrix(scores_cells_by_regulon)
    if (is.null(rownames(auc_mat)) || is.null(colnames(auc_mat))) {
      log_message(
        "{.arg scores_cells_by_regulon} must have cell row names and regulon column names",
        message_type = "error"
      )
    }
    return(t(auc_mat))
  }
  if (!is.null(result[["scores"]])) {
    scores <- as.matrix(result[["scores"]])
    if (!is.null(rownames(scores)) && !is.null(colnames(scores))) {
      if (all(colnames(scores) %in% colnames(srt))) {
        return(scores)
      }
      if (all(rownames(scores) %in% colnames(srt))) {
        return(t(scores))
      }
    }
  }
  if (!is.null(result[["auc"]])) {
    auc <- as.matrix(result[["auc"]])
    if (!is.null(rownames(auc)) && !is.null(colnames(auc))) {
      if (all(rownames(auc) %in% colnames(srt))) {
        return(t(auc))
      }
      if (all(colnames(auc) %in% colnames(srt))) {
        return(auc)
      }
    }
  }

  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "Cannot find SCENIC results in tools slot {.val {tool_name}} or assay {.val {assay}}",
      message_type = "error"
    )
  }
  auc_mat <- GetAssayData5(srt, assay = assay, layer = layer)
  if (nrow(auc_mat) == 0L && !identical(layer, "counts")) {
    auc_mat <- GetAssayData5(srt, assay = assay, layer = "counts")
  }
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
  1 - sqrt(pmax(jsd, 0))
}

scenic_calc_jsd <- function(p_regulon, p_cell_type) {
  scenic_entropy((p_regulon + p_cell_type) / 2) -
    ((scenic_entropy(p_regulon) + scenic_entropy(p_cell_type)) / 2)
}

scenic_entropy <- function(p_vector) {
  p_vector <- p_vector[p_vector > 0]
  -sum(p_vector * log2(p_vector))
}

scenic_continuous_colors <- function(palette = "RdYlBu", palcolor = NULL) {
  palette_colors(
    type = "continuous",
    palette = palette,
    palcolor = palcolor
  )
}

scenic_color_scale <- function(colors, name) {
  ggplot2::scale_color_gradientn(
    colours = colors,
    name = name,
    na.value = "grey80",
    guide = ggplot2::guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      title.hjust = 0
    )
  )
}

scenic_fill_scale <- function(colors, name) {
  ggplot2::scale_fill_gradientn(
    colours = colors,
    name = name,
    na.value = "grey80",
    guide = ggplot2::guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      title.hjust = 0
    )
  )
}

scenic_expression_assay <- function(
  srt,
  expression_assay = NULL,
  activity_assay = NULL
) {
  assays <- SeuratObject::Assays(srt)
  if (!is.null(expression_assay)) {
    if (!expression_assay %in% assays) {
      log_message(
        "{.arg expression_assay} {.val {expression_assay}} is not present in {.arg srt}",
        message_type = "error"
      )
    }
    return(expression_assay)
  }
  skip <- unique(c(activity_assay, "scenic", "scenicplus", "dorothea"))
  if ("RNA" %in% assays) {
    return("RNA")
  }
  remaining <- setdiff(assays, skip)
  if (length(remaining) > 0L) {
    return(remaining[[1L]])
  }
  SeuratObject::DefaultAssay(srt)
}

scenic_match_feature <- function(feature, rownames_vec) {
  if (length(rownames_vec) == 0L || is.na(feature) || !nzchar(feature)) {
    return(NA_character_)
  }
  if (feature %in% rownames_vec) {
    return(feature)
  }
  hit <- rownames_vec[tolower(rownames_vec) == tolower(feature)]
  if (length(hit) >= 1L) {
    return(hit[[1L]])
  }
  NA_character_
}

scenic_get_expression_matrix <- function(srt, assay, layer = "data") {
  out <- tryCatch(
    suppressWarnings(GetAssayData5(srt, assay = assay, layer = layer)),
    error = function(...) NULL
  )
  if (is.null(out) || nrow(out) == 0L || ncol(out) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }
  as.matrix(out)
}

scenic_resolve_expression_layer <- function(srt, assay, layer = "data") {
  expr <- scenic_get_expression_matrix(srt, assay = assay, layer = layer)
  empty <- nrow(expr) == 0L ||
    ncol(expr) == 0L ||
    !any(is.finite(expr) & expr != 0)
  if (isTRUE(empty) && !identical(layer, "counts")) {
    counts <- scenic_get_expression_matrix(srt, assay = assay, layer = "counts")
    if (nrow(counts) > 0L && any(is.finite(counts) & counts != 0)) {
      return("counts")
    }
  }
  layer
}

scenic_assay_feature <- function(feature, assay_features) {
  hit <- scenic_match_feature(feature, assay_features)
  if (!is.na(hit)) {
    return(hit)
  }
  scenic_match_feature(
    gsub("_", "-", as.character(feature), fixed = TRUE),
    assay_features
  )
}

scenic_assay_features <- function(features, assay_features) {
  vapply(
    as.character(features),
    scenic_assay_feature,
    character(1),
    assay_features = assay_features
  )
}

scenic_group_tf_expression <- function(
  srt,
  tfs,
  group_annotation,
  group_names,
  expression_assay = NULL,
  expression_layer = "data",
  activity_assay = "scenic",
  scale = TRUE,
  verbose = TRUE
) {
  expression_assay <- tryCatch(
    scenic_expression_assay(
      srt = srt,
      expression_assay = expression_assay,
      activity_assay = activity_assay
    ),
    error = function(...) NULL
  )
  if (is.null(expression_assay)) {
    return(NULL)
  }
  expression_layer <- scenic_resolve_expression_layer(
    srt = srt,
    assay = expression_assay,
    layer = expression_layer
  )
  expr <- scenic_get_expression_matrix(
    srt = srt,
    assay = expression_assay,
    layer = expression_layer
  )
  if (nrow(expr) == 0L) {
    return(NULL)
  }
  matched <- vapply(tfs, scenic_match_feature, character(1), rownames_vec = rownames(expr))
  keep <- !is.na(matched)
  if (!any(keep)) {
    return(NULL)
  }
  if (any(!keep)) {
    log_message(
      "Dropping {.val {sum(!keep)}} TFs missing from assay {.val {expression_assay}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  expr <- expr[matched[keep], , drop = FALSE]
  rownames(expr) <- tfs[keep]
  cells <- intersect(colnames(expr), names(group_annotation))
  if (length(cells) == 0L) {
    return(NULL)
  }
  expr <- expr[, cells, drop = FALSE]
  avg <- scenic_group_average_matrix(
    auc_mat = expr,
    group_annotation = group_annotation[cells],
    group_names = group_names
  )
  if (isTRUE(scale)) {
    avg <- scenic_scale_rows(avg)
  }
  list(matrix = avg, assay = expression_assay, tfs = tfs[keep])
}

scenic_get_region_auc_matrix <- function(srt, tool_name) {
  result <- srt@tools[[tool_name]]
  if (is.null(result)) {
    return(NULL)
  }
  for (nm in c("auc_regions", "region_auc", "auc_region")) {
    mat <- result[[nm]]
    if (is.null(mat)) {
      next
    }
    mat <- as.matrix(mat)
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      next
    }
    if (all(colnames(mat) %in% colnames(srt))) {
      return(mat)
    }
    if (all(rownames(mat) %in% colnames(srt))) {
      return(t(mat))
    }
  }
  NULL
}
