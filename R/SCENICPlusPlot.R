#' @title Plot SCENIC+ eRegulon activity and specificity
#'
#' @description
#' A [SCENICPlot()] wrapper with SCENIC+ defaults (`tool_name = "SCENICPlus"`,
#' `assay = "scenicplus"`). The default `"heatmap_dotplot"` matches the
#' official SCENIC+ signature figure: color is TF expression and size is
#' regulon specificity score (RSS). `"eregulon_dim"` compares TF expression
#' with gene-based AUC, and region-based AUC when it is stored in
#' `srt@tools$SCENICPlus`. `"coverage"` draws accessibility and region–gene
#' tracks for a target locus.
#'
#' @md
#' @inheritParams SCENICPlot
#' @param srt A Seurat object containing results from [RunSCENICPlus()].
#' @param tool_name Name of the `srt@tools` entry. Default `"SCENICPlus"`.
#' @param assay Assay storing eRegulon AUC scores. Default `"scenicplus"`.
#' @param plot_type Plot type. Defaults to `"heatmap_dotplot"`. See
#' [SCENICPlot()] for the full set of shared views.
#' @param ... Additional arguments passed to [SCENICPlot()], including
#' `features`, `compare_expression`, heatmap options, and network options.
#'
#' @return The same object returned by [SCENICPlot()].
#'
#' @seealso [RunSCENICPlus], [SCENICPlot], [FeatureDimPlot], [FeatureStatPlot]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#' pancreas_sub <- RunSCENICPlus(
#'   pancreas_sub,
#'   backend = "python",
#'   python_result_dir = "test/scenicplus"
#' )
#'
#' scenicplus_dot <- SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType"
#' )
#' scenicplus_dot$plot
#' example_regulons <- unique(scenicplus_dot$top_table$regulon)[1:2]
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = example_regulons,
#'   assay = "scenicplus",
#'   reduction = "StandardpcaUMAP2D"
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = example_regulons,
#'   group.by = "CellType",
#'   assay = "scenicplus"
#' )
#'
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "rss_rank"
#' )
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "eregulon_dim",
#'   features = example_regulons
#' )
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "network",
#'   features = example_regulons
#' )
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "coverage",
#'   features = example_regulons
#' )
#' }
SCENICPlusPlot <- function(
  srt,
  group.by,
  tool_name = "SCENICPlus",
  assay = "scenicplus",
  plot_type = c(
    "heatmap_dotplot",
    "rss_rank",
    "rss_heatmap",
    "rss_dotplot",
    "activity_heatmap",
    "activity_violin",
    "activity_dim",
    "eregulon_dim",
    "activity_cor_dumbbell",
    "regulon_size",
    "network_graph",
    "network",
    "target_bar",
    "coverage"
  ),
  ...
) {
  plot_type <- match.arg(plot_type)
  extra <- list(...)
  args <- c(
    list(
      srt = srt,
      group.by = group.by,
      tool_name = tool_name,
      assay = assay,
      plot_type = plot_type
    ),
    extra
  )
  do.call(SCENICPlot, args)
}
