#' @title Plot SCENIC+ eRegulon activity
#'
#' @description
#' [SCENICPlot()] with `tool_name = "SCENICPlus"` and `assay = "scenicplus"`.
#' Default `"heatmap_dotplot"`: color = TF expression, size = RSS.
#'
#' @md
#' @inheritParams SCENICPlot
#' @param srt A Seurat object from [RunSCENICPlus()].
#' @param tool_name Tools slot. Default `"SCENICPlus"`.
#' @param assay eRegulon AUC assay. Default `"scenicplus"`.
#' @param plot_type Plot type. Default `"heatmap_dotplot"`.
#' @param ... Passed to [SCENICPlot()].
#'
#' @return Same as [SCENICPlot()].
#'
#' @seealso [RunSCENICPlus], [SCENICPlot]
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
#' scenicplus_dot <- SCENICPlusPlot(pancreas_sub, group.by = "CellType")
#' example_tfs <- unique(scenicplus_dot$top_table$TF)[1:3]
#' SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "egrn", features = example_tfs)
#' SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "network_graph", features = example_tfs)
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
    "egrn",
    "overlap",
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
