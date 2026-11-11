#' @title Plot SCENIC+ eRegulon activity
#'
#' @description
#' [SCENICPlot()] with `tool_name = "SCENICPlus"` and `assay = "scenicplus"`.
#' Default `"heatmap_dotplot"`: color = TF expression, size = RSS.
#'
#' @md
#' @inheritParams SCENICPlot
#' @param srt Deprecated alias for `object`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param object A Seurat object from [RunSCENICPlus()].
#' @param tool_name Tools slot. Default `"SCENICPlus"`.
#' @param assay ERegulon AUC assay. Default `"scenicplus"`.
#' @param plot_type Plot type. Default `"heatmap_dotplot"`.
#' @param ... Arguments passed to other functions.
#' @return Same as [SCENICPlot()].
#'
#' @seealso [RunSCENICPlus], [SCENICPlot], [GroupHeatmap],
#' [FeatureStatPlot], [FeatureDimPlot]
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
#' regulons <- head(rownames(pancreas_sub[["scenicplus"]]), 6)
#' GroupHeatmap(
#'   pancreas_sub, group.by = "CellType", features = regulons,
#'   assay = "scenicplus", layer = "data", exp_method = "raw",
#'   lib_normalize = FALSE
#' )$plot
#' FeatureStatPlot(
#'   pancreas_sub, group.by = "CellType", stat.by = regulons,
#'   assay = "scenicplus", layer = "data", plot_type = "violin"
#' )
#' FeatureDimPlot(
#'   pancreas_sub, features = regulons, assay = "scenicplus",
#'   layer = "data", reduction = "UMAP"
#' )
#' SCENICPlusPlot(pancreas_sub, group.by = "CellType", plot_type = "egrn", features = example_tfs)
#' SCENICPlusPlot(
#'   pancreas_sub,
#'   group.by = "CellType", plot_type = "network_graph", features = example_tfs
#' )
#' }
SCENICPlusPlot <- function(
  object,
  group.by,
  tool_name = "SCENICPlus",
  assay = "scenicplus",
  plot_type = c("heatmap_dotplot", "rss_rank", "rss_heatmap", "rss_dotplot", "eregulon_dim", "activity_cor_dumbbell", "regulon_size", "network_graph", "network", "egrn", "overlap", "target_bar", "coverage"),
  srt = NULL,
  ...
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  plot_type <- match.arg(plot_type)
  extra <- list(...)
  args <- c(
    list(
      object = srt,
      group.by = group.by,
      tool_name = tool_name,
      assay = assay,
      plot_type = plot_type
    ),
    extra
  )
  do.call(SCENICPlot, args)
}
