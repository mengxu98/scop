#' @import ggplot2 Seurat thisutils thisplot
#' @importFrom ggplot2 %+replace%
#' @importFrom ComplexHeatmap %v%
#' @importFrom grDevices palette
#' @importFrom igraph layout_with_fr cluster_fast_greedy
#' @importFrom methods new
#' @importFrom stats median
#' @importFrom ggrepel GeomTextRepel
#' @importFrom Signac RunSVD
#' @importFrom thisutils wilkinsonp maximump minimump meanp votep sump log_message
#' @importFrom thisplot theme_this
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull geom_mark_rect geom_mark_circle
#' @importFrom dplyr "%>%" %>% .data
#' @export
dplyr::`%>%`

#' @importFrom rlang %||%
#' @export
rlang::`%||%`

log_message <- thisutils::log_message
theme_scop <- thisplot::theme_this

utils::globalVariables(
  c(
    ":=", ".", ".data",
    "avg_log2FC", "Axis_1", "Axis_2",
    "block_graphics", "boot_CI_2.5", "boot_CI_97.5",
    "celltype", "cluster", "clusters", "color", "colour",
    "combn", "conda_info", "count", "Count", "Database",
    "Description", "DescriptionP", "dim1", "dim2",
    "dx", "dy", "error", "fill", "Freq",
    "from_dim1", "from_dim2", "gene", "GeneName", "GeneRatio",
    "group", "group1", "groupn", "Groups",
    "id", "ID", "intersection", "is_branch",
    "keyword1", "keyword2", "label", "label_color", "Lineages",
    "name", "NES", "next_node", "next_x", "node", "nodes",
    "obs_log2FD", "palcolor", "raw_Axis_1", "raw_Axis_2", "runningScore",
    "score", "setSize", "significance", "Significant", "size",
    "status", "step", "to_dim1", "to_dim2",
    "value", "values", "variable", "weight",
    "x", "x_angle", "x_num", "X", "xend", "x_plot",
    "y", "y_radius", "Y", "yend", "y_plot"
  )
)
