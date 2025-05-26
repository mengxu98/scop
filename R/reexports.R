#' @import ggplot2 BiocParallel Seurat
#' @importFrom ggplot2 %+replace%
#' @importFrom ComplexHeatmap %v%
#' @importFrom grDevices palette
#' @importFrom igraph layout_with_fr cluster_fast_greedy
#' @importFrom methods new
#' @importFrom stats median
#' @importFrom ggrepel GeomTextRepel
#' @importFrom Signac RunSVD
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull geom_mark_rect geom_mark_circle
#' @importFrom dplyr "%>%" %>% .data
#' @export
dplyr::`%>%`

#' @importFrom rlang %||%
#' @export
rlang::`%||%`

utils::globalVariables(
  c(
    ":=", "Axis_1", "Axis_2", "X", "Y",
    "block_graphics", "capture.output", "celltype", "cluster",
    "color",
    "colour", "combn", "count", "Count", "Database",
    "Description", "DescriptionP", "dim1", "dim2",
    "dx", "dy", "error", "fill",
    "Freq", "from_dim1", "from_dim2", "gene", "GeneName",
    "GeneRatio", "groupn", "Groups",
    "id", "ID", "intersection", "is_branch",
    "keyword1", "keyword2", "label", "label_color", "Lineages",
    "name",
    "nodes", "NES", "palcolor",
    "raw_Axis_1",
    "raw_Axis_2", "runningScore",
    "score", "setSize", "Significant", "size",
    "status", "step",
    "to_dim1", "to_dim2", "value", "values", "variable",
    "weight", "xend", "y", "yend",
    ".",
    ".data",
    "x",
    "node",
    "next_node",
    "next_x",
    "..r"
  )
)
