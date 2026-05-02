#' @import ggplot2 Seurat thisutils thisplot
#' @importFrom ComplexHeatmap %v%
#' @importFrom ggrepel GeomTextRepel
#' @importFrom Signac RunSVD
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull geom_mark_rect geom_mark_circle
#' @importFrom dplyr "%>%" %>% .data
#' @export
dplyr::`%>%`

#' @importFrom rlang %||%
#' @export
rlang::`%||%`

theme_scop <- thisplot::theme_this

#' OmicVerse-inspired ggplot theme
#'
#' A clean white theme that can be used through `theme_use = "theme_omicverse"`
#' in scop plotting functions.
#'
#' @param aspect.ratio Optional fixed panel aspect ratio.
#' @param base_size Base text size.
#' @param base_family Base font family.
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme.
#' @export
theme_omicverse <- function(
  aspect.ratio = NULL,
  base_size = 12,
  base_family = "",
  ...
) {
  extra_theme_args <- list(...)
  theme_obj <- ggplot2::theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(
        colour = "#ECECEC",
        linewidth = 0.25
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#333333", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "#333333", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "#202020"),
      axis.title = ggplot2::element_text(colour = "#202020"),
      strip.background = ggplot2::element_rect(fill = "#F4F4F4", colour = NA),
      strip.text = ggplot2::element_text(colour = "#202020", face = "bold"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(colour = "#555555", hjust = 0)
    )

  if (!is.null(aspect.ratio)) {
    theme_obj <- theme_obj + ggplot2::theme(aspect.ratio = aspect.ratio)
  }
  if (length(extra_theme_args) > 0L) {
    theme_obj <- theme_obj + do.call(
      ggplot2::theme,
      c(extra_theme_args, list(validate = FALSE))
    )
  }

  theme_obj
}

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
    "y", "y_radius", "Y", "yend", "y_plot",
    "cell", "contribution", "dataset", "edge_colour", "edge_id", "incoming", "incoming_diff",
    "interaction_label", "interaction_name", "interaction_plot", "label_plot", "ligand",
    "magnitude", "object_names", "outgoing", "outgoing_diff", "pair", "pathway", "prob",
    "receiver", "reg", "sender", "sender_display", "size_scaled", "specificity", "target",
    "total_links", "value_1", "value_2", "x_from", "x_from_lr", "x_from_rr", "x_from_sl",
    "x_label", "x_loop_from", "x_loop_to", "x_to", "x_to_lr", "x_to_rr", "x_to_sl", "y_from",
    "y_from_lr", "y_from_rr", "y_from_sl", "y_label", "y_loop_from", "y_loop_to", "y_to",
    "y_to_lr", "y_to_rr", "y_to_sl",
    "credible_label", "direction", "effect_use", "enrich_key", "hdi_2.5", "hdi_97.5",
    "idx", "label_id", "minus_log10", "xmax", "xmin"
  )
)
