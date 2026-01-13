# Graph Plot

A function to plot a graph with nodes and edges.

## Usage

``` r
GraphPlot(
  node,
  edge,
  transition = NULL,
  node_coord = c("x", "y"),
  node_group = NULL,
  node_palette = "Paired",
  node_palcolor = NULL,
  node_size = 4,
  node_alpha = 1,
  node_highlight = NULL,
  node_highlight_color = "red",
  label = FALSE,
  label.size = 3.5,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1,
  label_insitu = FALSE,
  label_repel = FALSE,
  label_repulsion = 20,
  label_point_size = 1,
  label_point_color = "black",
  label_segment_color = "black",
  edge_threshold = 0.01,
  use_triangular = c("upper", "lower", "both"),
  edge_line = c("straight", "curved"),
  edge_line_curvature = 0.3,
  edge_line_angle = 90,
  edge_color = "grey40",
  edge_size = c(0.2, 1),
  edge_alpha = 0.5,
  edge_shorten = 0,
  edge_offset = 0,
  edge_highlight = NULL,
  edge_highlight_color = "red",
  transition_threshold = 0.01,
  transition_line = c("straight", "curved"),
  transition_line_curvature = 0.3,
  transition_line_angle = 90,
  transition_color = "black",
  transition_size = c(0.2, 1),
  transition_alpha = 1,
  transition_arrow_type = "closed",
  transition_arrow_angle = 20,
  transition_arrow_length = grid::unit(0.02, "npc"),
  transition_shorten = 0.05,
  transition_offset = 0,
  transition_highlight = NULL,
  transition_highlight_color = "red",
  aspect.ratio = 1,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE
)
```

## Arguments

- node:

  A data frame representing the nodes of the graph.

- edge:

  A matrix representing the edges of the graph.

- transition:

  A matrix representing the transitions between nodes.

- node_coord:

  A character vector specifying the names of the columns in `node` that
  represent the x and y coordinates.

- node_group:

  A character vector specifying the name of the column in `node` that
  represents the grouping of the nodes.

- node_palette:

  A character vector specifying the name of the color palette for node
  groups.

- node_palcolor:

  A character vector specifying the names of the colors for each node
  group.

- node_size:

  A numeric value or column name of `node` specifying the size of the
  nodes.

- node_alpha:

  A numeric value or column name of `node` specifying the transparency
  of the nodes.

- node_highlight:

  A character vector specifying the names of nodes to highlight.

- node_highlight_color:

  A character vector specifying the color for highlighting nodes.

- label:

  Whether to show labels for the nodes.

- label.size:

  The size of the labels.

- label.fg:

  A character vector specifying the foreground color of the labels.

- label.bg:

  A character vector specifying the background color of the labels.

- label.bg.r:

  The background color transparency of the labels.

- label_insitu:

  Whether to display the node group labels in situ or as numeric values.

- label_repel:

  Whether to use force-directed label repulsion.

- label_repulsion:

  The repulsion force for labels.

- label_point_size:

  The size of the label points.

- label_point_color:

  A character vector specifying the color of the label points.

- label_segment_color:

  A character vector specifying the color for the label segments.

- edge_threshold:

  The threshold for removing edges.

- use_triangular:

  A character vector specifying which part of the edge matrix to use
  (upper, lower, both).

- edge_line:

  A character vector specifying the type of line for edges (straight,
  curved).

- edge_line_curvature:

  The curvature of curved edges.

- edge_line_angle:

  The angle of curved edges.

- edge_color:

  A character vector specifying the color of the edges.

- edge_size:

  A numeric vector specifying the range of edge sizes.

- edge_alpha:

  The transparency of the edges.

- edge_shorten:

  The length of the edge shorten.

- edge_offset:

  The length of the edge offset.

- edge_highlight:

  A character vector specifying the names of edges to highlight.

- edge_highlight_color:

  A character vector specifying the color for highlighting edges.

- transition_threshold:

  The threshold for removing transitions.

- transition_line:

  A character vector specifying the type of line for transitions
  (straight, curved).

- transition_line_curvature:

  The curvature of curved transitions.

- transition_line_angle:

  The angle of curved transitions.

- transition_color:

  A character vector specifying the color of the transitions.

- transition_size:

  A numeric vector specifying the range of transition sizes.

- transition_alpha:

  The transparency of the transitions.

- transition_arrow_type:

  A character vector specifying the type of arrow for transitions
  (closed, open).

- transition_arrow_angle:

  The angle of the transition arrow.

- transition_arrow_length:

  The length of the transition arrow.

- transition_shorten:

  The length of the transition shorten.

- transition_offset:

  The length of the transition offset.

- transition_highlight:

  A character vector specifying the names of transitions to highlight.

- transition_highlight_color:

  A character vector specifying the color for highlighting transitions.

- aspect.ratio:

  Aspect ratio of the panel. Default is `1`.

- title:

  The text for the title. Default is `NULL`.

- subtitle:

  The text for the subtitle for the plot which will be displayed below
  the title. Default is `NULL`.

- xlab:

  The x-axis label of the plot. Default is `NULL`.

- ylab:

  The y-axis label of the plot. Default is `NULL`.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

- legend.direction:

  The direction of the legend in the plot. Can be one of `"vertical"` or
  `"horizontal"`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- return_layer:

  Whether to return the plot layers as a list. Defaults is `FALSE`.

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
