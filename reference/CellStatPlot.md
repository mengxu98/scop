# Statistical plot of cells

Statistical plot of cells

## Usage

``` r
CellStatPlot(
  srt,
  stat.by,
  group.by = NULL,
  split.by = NULL,
  bg.by = NULL,
  cells = NULL,
  flip = FALSE,
  NA_color = "grey",
  NA_stat = TRUE,
  keep_empty = FALSE,
  individual = FALSE,
  stat_level = NULL,
  plot_type = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord",
    "venn", "upset"),
  stat_type = c("percent", "count"),
  position = c("stack", "dodge"),
  palette = "Chinese",
  palcolor = NULL,
  alpha = 1,
  bg_palette = "Chinese",
  bg_palcolor = NULL,
  bg_alpha = 0.2,
  label = FALSE,
  label.size = 3.5,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- stat.by:

  A character vector specifying the features to plot.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- bg.by:

  A character vector specifying the variable to use as the background
  color. Default is `NULL`.

- cells:

  A character vector of cell names to use. Default is `NULL`.

- flip:

  A logical specifying whether to flip the plot vertically. Default is
  `FALSE`.

- NA_color:

  The color to use for missing values.

- NA_stat:

  Whether to include missing values in the plot. Default is `TRUE`.

- keep_empty:

  Whether to keep empty levels in the plot. Default is `FALSE`.

- individual:

  Whether to create individual plots for each group. Default is `FALSE`.

- stat_level:

  The level(s) of the variable(s) specified in `stat.by` to include in
  the plot. Default is `NULL`.

- plot_type:

  A string specifying the type of plot to create. Possible values are
  `"violin"`, `"box"`, `"bar"`, `"dot"`, or `"col"`. Default is
  `"violin"`.

- stat_type:

  The type of statistic to compute for the plot. Can be one of
  `"percent"` or `"count"`.

- position:

  The position adjustment for the plot. Can be one of `"stack"` or
  `"dodge"`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- alpha:

  The transparency of the plot. Default is `1`.

- bg_palette:

  A string specifying the color palette to use for the background.
  Default is `"Chinese"`.

- bg_palcolor:

  A character vector specifying specific colors to use for the
  background. Default is `NULL`.

- bg_alpha:

  The transparency of the background. Default is `0.2`.

- label:

  Whether to add labels on the plot. Default is `FALSE`.

- label.size:

  The size of the labels.

- label.fg:

  The foreground color of the labels.

- label.bg:

  The background color of the labels.

- label.bg.r:

  The radius of the rounded corners of the label background.

- aspect.ratio:

  Aspect ratio of the panel. Default is `NULL`.

- title:

  The text for the title. Default is `NULL`.

- subtitle:

  The text for the subtitle for the plot which will be displayed below
  the title. Default is `NULL`.

- xlab:

  The x-axis label of the plot. Default is `NULL`.

- ylab:

  A string specifying the label of the y-axis. Default is
  `"Expression level"`.

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

- grid_major:

  Whether to show major panel grid lines. Default is `TRUE`.

- grid_major_colour:

  Color of major panel grid lines.

- grid_major_linetype:

  Linetype of major panel grid lines.

- grid_major_linewidth:

  Line width of major panel grid lines.

- combine:

  Combine plots into a single `patchwork` object. If `FALSE`, return a
  list of ggplot objects.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

- force:

  Whether to force drawing regardless of maximum levels in any cell
  group is greater than 100. Default is `FALSE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## See also

[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-26 00:57:07] Start standard processing workflow...
#> ℹ [2026-04-26 00:57:08] Checking a list of <Seurat>...
#> ! [2026-04-26 00:57:08] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-26 00:57:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-26 00:57:10] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-26 00:57:10] Use the separate HVF from `srt_list`
#> ℹ [2026-04-26 00:57:10] Number of available HVF: 2000
#> ℹ [2026-04-26 00:57:11] Finished check
#> ℹ [2026-04-26 00:57:11] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-26 00:57:11] Perform pca linear dimension reduction
#> ℹ [2026-04-26 00:57:12] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-26 00:57:12] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-26 00:57:12] Reorder clusters...
#> ℹ [2026-04-26 00:57:13] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-26 00:57:13] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-26 00:57:13] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-26 00:57:16] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-26 00:57:20] Standard processing workflow completed
p1 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  label = TRUE
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)
p1
#> Error: object 'p1' not found

thisplot::panel_fix(
  p1, height = 2, width = 3
)
#> Error: object 'p1' not found

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  stat_type = "count",
  position = "dodge",
  label = TRUE
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  bg.by = "CellType",
  palette = "Set1",
  stat_type = "count",
  position = "dodge"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "bar"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "rose"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "ring"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "pie"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "dot"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "rose"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "ring"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "area"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "dot"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "trend"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar",
  individual = TRUE
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "bar"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "area"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "dot"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "trend"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "bar",
  position = "dodge",
  label = TRUE
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose",
  position = "dodge",
  label = TRUE
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring",
  position = "dodge",
  label = TRUE
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "sankey"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "chord"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "venn",
  stat_level = list(
    CellType = c("Ductal", "Ngn3-low-EP"),
    Phase = "S"
  )
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

pancreas_sub$Progenitor <- pancreas_sub$CellType %in% c("Ngn3-low-EP", "Ngn3-high-EP")
pancreas_sub$G2M <- pancreas_sub$Phase == "G2M"
pancreas_sub$Fancb_Expressed <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)["Fancb", ] > 0
pancreas_sub$Dlg3_Expressed <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)["Dlg3", ] > 0
CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
  ),
  plot_type = "venn",
  stat_level = "TRUE"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
  ),
  plot_type = "upset",
  stat_level = "TRUE"
)
#> Error in StatPlot(meta_data, stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, flip = flip, NA_color = NA_color, NA_stat = NA_stat,     keep_empty = keep_empty, individual = individual, stat_level = stat_level,     plot_type = plot_type, stat_type = stat_type, position = position,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, label = label,     label.size = label.size, label.fg = label.fg, label.bg = label.bg,     label.bg.r = label.bg.r, aspect.ratio = aspect.ratio, title = title,     subtitle = subtitle, xlab = xlab, ylab = ylab, legend.position = legend.position,     legend.direction = legend.direction, theme_use = theme_use,     theme_args = theme_args, grid_major = grid_major, grid_major_colour = grid_major_colour,     grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth,     combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,     force = force, seed = seed): unused arguments (grid_major = grid_major, grid_major_colour = grid_major_colour, grid_major_linetype = grid_major_linetype, grid_major_linewidth = grid_major_linewidth)

sum(
  pancreas_sub$Progenitor == "FALSE" &
    pancreas_sub$G2M == "FALSE" &
    pancreas_sub$Fancb_Expressed == "TRUE" &
    pancreas_sub$Dlg3_Expressed == "FALSE"
)
#> [1] 6
```
