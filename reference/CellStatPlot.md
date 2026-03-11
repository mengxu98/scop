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
#> ℹ [2026-03-11 16:43:00] Start standard scop workflow...
#> ℹ [2026-03-11 16:43:00] Checking a list of <Seurat>...
#> ! [2026-03-11 16:43:01] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-11 16:43:01] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-11 16:43:02] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-11 16:43:03] Use the separate HVF from `srt_list`
#> ℹ [2026-03-11 16:43:03] Number of available HVF: 2000
#> ℹ [2026-03-11 16:43:03] Finished check
#> ℹ [2026-03-11 16:43:03] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-11 16:43:03] Perform pca linear dimension reduction
#> ℹ [2026-03-11 16:43:04] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-11 16:43:04] Reorder clusters...
#> ℹ [2026-03-11 16:43:05] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-11 16:43:05] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-11 16:43:08] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-11 16:43:11] Run scop standard workflow completed
p1 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  label = TRUE
)
p1


thisplot::panel_fix(
  p1, height = 2, width = 3
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  stat_type = "count",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  bg.by = "CellType",
  palette = "Set1",
  stat_type = "count",
  position = "dodge"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "bar"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "rose"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "ring"
)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "pie"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "dot"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "rose"
)
#> Error in loadNamespace(x): there is no package called ‘geomtextpath’

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "ring"
)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "area"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "dot"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "trend"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar",
  individual = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "bar"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose"
)
#> Error in loadNamespace(x): there is no package called ‘geomtextpath’

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring"
)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "area"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "dot"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "trend"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "bar",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose",
  position = "dodge",
  label = TRUE
)
#> Error in loadNamespace(x): there is no package called ‘geomtextpath’

CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring",
  position = "dodge",
  label = TRUE
)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "sankey"
)
#> ! [2026-03-11 16:43:17] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "chord"
)
#> ! [2026-03-11 16:43:18] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "venn",
  stat_level = list(
    CellType = c("Ductal", "Ngn3-low-EP"),
    Phase = "S"
  )
)
#> ! [2026-03-11 16:43:19] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"
#> Error in loadNamespace(x): there is no package called ‘ggVennDiagram’

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
#> ! [2026-03-11 16:43:19] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"
#> Error in loadNamespace(x): there is no package called ‘ggVennDiagram’

CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
  ),
  plot_type = "upset",
  stat_level = "TRUE"
)
#> ! [2026-03-11 16:43:19] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"
#> Error in loadNamespace(x): there is no package called ‘ggupset’

sum(
  pancreas_sub$Progenitor == "FALSE" &
    pancreas_sub$G2M == "FALSE" &
    pancreas_sub$Fancb_Expressed == "TRUE" &
    pancreas_sub$Dlg3_Expressed == "FALSE"
)
#> [1] 6
```
