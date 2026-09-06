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
  plot_type = c("bar", "rose", "ring", "pie", "trend", "trend_alluvial", "area", "dot",
    "sankey", "chord", "venn", "upset"),
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
  x_text_angle = 45,
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- stat.by:

  Features to plot.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column to facet by.

- bg.by:

  Metadata column used as background color.

- cells:

  Cell names to include.

- flip:

  Flip x and y.

- NA_color:

  The color to use for missing values.

- NA_stat:

  Whether to include missing values in the plot.

- keep_empty:

  Keep empty factor levels.

- individual:

  One plot per group.

- stat_level:

  The level(s) of the variable(s) specified in `stat.by` to include in
  the plot.

- plot_type:

  The type of plot to create. Can be one of `"bar"`, `"rose"`, `"ring"`,
  `"pie"`, `"trend"`, `"trend_alluvial"`, `"area"`, `"dot"`, `"sankey"`,
  `"chord"`, `"venn"`, or `"upset"`.

- stat_type:

  The type of statistic to compute for the plot. Can be one of
  `"percent"` or `"count"`.

- position:

  The position adjustment for the plot. Can be one of `"stack"` or
  `"dodge"`.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- alpha:

  Plot transparency.

- bg_palette, bg_palcolor, bg_alpha:

  Background palette and transparency.

- label:

  Whether to add labels on the plot.

- label.size:

  The size of the labels.

- label.fg:

  The foreground color of the labels.

- label.bg:

  The background color of the labels.

- label.bg.r:

  The radius of the rounded corners of the label background.

- aspect.ratio:

  Panel aspect ratio.

- title:

  Plot title. `NULL` hides the title for merged/single panels. When
  multiple lineages are plotted and `title` is `NULL`, each panel is
  titled with its lineage column.

- subtitle:

  Plot subtitle.

- xlab:

  Plot labels.

- ylab:

  Y-axis label.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- x_text_angle:

  Rotation angle for x-axis labels.

- grid_major:

  Whether to show major panel grid lines.

- grid_major_colour:

  Color of major panel grid lines.

- grid_major_linetype:

  Linetype of major panel grid lines.

- grid_major_linewidth:

  Line width of major panel grid lines.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- force:

  Draw even when a grouping has more than 100 levels.

- seed:

  Random seed.

- ...:

  Additional arguments passed to the plotting helpers.

## See also

[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:04:02] Start standard processing workflow...
#> ℹ [2026-09-06 21:04:03] Checking a list of <Seurat>...
#> ! [2026-09-06 21:04:03] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:04:03] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:04:03] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:04:03] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:04:03] Number of available HVF: 2000
#> ℹ [2026-09-06 21:04:03] Finished check
#> ℹ [2026-09-06 21:04:03] Perform `ScaleData()`
#> ℹ [2026-09-06 21:04:03] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:04:04] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:04:04] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:04:04] Reorder clusters...
#> ℹ [2026-09-06 21:04:04] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:04:04] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:04:10] Standard processing workflow completed
p1 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  label = TRUE
)
p1


thisplot::panel_fix(
  p1,
  height = 2, width = 3
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


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "ring"
)


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
  plot_type = "trend_alluvial"
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


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring"
)


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


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "sankey"
)
#> ! [2026-09-06 21:04:24] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "chord"
)
#> ! [2026-09-06 21:04:25] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "venn",
  stat_level = list(
    CellType = c("Ductal", "Ngn3-low-EP"),
    Phase = "S"
  )
)
#> ! [2026-09-06 21:04:25] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


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
#> ! [2026-09-06 21:04:25] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
  ),
  plot_type = "upset",
  stat_level = "TRUE"
)
#> ! [2026-09-06 21:04:26] `stat_type` is forcibly set to "count" when plot "sankey", "chord", "venn", and "upset"


sum(
  pancreas_sub$Progenitor == "FALSE" &
    pancreas_sub$G2M == "FALSE" &
    pancreas_sub$Fancb_Expressed == "TRUE" &
    pancreas_sub$Dlg3_Expressed == "FALSE"
)
#> [1] 6
```
