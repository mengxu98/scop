# Cell density plot

Feature density by group.

## Usage

``` r
CellDensityPlot(
  srt,
  features,
  group.by = NULL,
  split.by = NULL,
  assay = NULL,
  layer = "data",
  flip = FALSE,
  reverse = FALSE,
  x_order = c("value", "rank"),
  decreasing = NULL,
  palette = "Chinese",
  palcolor = NULL,
  cells = NULL,
  keep_empty = FALSE,
  y.nbreaks = 4,
  y.min = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- features:

  Features to plot.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column to facet by.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- flip, reverse:

  Flip the x-axis or reverse the y-axis.

- x_order:

  `"value"` or `"rank"`.

- decreasing:

  Order groups decreasingly.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- cells:

  Cell names to include. `NULL` uses all cells.

- keep_empty:

  Keep empty groups.

- y.nbreaks, y.min, y.max:

  Y-axis breaks and limits. `NULL` limits are automatic.

- same.y.lims:

  Share y-axis limits across panels.

- aspect.ratio:

  Panel aspect ratio.

- title:

  Plot title. `NULL` hides the title for merged/single panels. When
  multiple lineages are plotted and `title` is `NULL`, each panel is
  titled with its lineage column.

- subtitle:

  Plot subtitle.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- force:

  Draw even when there are more than 50 features.

- verbose:

  Whether to print messages.

## See also

[CellStatPlot](https://mengxu98.github.io/scop/reference/CellStatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:01:05] Start standard processing workflow...
#> ℹ [2026-09-06 21:01:06] Checking a list of <Seurat>...
#> ! [2026-09-06 21:01:06] Data 1/1 of the `srt_list` is "unknown"
#> Warning: Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:01:06] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:01:06] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:01:06] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:01:06] Number of available HVF: 2000
#> ℹ [2026-09-06 21:01:06] Finished check
#> ℹ [2026-09-06 21:01:06] Perform `ScaleData()`
#> ℹ [2026-09-06 21:01:07] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:01:07] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:01:07] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:01:07] Reorder clusters...
#> ℹ [2026-09-06 21:01:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:01:07] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:01:11] Standard processing workflow completed
CellDensityPlot(
  pancreas_sub,
  features = "Sox9",
  group.by = "SubCellType"
)
#> Picking joint bandwidth of 0.209


pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_path()`).


CellDensityPlot(
  pancreas_sub,
  features = "Lineage1",
  group.by = "SubCellType",
  aspect.ratio = 1
)
#> Picking joint bandwidth of 0.573


CellDensityPlot(
  pancreas_sub,
  features = "Lineage1",
  group.by = "SubCellType",
  flip = TRUE
)
#> Picking joint bandwidth of 0.573
```
