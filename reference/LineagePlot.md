# Lineage Plot

Generate a lineage plot based on the pseudotime.

## Usage

``` r
LineagePlot(
  srt,
  lineages,
  reduction = NULL,
  dims = c(1, 2),
  cells = NULL,
  trim = c(0.01, 0.99),
  span = 0.75,
  palette = "Dark2",
  palcolor = NULL,
  lineages_arrow = grid::arrow(length = grid::unit(0.1, "inches")),
  linewidth = 1,
  line_bg = "white",
  line_bg_stroke = 0.5,
  whiskers = FALSE,
  whiskers_linewidth = 0.5,
  whiskers_alpha = 0.5,
  aspect.ratio = 1,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- lineages:

  A character vector that specifies the lineages to be included.
  Typically, use the pseudotime of cells.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Length-2 vector of dimensions to plot.

- cells:

  Cell names to include.

- trim:

  A numeric vector of length 2 specifying the quantile range of lineages
  to include in the plot.

- span:

  The span of the loess smoother.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).

- palcolor:

  Custom colors used to create a color palette.

- lineages_arrow:

  An arrow object specifying the arrow for lineages.

- linewidth:

  The linewidth for the lineages.

- line_bg:

  Color for the background lines.

- line_bg_stroke:

  The stroke width for the background lines.

- whiskers:

  Whether to include whiskers in the plot.

- whiskers_linewidth:

  The linewidth for the whiskers.

- whiskers_alpha:

  The transparency for the whiskers.

- aspect.ratio:

  Panel aspect ratio.

- title, subtitle, xlab, ylab:

  Plot labels.

- legend.position:

  Legend position passed to `theme()`.

- legend.direction:

  Legend direction passed to `theme()`.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- return_layer:

  Whether to return the plot layers as a list. Defaults is `FALSE`.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[RunSlingshot](https://mengxu98.github.io/scop/reference/RunSlingshot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:37:35] Start standard processing workflow...
#> ℹ [2026-09-06 21:37:35] Checking a list of <Seurat>...
#> ! [2026-09-06 21:37:36] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:37:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:37:36] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:37:36] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:37:36] Number of available HVF: 2000
#> ℹ [2026-09-06 21:37:36] Finished check
#> ℹ [2026-09-06 21:37:36] Perform `ScaleData()`
#> ℹ [2026-09-06 21:37:36] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:37:36] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:37:37] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:37:37] Reorder clusters...
#> ℹ [2026-09-06 21:37:37] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:37:37] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:37:43] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  show_plot = FALSE
)
LineagePlot(
  pancreas_sub,
  lineages = paste0("Lineage", 1:2)
)

LineagePlot(
  pancreas_sub,
  lineages = paste0("Lineage", 1:2),
  whiskers = TRUE
)
```
