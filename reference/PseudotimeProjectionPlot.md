# Pseudotime Projection Plot

Creates a projection plot similar to
[VelocityPlot](https://mengxu98.github.io/scop/reference/VelocityPlot.md),
but uses pseudotime data instead of RNA velocity analysis results.

## Usage

``` r
PseudotimeProjectionPlot(
  srt,
  reduction,
  time_key,
  dims = c(1, 2),
  cells = NULL,
  method = c("knn", "gradient"),
  k = 30,
  graph_name = NULL,
  plot_type = c("raw", "grid", "stream"),
  group.by = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  n_neighbors = ceiling(ncol(srt@assays[[1]])/50),
  density = 2,
  smooth = 0.5,
  scale = 1,
  min_mass = 1,
  cutoff_perc = 5,
  arrow_angle = 20,
  arrow_color = "black",
  streamline_L = 5,
  streamline_minL = 1,
  streamline_res = 1,
  streamline_n = 15,
  streamline_width = c(0, 0.8),
  streamline_alpha = 1,
  streamline_color = NULL,
  streamline_palette = "RdYlBu",
  streamline_palcolor = NULL,
  streamline_bg_color = "white",
  streamline_bg_stroke = 0.5,
  aspect.ratio = 1,
  title = "Pseudotime projection",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  palette = NULL,
  palcolor = NULL,
  show_cells = TRUE,
  pt.size = 2,
  pt.alpha = 0.3,
  label = NULL,
  label.size = 4,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- time_key:

  Name of the column in the Seurat object metadata containing pseudotime
  values.

- dims:

  Length-2 vector of dimensions to plot.

- cells:

  Cell names to include.

- method:

  Method to compute velocity vectors from pseudotime. Can be
  `"gradient"` or `"knn"`.

- k:

  Number of nearest neighbors to use when `method = "knn"`.

- graph_name:

  Name of the KNN graph in the Seurat object to use. If `NULL`, a new
  graph will be computed. Default is `NULL`.

- plot_type:

  Type of plot to create. Can be `"raw"`, `"grid"`, or `"stream"`.

- group.by:

  Metadata column(s) used to color cells.

- group_palette:

  Name of the palette to use for coloring the groups. Defaults is
  `"Chinese"`.

- group_palcolor:

  Colors to use for coloring the groups. Defaults is `NULL`.

- n_neighbors:

  Number of neighbors to include for the density estimation. Defaults is
  `ceiling(ncol(srt@assays[[1]]) / 50)`.

- density:

  Scale for the streamline grid (number of points per axis is
  `ceiling(50 * density)`). Default is `2` (CellRank/scvelo-style).

- smooth:

  Smoothing parameter for density estimation. Defaults is `0.5`.

- scale:

  Scaling factor for the velocity vectors. Defaults is `1`.

- min_mass:

  Minimum mass value for the density-based cutoff. Defaults is `1`.

- cutoff_perc:

  Percentile value for the density-based cutoff. Defaults is `5`.

- arrow_angle:

  Angle of the arrowheads. Defaults is `20`.

- arrow_color:

  Color of the arrowheads. Defaults is `"black"`.

- streamline_L, streamline_minL, streamline_res, streamline_n,
  streamline_width, streamline_alpha, streamline_color,
  streamline_palette, streamline_palcolor, streamline_bg_color,
  streamline_bg_stroke:

  Streamline appearance for velocity plots.

- aspect.ratio:

  Panel aspect ratio.

- title:

  The text for the title. Defaults is `"Pseudotime projection"`.

- subtitle:

  Plot subtitle.

- xlab:

  Label for the x axis.

- ylab:

  Label for the y axis.

- legend.position:

  Legend position passed to `theme()`.

- legend.direction:

  Legend direction passed to `theme()`.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- return_layer:

  Whether to return the plot layers as a list. Defaults is `FALSE`.

- palette:

  Deprecated alias of `group_palette`.

- palcolor:

  Deprecated alias of `group_palcolor`.

- show_cells:

  Whether to show cell points on the plot. Defaults is `TRUE`.

- pt.size:

  Size of cell points. Defaults is `2` (CellRank-style overlapping
  patches).

- pt.alpha:

  The transparency of the data points. Default is `0.3`
  (CellRank/scvelo-style).

- label:

  Whether to label the cell groups. Defaults is `TRUE` when `group.by`
  is specified.

- label.size:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- label.fg:

  Foreground color of labels. Defaults is `"black"`.

- label.bg:

  Background color of labels. Defaults is `"white"`.

- label.bg.r:

  Group labels. `label_insitu = FALSE` uses numbers instead of group
  names.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[VelocityPlot](https://mengxu98.github.io/scop/reference/VelocityPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:50:39] Start standard processing workflow...
#> ℹ [2026-09-06 21:50:39] Checking a list of <Seurat>...
#> ! [2026-09-06 21:50:40] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:50:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:50:40] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:50:40] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:50:40] Number of available HVF: 2000
#> ℹ [2026-09-06 21:50:40] Finished check
#> ℹ [2026-09-06 21:50:40] Perform `ScaleData()`
#> ℹ [2026-09-06 21:50:40] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:50:40] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:50:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:50:41] Reorder clusters...
#> ℹ [2026-09-06 21:50:41] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:50:41] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:50:47] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  reduction = "UMAP",
  group.by = "SubCellType"
)


PseudotimeProjectionPlot(
  pancreas_sub,
  reduction = "UMAP",
  group.by = "SubCellType",
  time_key = "Lineage1",
  method = "gradient",
  plot_type = "raw"
)
#> ! [2026-09-06 21:50:49] Removed 312 cells with NA pseudotime values


PseudotimeProjectionPlot(
  pancreas_sub,
  reduction = "UMAP",
  time_key = "Lineage1",
  group.by = "SubCellType",
  plot_type = "stream",
  show_cells = TRUE,
  label = TRUE
)
#> ! [2026-09-06 21:50:49] Removed 312 cells with NA pseudotime values
#> ℹ [2026-09-06 21:50:49] Computing KNN graph from embedding...


PseudotimeProjectionPlot(
  pancreas_sub,
  reduction = "UMAP",
  time_key = "Lineage2",
  plot_type = "grid"
)
#> ! [2026-09-06 21:50:51] Removed 354 cells with NA pseudotime values
#> ℹ [2026-09-06 21:50:51] Computing KNN graph from embedding...


PseudotimeProjectionPlot(
  pancreas_sub,
  reduction = "UMAP",
  time_key = "Lineage1",
  method = "gradient",
  plot_type = "raw"
)
#> ! [2026-09-06 21:50:53] Removed 312 cells with NA pseudotime values
```
