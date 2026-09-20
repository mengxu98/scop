# Velocity Plot

Creates a velocity plot for a given Seurat object. The plot shows the
velocity vectors of the cells in a specified reduction space.

## Usage

``` r
VelocityPlot(
  object,
  reduction,
  dims = c(1, 2),
  cells = NULL,
  velocity = "stochastic",
  plot_type = c("raw", "grid", "stream"),
  group.by = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  n_neighbors = ceiling(ncol(srt@assays[[1]])/50),
  density = 1,
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
  title = "Cell velocity",
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  return_layer = FALSE,
  seed = 11,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Length-2 vector of dimensions to plot.

- cells:

  Cell names to include.

- velocity:

  Name of the velocity to use for plotting. Default is `"stochastic"`.
  If the corresponding velocity embedding (e.g. `stochastic_umap`) is
  not found, falls back to the scVelo convention (`velocity_umap`),
  which is what an object converted from an AnnData (via
  [adata_to_srt](https://mengxu98.github.io/scop/reference/adata_to_srt.md))
  contains.

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

  Proportion of cells to plot. Defaults is `1` (plot all cells).

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

- streamline_L:

  Length of the streamlines. Defaults is `5`.

- streamline_minL:

  Minimum length of the streamlines. Defaults is `1`.

- streamline_res:

  Resolution of the streamlines. Defaults is `1`.

- streamline_n:

  Number of streamlines to plot. Defaults is `15`.

- streamline_width:

  Width of the streamlines. Defaults is `c(0, 0.8)`.

- streamline_alpha:

  Alpha transparency of the streamlines. Defaults is `1`.

- streamline_color:

  Color of the streamlines. Defaults is `NULL`.

- streamline_palette:

  Name of the palette to use for coloring the streamlines. Defaults is
  `"RdYlBu"`.

- streamline_palcolor:

  Colors to use for coloring the streamlines. Defaults is `NULL`.

- streamline_bg_color:

  Background color of the streamlines. Defaults is `"white"`.

- streamline_bg_stroke:

  Stroke width of the streamlines background. Defaults is `0.5`.

- aspect.ratio:

  Panel aspect ratio.

- title:

  The text for the title. Defaults is `"Cell velocity"`.

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

- seed:

  Random seed.

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## See also

[RunSCVELO](https://mengxu98.github.io/scop/reference/RunSCVELO.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 22:58:32] Start standard processing workflow...
#> ℹ [2026-09-20 22:58:32] Checking a list of <Seurat>...
#> ! [2026-09-20 22:58:32] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 22:58:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:58:32] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:58:33] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:58:33] Number of available HVF: 2000
#> ℹ [2026-09-20 22:58:33] Finished check
#> ℹ [2026-09-20 22:58:33] Perform `ScaleData()`
#> ℹ [2026-09-20 22:58:33] Perform pca linear dimension reduction
#> ℹ [2026-09-20 22:58:33] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 22:58:34] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 22:58:34] Reorder clusters...
#> ℹ [2026-09-20 22:58:34] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:58:34] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 22:58:42] Standard processing workflow completed
pancreas_sub <- RunSCVELO(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "pca",
  nonlinear_reduction = "umap",
  backend = "cpp",
  show_plot = FALSE,
  return_seurat = TRUE
)
#> ℹ [2026-09-20 22:58:42] Running scanpy-compatible preprocessing (15998 features -> filter + normalize)...
#> ℹ [2026-09-20 22:58:45] Running scVelo "stochastic" mode with `backend = 'cpp'` (9699 features)
#> ✔ [2026-09-20 22:58:48] scVelo "stochastic" mode completed
#> ✔ [2026-09-20 22:58:48] scVelo cpp backend completed
VelocityPlot(
  pancreas_sub,
  reduction = "umap"
)


VelocityPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "SubCellType"
)


VelocityPlot(
  pancreas_sub,
  reduction = "umap",
  plot_type = "grid"
)


VelocityPlot(
  pancreas_sub,
  reduction = "umap",
  plot_type = "stream"
)


VelocityPlot(
  pancreas_sub,
  reduction = "umap",
  plot_type = "stream",
  streamline_color = "black"
)


VelocityPlot(
  pancreas_sub,
  reduction = "umap",
  plot_type = "stream",
  streamline_color = "black",
  arrow_color = "red"
)
```
