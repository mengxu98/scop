# 3D-Dimensional reduction plot for cell classification visualization.

Plotting cell points on a reduced 3D space and coloring according to the
groups of the cells.

## Usage

``` r
CellDimPlot3D(
  srt,
  group.by,
  reduction = NULL,
  dims = c(1, 2, 3),
  axis_labs = NULL,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey80",
  pt.size = 1.5,
  cells.highlight = NULL,
  cols.highlight = "black",
  shape.highlight = "circle-open",
  sizes.highlight = 2,
  lineages = NULL,
  lineages_palette = "Dark2",
  span = 0.75,
  width = NULL,
  height = NULL,
  save = NULL,
  force = FALSE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a three-length numeric vector specifying
  x-, y- and z-dimensions

- axis_labs:

  A character vector of length 3 indicating the labels for the axes.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- bg_color:

  Color value for background(NA) points.

- pt.size:

  The size of the points in the plot.

- cells.highlight:

  A logical or character vector specifying the cells to highlight in the
  plot. If `TRUE`, all cells are highlighted. If `FALSE`, no cells are
  highlighted. Default is `NULL`.

- cols.highlight:

  Color used to highlight the cells.

- shape.highlight:

  Shape of the cell to highlight. See
  [scattergl-marker-symbol](https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol)

- sizes.highlight:

  Size of highlighted cell points.

- lineages:

  Lineages/pseudotime to add to the plot. If specified, curves will be
  fitted using [stats::loess](https://rdrr.io/r/stats/loess.html)
  method.

- lineages_palette:

  Color palette used for lineages.

- span:

  The span of the loess smoother for lineages line.

- width:

  Width in pixels, defaults to automatic sizing.

- height:

  Height in pixels, defaults to automatic sizing.

- save:

  The name of the file to save the plot to. Must end in ".html".

- force:

  Whether to force drawing regardless of maximum levels in any cell
  group is greater than 100. Default is `FALSE`.

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[FeatureDimPlot3D](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 15:27:16] Start standard processing workflow...
#> ℹ [2026-04-02 15:27:17] Checking a list of <Seurat>...
#> ! [2026-04-02 15:27:17] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 15:27:17] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:27:19] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:27:19] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 15:27:19] Number of available HVF: 2000
#> ℹ [2026-04-02 15:27:19] Finished check
#> ℹ [2026-04-02 15:27:20] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 15:27:20] Perform pca linear dimension reduction
#> ℹ [2026-04-02 15:27:24] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 15:27:25] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T15:27:25. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-7); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-7 (‘future_lapply-1’)
#> !                       
#> !                       DEBUG: BEGIN TROUBLESHOOTING HELP
#> !                       SequentialFuture:
#> !                       Label: ‘future_lapply-1’
#> !                       Expression:
#> Error in glue(str, .envir = .envir, .transformer = transformer, .cli = TRUE,     .trim = .trim): Expecting '}'
CellDimPlot3D(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D"
)
#> Error in DefaultReduction(srt, pattern = reduction, min_dim = 3): Unable to find any reductions

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D",
  show_plot = FALSE
)
#> Error in DefaultReduction(srt, pattern = reduction): Unable to find any reductions
CellDimPlot3D(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D",
  lineages = "Lineage1"
)
#> Error in CellDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D",     lineages = "Lineage1"): Lineage1 is not in the meta.data of srt object.
```
