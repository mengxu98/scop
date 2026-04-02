# 3D-Dimensional reduction plot for gene expression visualization.

3D-Dimensional reduction plot for gene expression visualization.

## Usage

``` r
FeatureDimPlot3D(
  srt,
  features,
  reduction = NULL,
  dims = c(1, 2, 3),
  axis_labs = NULL,
  split.by = NULL,
  layer = "data",
  assay = NULL,
  calculate_coexp = FALSE,
  pt.size = 1.5,
  cells.highlight = NULL,
  cols.highlight = "black",
  shape.highlight = "circle-open",
  sizes.highlight = 2,
  width = NULL,
  height = NULL,
  save = NULL,
  force = FALSE
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A character vector or a named list of features to plot. Features can
  be gene names in Assay or names of numeric columns in meta.data.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a two-length numeric vector specifying x-
  and y-dimensions

- axis_labs:

  A character vector of length 3 indicating the labels for the axes.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- calculate_coexp:

  Whether to calculate the co-expression value (geometric mean) of the
  features.

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

- width:

  Width in pixels, defaults to automatic sizing.

- height:

  Height in pixels, defaults to automatic sizing.

- save:

  The name of the file to save the plot to. Must end in ".html".

- force:

  Whether to force drawing regardless of the number of features greater
  than 100. Default is `FALSE`.

## See also

[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[CellDimPlot3D](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 16:13:28] Start standard processing workflow...
#> ℹ [2026-04-02 16:13:29] Checking a list of <Seurat>...
#> ! [2026-04-02 16:13:29] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:13:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:13:30] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:13:31] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:13:31] Number of available HVF: 2000
#> ℹ [2026-04-02 16:13:31] Finished check
#> ℹ [2026-04-02 16:13:32] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:13:32] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:13:36] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 16:13:36] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T16:13:36. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-19); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-19 (‘future_lapply-1’)
#> !                       
#> !                       DEBUG: BEGIN TROUBLESHOOTING HELP
#> !                       SequentialFuture:
#> !                       Label: ‘future_lapply-1’
#> !                       Expression:
#> Error in glue(str, .envir = .envir, .transformer = transformer, .cli = TRUE,     .trim = .trim): Expecting '}'
FeatureDimPlot3D(
  pancreas_sub,
  features = c("Ghrl", "Ins1", "Gcg", "Ins2"),
  reduction = "StandardpcaUMAP3D"
)
#> Error in DefaultReduction(srt, pattern = reduction, min_dim = 3): Unable to find any reductions
```
