# Projection Plot

This function generates a projection plot, which can be used to compare
two groups of cells in a dimensionality reduction space.

## Usage

``` r
ProjectionPlot(
  srt_query,
  srt_ref,
  query_group = NULL,
  ref_group = NULL,
  query_reduction = "ref.embeddings",
  ref_reduction = srt_query[[query_reduction]]@misc[["reduction.model"]] %||% NULL,
  query_param = list(palette = "Set1", cells.highlight = TRUE),
  ref_param = list(palette = "Chinese"),
  xlim = NULL,
  ylim = NULL,
  pt.size = 0.8,
  stroke.highlight = 0.5,
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  A Seurat object or count matrix representing the reference object. If
  provided, the similarities will be calculated between cells from the
  query and reference objects. If not provided, the similarities will be
  calculated within the query object.

- query_group:

  The grouping variable for the query group cells.

- ref_group:

  The grouping variable for the reference group cells.

- query_reduction:

  The name of the reduction in the query group cells.

- ref_reduction:

  The name of the reduction in the reference group cells.

- query_param:

  A list of parameters for customizing the query group plot. Available
  parameters: palette (color palette for groups) and cells.highlight
  (whether to highlight cells).

- ref_param:

  A list of parameters for customizing the reference group plot.
  Available parameters: palette (color palette for groups) and
  cells.highlight (whether to highlight cells).

- xlim:

  The x-axis limits for the plot. If not provided, the limits will be
  calculated based on the data.

- ylim:

  The y-axis limits for the plot. If not provided, the limits will be
  calculated based on the data.

- pt.size:

  The size of the points in the plot.

- stroke.highlight:

  The size of the stroke highlight for cells.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-06-24 18:32:23] Start standard processing workflow...
#> ℹ [2026-06-24 18:32:23] Checking a list of <Seurat>...
#> ! [2026-06-24 18:32:23] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 18:32:23] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 18:32:23] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 18:32:24] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 18:32:24] Number of available HVF: 2000
#> ℹ [2026-06-24 18:32:24] Finished check
#> ℹ [2026-06-24 18:32:24] Perform `ScaleData()`
#> ℹ [2026-06-24 18:32:24] Perform pca linear dimension reduction
#> ℹ [2026-06-24 18:32:24] Use stored estimated dimensions 1:27 for Standardpca
#> ℹ [2026-06-24 18:32:25] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-24 18:32:25] Reorder clusters...
#> ℹ [2026-06-24 18:32:25] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-24 18:32:25] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-24 18:32:31] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-06-24 18:32:32] Run integration workflow...
#> ℹ [2026-06-24 18:32:32] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-06-24 18:32:33] Checking a list of <Seurat>...
#> ℹ [2026-06-24 18:32:33] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-24 18:32:33] Perform `FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-06-24 18:32:33] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-24 18:32:33] Perform `FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-06-24 18:32:33] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-24 18:32:34] Perform `FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-06-24 18:32:34] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-06-24 18:32:34] Perform `FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-06-24 18:32:34] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 18:32:35] Number of available HVF: 2000
#> ℹ [2026-06-24 18:32:35] Finished check
#> ℹ [2026-06-24 18:32:36] Perform Uncorrected integration
#> ℹ [2026-06-24 18:32:37] Perform `Seurat::ScaleData()`
#> Error: ScaleData.Seurat requires an Assay5 object with a data layer.
CellDimPlot(
  srt_ref,
  group.by = c("celltype", "tech")
)


# Projection
srt_query <- RunKNNMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_umap = "UncorrectedUMAP2D"
)
#> Error in srt_ref[[ref_umap]]: ‘UncorrectedUMAP2D’ not found in this Seurat object
#>  
ProjectionPlot(
  srt_query = srt_query,
  srt_ref = srt_ref,
  query_group = "celltype",
  ref_group = "celltype"
)
#> Error in srt_query[[query_reduction]]: ‘ref.embeddings’ not found in this Seurat object
#>  
```
