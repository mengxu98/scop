# Projection Plot

Generates a projection plot, which can be used to compare two groups of
cells in a dimensionality reduction space.

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

  Reduction in the query group cells.

- ref_reduction:

  Reduction in the reference group cells.

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

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- stroke.highlight:

  The size of the stroke highlight for cells.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- RunStandardWorkflow(panc8_sub)
#> ℹ [2026-10-04 22:09:04] Start standard processing workflow...
#> ℹ [2026-10-04 22:09:04] Checking a list of <Seurat>...
#> ! [2026-10-04 22:09:04] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-04 22:09:04] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:09:04] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:09:04] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:09:04] Number of available HVF: 2000
#> ℹ [2026-10-04 22:09:04] Finished check
#> ℹ [2026-10-04 22:09:04] Perform `ScaleData()`
#> ℹ [2026-10-04 22:09:04] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:09:05] Use stored estimated dimensions 1:26 for Standardpca
#> ℹ [2026-10-04 22:09:05] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:09:05] Reorder clusters...
#> ℹ [2026-10-04 22:09:05] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:09:05] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-04 22:09:13] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- RunIntegration(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-10-04 22:09:13] Run integration workflow...
#> ℹ [2026-10-04 22:09:14] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-04 22:09:14] Checking a list of <Seurat>...
#> ℹ [2026-10-04 22:09:14] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:09:14] Perform `FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-10-04 22:09:15] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:09:15] Perform `FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-10-04 22:09:15] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:09:15] Perform `FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-10-04 22:09:15] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:09:15] Perform `FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-10-04 22:09:15] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:09:15] Number of available HVF: 2000
#> ℹ [2026-10-04 22:09:15] Finished check
#> ℹ [2026-10-04 22:09:15] Perform Uncorrected integration
#> ℹ [2026-10-04 22:09:15] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-04 22:09:15] Perform "pca" linear dimension reduction
#> ℹ [2026-10-04 22:09:16] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-04 22:09:16] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-04 22:09:16] Reorder clusters...
#> ℹ [2026-10-04 22:09:16] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:09:16] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-10-04 22:09:21] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-10-04 22:09:26] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ✔ [2026-10-04 22:09:32] Uncorrected integration completed
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
#> ℹ [2026-10-04 22:09:32] Use the features to calculate distance metric
#> ℹ [2026-10-04 22:09:32] Data type is log-normalized
#> ℹ [2026-10-04 22:09:33] Data type is log-normalized
#> ℹ [2026-10-04 22:09:33] Use 675 features to calculate distance
#> ℹ [2026-10-04 22:09:33] Use raw method to find neighbors
#> ℹ [2026-10-04 22:09:33] Running UMAP projection
ProjectionPlot(
  srt_query = srt_query,
  srt_ref = srt_ref,
  query_group = "celltype",
  ref_group = "celltype"
)
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
```
