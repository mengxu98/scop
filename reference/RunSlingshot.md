# RunSlingshot

RunSlingshot

## Usage

``` r
RunSlingshot(
  srt,
  group.by,
  reduction = NULL,
  dims = NULL,
  start = NULL,
  end = NULL,
  prefix = NULL,
  reverse = FALSE,
  align_start = FALSE,
  show_plot = TRUE,
  lineage_palette = "Dark2",
  seed = 11,
  ...,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- group.by:

  Metadata column(s) used to color cells.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  The dimensions to use for the Slingshot algorithm. Default is `NULL`,
  which uses first two dimensions.

- start:

  The starting group for the Slingshot algorithm.

- end:

  The ending group for the Slingshot algorithm.

- prefix:

  The prefix to add to the column names of the resulting pseudotime
  variable.

- reverse:

  Whether to reverse the pseudotime variable.

- align_start:

  Whether to align the starting pseudotime values at the maximum
  pseudotime.

- show_plot:

  Whether to show the dimensionality plot.

- lineage_palette:

  The color palette to use for the lineages in the plot.

- seed:

  Random seed.

- ...:

  Passed to the
  [slingshot::slingshot](https://rdrr.io/pkg/slingshot/man/slingshot.html)
  function.

- verbose:

  Whether to print messages.

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:34:18] Start standard processing workflow...
#> ℹ [2026-09-06 22:34:19] Checking a list of <Seurat>...
#> ! [2026-09-06 22:34:19] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:34:19] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:34:19] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:34:19] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:34:19] Number of available HVF: 2000
#> ℹ [2026-09-06 22:34:19] Finished check
#> ℹ [2026-09-06 22:34:19] Perform `ScaleData()`
#> ℹ [2026-09-06 22:34:20] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:34:20] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:34:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:34:20] Reorder clusters...
#> ℹ [2026-09-06 22:34:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:34:20] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:34:28] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "PCA"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1
)


# 3D lineage
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1,
  lineages_trim = c(0.05, 0.95)
)
```
