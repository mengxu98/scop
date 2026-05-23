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
  ...
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

  The dimensions to use for the Slingshot algorithm. Default is `NULL`,
  which uses first two dimensions.

- start:

  The starting group for the Slingshot algorithm. Default is `NULL`.

- end:

  The ending group for the Slingshot algorithm. Default is `NULL`.

- prefix:

  The prefix to add to the column names of the resulting pseudotime
  variable. Default is `NULL`.

- reverse:

  Logical value indicating whether to reverse the pseudotime variable.
  Default is `FALSE`.

- align_start:

  Logical value indicating whether to align the starting pseudotime
  values at the maximum pseudotime. Default is `FALSE`.

- show_plot:

  Logical value indicating whether to show the dimensionality plot.
  Default is `TRUE`.

- lineage_palette:

  The color palette to use for the lineages in the plot. Default is
  `"Dark2"`.

- seed:

  Random seed for reproducibility. Default is `11`.

- ...:

  Additional arguments to be passed to the
  [slingshot::slingshot](https://rdrr.io/pkg/slingshot/man/slingshot.html)
  function.

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-23 09:18:13] Start standard processing workflow...
#> ℹ [2026-05-23 09:18:14] Checking a list of <Seurat>...
#> ! [2026-05-23 09:18:14] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-23 09:18:14] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-23 09:18:16] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-23 09:18:17] Use the separate HVF from `srt_list`
#> ℹ [2026-05-23 09:18:17] Number of available HVF: 2000
#> ℹ [2026-05-23 09:18:17] Finished check
#> ℹ [2026-05-23 09:18:17] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-23 09:18:18] Perform pca linear dimension reduction
#> ℹ [2026-05-23 09:18:18] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-23 09:18:18] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-23 09:18:18] Reorder clusters...
#> ℹ [2026-05-23 09:18:19] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-23 09:18:19] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-23 09:18:19] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-23 09:18:25] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-23 09:18:31] Standard processing workflow completed
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
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_path()`).

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
