# Run VECTOR developmental direction inference

Run VECTOR developmental direction inference

## Usage

``` r
RunVECTOR(
  object,
  reduction = NULL,
  pca.reduction = "pca",
  dims = 1:2,
  pca.dims = 1:30,
  grid.n = 30,
  arrow.p = 0.9,
  arrow.ol = 1.5,
  score.name = "VECTOR_Score",
  tool_name = "VECTOR",
  verbose = TRUE,
  backend = c("cpp", "r")
)
```

## Arguments

- object:

  A `Seurat` object or expression matrix with genes in rows and cells in
  columns.

- reduction:

  Embedding reduction used for the direction field.

- pca.reduction:

  PCA-like reduction used to score local polarization.

- dims:

  Dimensions from `reduction`.

- pca.dims:

  Dimensions from `pca.reduction`.

- grid.n:

  Number of bins per embedding axis.

- arrow.p:

  Distance decay parameter used by the original VECTOR arrow weighting.
  Larger values keep more distant grid centers influential.

- arrow.ol:

  Arrow vector length multiplier relative to grid spacing.

- score.name:

  Metadata column for the VECTOR score.

- tool_name:

  Name used in `srt@tools`.

- verbose:

  Whether to print progress messages.

- backend:

  Numerical backend for grid-arrow aggregation. `"cpp"` avoids
  materializing the full grid-to-grid distance matrix; `"r"` is retained
  as the reference implementation.

## Value

A modified `Seurat` object.

## References

Zhang F, Li X, Tian W. Unsupervised inference of developmental
directions for single cells using VECTOR. Cell Reports, 2020.
doi:10.1016/j.celrep.2020.108069.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:38:30] Start standard processing workflow...
#> ℹ [2026-09-06 22:38:30] Checking a list of <Seurat>...
#> ! [2026-09-06 22:38:30] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:38:30] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:38:30] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:38:31] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:38:31] Number of available HVF: 2000
#> ℹ [2026-09-06 22:38:31] Finished check
#> ℹ [2026-09-06 22:38:31] Perform `ScaleData()`
#> ℹ [2026-09-06 22:38:31] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:38:31] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:38:32] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:38:32] Reorder clusters...
#> ℹ [2026-09-06 22:38:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:38:32] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:38:39] Standard processing workflow completed
pancreas_sub <- RunVECTOR(
  pancreas_sub,
  reduction = "umap",
  pca.reduction = "pca",
  verbose = FALSE
)
FeatureDimPlot(pancreas_sub, features = "VECTOR_Score")


VECTORPlot(pancreas_sub, plot_type = "grid")


VECTORPlot(
  pancreas_sub,
  plot_type = "raw",
  group.by = "SubCellType",
  background = "none"
)
```
