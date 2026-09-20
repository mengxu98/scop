# Run doublet-calling with scds

Run doublet-calling with scds

## Usage

``` r
Runscds(
  object,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  method = c("hybrid", "cxds", "bcds"),
  data_type = NULL,
  ...,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay used for doublet calling.

- db_rate:

  Expected doublet rate.

- method:

  The method to be used for doublet-calling. Options are `"hybrid"`,
  `"cxds"`, or `"bcds"`.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- ...:

  Additional arguments passed to the selected `scds` method.

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 22:56:55] Start standard processing workflow...
#> ℹ [2026-09-20 22:56:55] Checking a list of <Seurat>...
#> ! [2026-09-20 22:56:55] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 22:56:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:56:55] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:56:55] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:56:55] Number of available HVF: 2000
#> ℹ [2026-09-20 22:56:55] Finished check
#> ℹ [2026-09-20 22:56:55] Perform `ScaleData()`
#> ℹ [2026-09-20 22:56:55] Perform pca linear dimension reduction
#> ℹ [2026-09-20 22:56:56] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 22:56:56] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 22:56:56] Reorder clusters...
#> ℹ [2026-09-20 22:56:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:56:56] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 22:57:04] Standard processing workflow completed
pancreas_sub <- Runscds(pancreas_sub, method = "hybrid")
#> ℹ [2026-09-20 22:57:04] Running scds with method "hybrid"
#> ℹ [2026-09-20 22:57:05] Data type is raw counts
#> Registered S3 method overwritten by 'pROC':
#>   method   from            
#>   plot.roc spatstat.explore
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.scds_hybrid_class"
)


FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.scds_hybrid_score"
)
```
