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
  verbose = TRUE,
  srt = NULL,
  ...
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

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- ...:

  Arguments passed to other functions.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-11-11 05:44:20] Start standard processing workflow...
#> ℹ [2026-11-11 05:44:20] Checking a list of <Seurat>...
#> ! [2026-11-11 05:44:20] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-11-11 05:44:20] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-11-11 05:44:20] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-11-11 05:44:20] Use the separate HVF from `srt_list`
#> ℹ [2026-11-11 05:44:20] Number of available HVF: 2000
#> ℹ [2026-11-11 05:44:20] Finished check
#> ℹ [2026-11-11 05:44:20] Perform `ScaleData()`
#> ℹ [2026-11-11 05:44:20] Perform pca linear dimension reduction
#> ℹ [2026-11-11 05:44:21] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-11-11 05:44:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-11-11 05:44:21] Reorder clusters...
#> ℹ [2026-11-11 05:44:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-11 05:44:21] Perform umap nonlinear dimension reduction
#> ✔ [2026-11-11 05:44:30] Standard processing workflow completed
pancreas_sub <- Runscds(pancreas_sub, method = "hybrid")
#> ℹ [2026-11-11 05:44:30] Running scds with method "hybrid"
#> ℹ [2026-11-11 05:44:30] Data type is raw counts
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
