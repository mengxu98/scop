# Run doublet-calling with scds

Run doublet-calling with scds

## Usage

``` r
db_scds(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  method = c("hybrid", "cxds", "bcds"),
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for doublet-calling. Default is
  `"RNA"`.

- db_rate:

  The expected doublet rate. Default is calculated as
  `ncol(srt) / 1000 * 0.01`.

- method:

  The method to be used for doublet-calling. Options are `"hybrid"`,
  `"cxds"`, or `"bcds"`.

- ...:

  Additional arguments to be passed to
  [`scds::cxds_bcds_hybrid()`](https://rdrr.io/pkg/scds/man/cxds_bcds_hybrid.html).

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-03 04:45:51] Start standard processing workflow...
#> ℹ [2026-04-03 04:45:52] Checking a list of <Seurat>...
#> ! [2026-04-03 04:45:52] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:45:52] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 04:45:55] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 04:45:55] Use the separate HVF from `srt_list`
#> ℹ [2026-04-03 04:45:55] Number of available HVF: 2000
#> ℹ [2026-04-03 04:45:55] Finished check
#> ℹ [2026-04-03 04:45:56] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-03 04:45:56] Perform pca linear dimension reduction
#> ℹ [2026-04-03 04:45:56] Use stored estimated dimensions 1:12 for Standardpca
#> ℹ [2026-04-03 04:45:57] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-03 04:45:57] Reorder clusters...
#> ℹ [2026-04-03 04:45:57] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-03 04:45:57] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-03 04:45:57] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ℹ [2026-04-03 04:46:02] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ✔ [2026-04-03 04:46:07] Standard processing workflow completed
pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#> ℹ [2026-04-03 04:46:07] Data type is raw counts
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
