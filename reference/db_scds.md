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
#> ℹ [2026-04-06 05:20:02] Start standard processing workflow...
#> ℹ [2026-04-06 05:20:03] Checking a list of <Seurat>...
#> ! [2026-04-06 05:20:03] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 05:20:03] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 05:20:05] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 05:20:06] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 05:20:06] Number of available HVF: 2000
#> ℹ [2026-04-06 05:20:06] Finished check
#> ℹ [2026-04-06 05:20:07] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-06 05:20:07] Perform pca linear dimension reduction
#> ℹ [2026-04-06 05:20:08] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-06 05:20:08] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-06 05:20:08] Reorder clusters...
#> ℹ [2026-04-06 05:20:08] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 05:20:08] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-06 05:20:08] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-06 05:20:13] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-06 05:20:17] Standard processing workflow completed
pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#> ℹ [2026-04-06 05:20:18] Data type is raw counts
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
