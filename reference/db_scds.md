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
#> ℹ [2026-01-27 08:36:35] Start standard scop workflow...
#> ℹ [2026-01-27 08:36:35] Checking a list of <Seurat>...
#> ! [2026-01-27 08:36:35] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-27 08:36:35] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 08:36:38] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 08:36:38] Use the separate HVF from srt_list
#> ℹ [2026-01-27 08:36:38] Number of available HVF: 2000
#> ℹ [2026-01-27 08:36:38] Finished check
#> ℹ [2026-01-27 08:36:39] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-27 08:36:39] Perform pca linear dimension reduction
#> ℹ [2026-01-27 08:36:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-27 08:36:40] Reorder clusters...
#> ℹ [2026-01-27 08:36:40] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-27 08:36:40] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-27 08:36:45] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-27 08:36:49] Run scop standard workflow completed
pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#> ℹ [2026-01-27 08:36:49] Data type is raw counts
#> Registered S3 method overwritten by 'pROC':
#>   method   from            
#>   plot.roc spatstat.explore
#> Error in xgboost(mm, nrounds = nmax, tree_method = "hist", nthread = 2,     early_stopping_rounds = 2, subsample = 0.5, objective = "binary:logistic",     verbose = 0): argument "y" is missing, with no default
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.scds_hybrid_class"
)
#> Error in CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.scds_hybrid_class"): "db.scds_hybrid_class" is not in the meta.data of srt object

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.scds_hybrid_score"
)
#> ! [2026-01-27 08:37:11] "db.scds_hybrid_score" are not in the features of <Seurat>
#> Error in FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.scds_hybrid_score"): There are no valid features present.
```
