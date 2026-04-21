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
#> ℹ [2026-04-21 07:58:27] Start standard processing workflow...
#> ℹ [2026-04-21 07:58:27] Checking a list of <Seurat>...
#> ! [2026-04-21 07:58:27] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-21 07:58:27] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:58:30] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:58:30] Use the separate HVF from `srt_list`
#> ℹ [2026-04-21 07:58:30] Number of available HVF: 2000
#> ℹ [2026-04-21 07:58:31] Finished check
#> ℹ [2026-04-21 07:58:31] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-21 07:58:31] Perform pca linear dimension reduction
#> ℹ [2026-04-21 07:58:32] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-21 07:58:32] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-21 07:58:32] Reorder clusters...
#> ℹ [2026-04-21 07:58:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-21 07:58:32] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-21 07:58:32] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-21 07:58:37] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-21 07:58:42] Standard processing workflow completed
pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#> ℹ [2026-04-21 07:58:42] Data type is raw counts
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
