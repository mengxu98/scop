# Run doublet-calling with scds

Run doublet-calling with scds

## Usage

``` r
db_scds(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  method = c("hybrid", "cxds", "bcds"),
  data_type = NULL,
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

- data_type:

  Optional precomputed result from
  [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  for the input assay. Primarily used internally to avoid repeated scans
  of the same count matrix across nested QC calls.

- ...:

  Additional arguments to be passed to
  [`scds::cxds_bcds_hybrid()`](https://rdrr.io/pkg/scds/man/cxds_bcds_hybrid.html).

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-11 16:34:28] Start standard processing workflow...
#> ℹ [2026-05-11 16:34:28] Checking a list of <Seurat>...
#> ! [2026-05-11 16:34:28] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:34:28] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 16:34:30] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 16:34:31] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 16:34:31] Number of available HVF: 2000
#> ℹ [2026-05-11 16:34:31] Finished check
#> ℹ [2026-05-11 16:34:31] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 16:34:31] Perform pca linear dimension reduction
#> ℹ [2026-05-11 16:34:32] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-11 16:34:32] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-11 16:34:32] Reorder clusters...
#> ℹ [2026-05-11 16:34:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-11 16:34:32] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-11 16:34:32] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-11 16:34:37] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-11 16:34:42] Standard processing workflow completed
pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#> ℹ [2026-05-11 16:34:42] Running scds with method "hybrid"
#> ℹ [2026-05-11 16:34:43] Data type is raw counts
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
