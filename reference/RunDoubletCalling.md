# Run doublet-calling for single cell RNA-seq data.

Run doublet-calling for single cell RNA-seq data.

## Usage

``` r
RunDoubletCalling(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  db_method = "scDblFinder",
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

- db_method:

  Method used for doublet-calling. Can be one of `"scDblFinder"`,
  `"Scrublet"`, `"DoubletDetection"`, `"scds_cxds"`, `"scds_bcds"`,
  `"scds_hybrid"`.

- ...:

  Additional arguments to be passed to the corresponding doublet-calling
  method.

## Value

Returns a Seurat object with the doublet prediction results and
prediction scores stored in the meta.data.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-03 09:55:51] Start standard processing workflow...
#> ℹ [2026-04-03 09:55:51] Checking a list of <Seurat>...
#> ! [2026-04-03 09:55:51] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 09:55:51] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 09:55:54] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 09:55:54] Use the separate HVF from `srt_list`
#> ℹ [2026-04-03 09:55:54] Number of available HVF: 2000
#> ℹ [2026-04-03 09:55:54] Finished check
#> ℹ [2026-04-03 09:55:55] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-03 09:55:55] Perform pca linear dimension reduction
#> ℹ [2026-04-03 09:55:56] Use stored estimated dimensions 1:12 for Standardpca
#> ℹ [2026-04-03 09:55:56] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-03 09:55:56] Reorder clusters...
#> ℹ [2026-04-03 09:55:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-03 09:55:56] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-03 09:55:56] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ℹ [2026-04-03 09:56:01] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ✔ [2026-04-03 09:56:06] Standard processing workflow completed
pancreas_sub <- RunDoubletCalling(
  pancreas_sub,
  db_method = "scDblFinder"
)
#> ℹ [2026-04-03 09:56:07] Data type is raw counts
#> ℹ [2026-04-03 09:56:07] Data type is raw counts
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘scale.data’ is empty
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.scDblFinder_class"
)
#> Error in DefaultReduction(srt, pattern = reduction): Unable to find any reductions

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.scDblFinder_score"
)
#> Error in DefaultReduction(srt, pattern = reduction): Unable to find any reductions
```
