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
#> ℹ [2026-01-29 13:17:36] Start standard scop workflow...
#> ℹ [2026-01-29 13:17:36] Checking a list of <Seurat>...
#> ! [2026-01-29 13:17:36] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-29 13:17:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-29 13:17:38] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-29 13:17:39] Use the separate HVF from srt_list
#> ℹ [2026-01-29 13:17:39] Number of available HVF: 2000
#> ℹ [2026-01-29 13:17:39] Finished check
#> ℹ [2026-01-29 13:17:39] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-29 13:17:40] Perform pca linear dimension reduction
#> ℹ [2026-01-29 13:17:41] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-29 13:17:41] Reorder clusters...
#> ℹ [2026-01-29 13:17:41] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-29 13:17:41] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-29 13:17:45] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-29 13:17:49] Run scop standard workflow completed
pancreas_sub <- RunDoubletCalling(
  pancreas_sub,
  db_method = "scDblFinder"
)
#> ℹ [2026-01-29 13:17:50] Data type is raw counts
#> ℹ [2026-01-29 13:17:50] Data type is raw counts
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
