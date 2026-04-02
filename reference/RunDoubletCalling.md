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
#> ℹ [2026-04-02 16:39:06] Start standard processing workflow...
#> ℹ [2026-04-02 16:39:06] Checking a list of <Seurat>...
#> ! [2026-04-02 16:39:06] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:39:06] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:39:08] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:39:08] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:39:09] Number of available HVF: 2000
#> ℹ [2026-04-02 16:39:09] Finished check
#> ℹ [2026-04-02 16:39:09] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:39:09] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:39:13] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:39:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:39:13] Reorder clusters...
#> ℹ [2026-04-02 16:39:14] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:39:14] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:39:14] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:39:14] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:39:14] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:39:17] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:39:20] Standard processing workflow completed
pancreas_sub <- RunDoubletCalling(
  pancreas_sub,
  db_method = "scDblFinder"
)
#> ℹ [2026-04-02 16:39:20] Data type is raw counts
#> ℹ [2026-04-02 16:39:20] Data type is raw counts
#> Warning: Layer ‘data’ is empty
#> Warning: Layer ‘scale.data’ is empty
#> Error in tryCatchOne(expr, names, parentenv, handlers[[1L]]): <packageNotFoundError in loadNamespace(x): there is no package called
#> ‘scDblFinder’>
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.scDblFinder_class"
)
#> Error in CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.scDblFinder_class"): "db.scDblFinder_class" is not in the meta.data of srt object

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.scDblFinder_score"
)
#> ! [2026-04-02 16:39:27] "db.scDblFinder_score" are not in the features of <Seurat>
#> Error in FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.scDblFinder_score"): There are no valid features present.
```
