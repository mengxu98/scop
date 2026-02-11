# Find the default reduction name in a Seurat object

Find the default reduction name in a Seurat object

## Usage

``` r
DefaultReduction(srt, pattern = NULL, min_dim = 2, max_distance = 0.1)
```

## Arguments

- srt:

  A Seurat object.

- pattern:

  Character string containing a regular expression to search for.

- min_dim:

  Minimum dimension threshold.

- max_distance:

  Maximum distance allowed for a match.

## Value

Default reduction name.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-02-11 03:16:48] Start standard scop workflow...
#> ℹ [2026-02-11 03:16:49] Checking a list of <Seurat>...
#> ! [2026-02-11 03:16:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-02-11 03:16:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-02-11 03:16:50] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-02-11 03:16:51] Use the separate HVF from srt_list
#> ℹ [2026-02-11 03:16:51] Number of available HVF: 2000
#> ℹ [2026-02-11 03:16:51] Finished check
#> ℹ [2026-02-11 03:16:51] Perform `Seurat::ScaleData()`
#> ℹ [2026-02-11 03:16:52] Perform pca linear dimension reduction
#> ℹ [2026-02-11 03:16:53] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-02-11 03:16:53] Reorder clusters...
#> ℹ [2026-02-11 03:16:53] Perform umap nonlinear dimension reduction
#> ℹ [2026-02-11 03:16:53] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-02-11 03:16:56] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-02-11 03:16:59] Run scop standard workflow completed
names(pancreas_sub@reductions)
#> [1] "Standardpca"       "StandardpcaUMAP2D" "StandardpcaUMAP3D"
#> [4] "StandardUMAP2D"    "StandardUMAP3D"   

DefaultReduction(pancreas_sub)
#> [1] "StandardUMAP2D"

DefaultReduction(pancreas_sub, pattern = "pca")
#> [1] "Standardpca"

DefaultReduction(pancreas_sub, pattern = "umap")
#> [1] "StandardUMAP2D"
```
