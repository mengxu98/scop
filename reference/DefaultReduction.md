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
#> ℹ [2026-01-30 16:31:28] Start standard scop workflow...
#> ℹ [2026-01-30 16:31:29] Checking a list of <Seurat>...
#> ! [2026-01-30 16:31:29] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 16:31:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 16:31:31] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 16:31:31] Use the separate HVF from srt_list
#> ℹ [2026-01-30 16:31:32] Number of available HVF: 2000
#> ℹ [2026-01-30 16:31:32] Finished check
#> ℹ [2026-01-30 16:31:32] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 16:31:32] Perform pca linear dimension reduction
#> ℹ [2026-01-30 16:31:33] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 16:31:33] Reorder clusters...
#> ℹ [2026-01-30 16:31:33] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 16:31:33] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 16:31:37] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 16:31:40] Run scop standard workflow completed
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
