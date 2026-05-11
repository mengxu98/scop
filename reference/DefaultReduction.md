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
#> ℹ [2026-05-11 14:50:59] Start standard processing workflow...
#> ℹ [2026-05-11 14:51:00] Checking a list of <Seurat>...
#> ! [2026-05-11 14:51:00] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 14:51:00] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 14:51:01] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 14:51:02] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 14:51:02] Number of available HVF: 2000
#> ℹ [2026-05-11 14:51:02] Finished check
#> ℹ [2026-05-11 14:51:02] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 14:51:02] Perform pca linear dimension reduction
#> ℹ [2026-05-11 14:51:03] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-11 14:51:03] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-11 14:51:03] Reorder clusters...
#> ℹ [2026-05-11 14:51:03] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-11 14:51:03] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-11 14:51:03] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-11 14:51:07] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-11 14:51:10] Standard processing workflow completed
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
