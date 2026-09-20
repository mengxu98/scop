# Find the default reduction name in a Seurat object

Find the default reduction name in a Seurat object

## Usage

``` r
DefaultReduction(
  object,
  pattern = NULL,
  min_dim = 2,
  max_distance = 0.1,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- pattern:

  Character string containing a regular expression to search for.

- min_dim:

  Minimum dimension threshold.

- max_distance:

  Maximum distance allowed for a match.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

Default reduction name.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 21:31:54] Start standard processing workflow...
#> ℹ [2026-09-20 21:31:54] Checking a list of <Seurat>...
#> ! [2026-09-20 21:31:54] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 21:31:54] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:31:54] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:31:55] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 21:31:55] Number of available HVF: 2000
#> ℹ [2026-09-20 21:31:55] Finished check
#> ℹ [2026-09-20 21:31:55] Perform `ScaleData()`
#> ℹ [2026-09-20 21:31:55] Perform pca linear dimension reduction
#> ℹ [2026-09-20 21:31:55] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 21:31:55] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 21:31:55] Reorder clusters...
#> ℹ [2026-09-20 21:31:55] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 21:31:55] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 21:32:01] Standard processing workflow completed
names(pancreas_sub@reductions)
#> [1] "Standardpca"       "StandardpcaUMAP2D" "StandardUMAP2D"   

DefaultReduction(pancreas_sub)
#> [1] "StandardUMAP2D"

DefaultReduction(pancreas_sub, pattern = "pca")
#> [1] "Standardpca"

DefaultReduction(pancreas_sub, pattern = "umap")
#> [1] "StandardUMAP2D"
```
