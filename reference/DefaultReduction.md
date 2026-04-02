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
#> ℹ [2026-04-02 15:30:33] Start standard processing workflow...
#> ℹ [2026-04-02 15:30:34] Checking a list of <Seurat>...
#> ! [2026-04-02 15:30:34] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 15:30:34] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:30:36] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 15:30:36] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 15:30:36] Number of available HVF: 2000
#> ℹ [2026-04-02 15:30:36] Finished check
#> ℹ [2026-04-02 15:30:37] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 15:30:38] Perform pca linear dimension reduction
#> ℹ [2026-04-02 15:30:42] Use stored estimated dimensions 1:50 for Standardpca
#> Warning: Caught FutureLaunchError. Canceling all iterations ...
#> ! [2026-04-02 15:30:42] <FutureLaunchError: Caught an unexpected error of class FutureLaunchError when trying to launch future (‘future_lapply-1’) on backend of class SequentialFutureBackend. The reason was: future::evalFuture() failed on runnervmrg6be (pid 85355) at 2026-04-02T15:30:42. Using package 'future' v1.70.0. Possible other reasons: Failed to attach one or more future-backend packages: there is no package called ‘future’ [future <unnamed>; on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>] [future ‘future_lapply-1’ (4a75d434f7a9a2903adedbeee3372830-14); on 4a75d434f7a9a2903adedbeee3372830@runnervmrg6be<85355>]>
#> !                       
#> !                       Occurred on: 4a75d434f7a9a2903adedbeee3372830 [runnervmrg6be; pid 85355]
#> !                       Future: 4a75d434f7a9a2903adedbeee3372830-14 (‘future_lapply-1’)
#> !                       
#> !                       DEBUG: BEGIN TROUBLESHOOTING HELP
#> !                       SequentialFuture:
#> !                       Label: ‘future_lapply-1’
#> !                       Expression:
#> Error in glue(str, .envir = .envir, .transformer = transformer, .cli = TRUE,     .trim = .trim): Expecting '}'
names(pancreas_sub@reductions)
#> character(0)

DefaultReduction(pancreas_sub)
#> Error in DefaultReduction(pancreas_sub): Unable to find any reductions

DefaultReduction(pancreas_sub, pattern = "pca")
#> Error in DefaultReduction(pancreas_sub, pattern = "pca"): Unable to find any reductions

DefaultReduction(pancreas_sub, pattern = "umap")
#> Error in DefaultReduction(pancreas_sub, pattern = "umap"): Unable to find any reductions
```
