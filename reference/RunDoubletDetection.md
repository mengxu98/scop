# Run doublet-calling with DoubletDetection

Run doublet-calling with DoubletDetection

## Usage

``` r
RunDoubletDetection(
  object,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  cores = 1,
  data_type = NULL,
  verbose = TRUE,
  srt = NULL,
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay used for doublet calling.

- db_rate:

  Expected doublet rate.

- cores:

  The number of CPU cores to use for `doubletdetection`.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- ...:

  Arguments passed to other functions.

## Examples

``` r
PrepareEnv(modules = "doubletdetection")
#> ◌ [2026-11-12 03:51:10] Preparing python environment...
#> ℹ [2026-11-12 03:51:10] Environment name: scop_env and python version: 3.10-1
#> ℹ [2026-11-12 03:51:10] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-11-12 03:51:10] Using micromamba executable: /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#> ℹ [2026-11-12 03:51:10] Using existing environment: /home/runner/work/_temp/pkgdown-python/envs/scop_env
#> ℹ [2026-11-12 03:51:10] 14 requirements to check, 3 conda requirements and 11 pip requirements
#> ℹ [2026-11-12 03:51:10] Checking 3 conda requirements
#> ℹ [2026-11-12 03:51:11] Checking 13 requirements using uv
#> + /home/runner/work/_temp/pkgdown-python/envs/scop_env/bin/python3.10 -m pip check
#> ✔ [2026-11-12 03:51:13] Python environment ready
#> ✔                       Until next loading, the environment will be cached
#> Environment manager config:
#>  Manager:       micromamba
#>  Executable:    /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#>  Environment:   /home/runner/work/_temp/pkgdown-python/envs/scop_env
#> Python config:
#>  python:        /home/runner/work/_temp/pkgdown-python/envs/scop_env/bin/python3.10
#>  version:       Python 3.10.1
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-11-12 03:51:14] Start standard processing workflow...
#> ℹ [2026-11-12 03:51:14] Checking a list of <Seurat>...
#> ! [2026-11-12 03:51:14] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-11-12 03:51:14] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-11-12 03:51:14] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-11-12 03:51:14] Use the separate HVF from `srt_list`
#> ℹ [2026-11-12 03:51:14] Number of available HVF: 2000
#> ℹ [2026-11-12 03:51:14] Finished check
#> ℹ [2026-11-12 03:51:14] Perform `ScaleData()`
#> ℹ [2026-11-12 03:51:14] Perform pca linear dimension reduction
#> ℹ [2026-11-12 03:51:15] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-11-12 03:51:15] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-11-12 03:51:15] Reorder clusters...
#> ℹ [2026-11-12 03:51:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-12 03:51:15] Perform umap nonlinear dimension reduction
#> ✔ [2026-11-12 03:51:24] Standard processing workflow completed
pancreas_sub <- RunDoubletDetection(pancreas_sub)
#> ℹ [2026-11-12 03:51:24] Running DoubletDetection
#> ✔ [2026-11-12 03:51:24] Python environment completed
#> ✔                       Until next loading, the environment will be cached
#> ℹ [2026-11-12 03:51:24] Data type is raw counts
#> ℹ [2026-11-12 03:51:25] Checking 1 package in environment: scop_env
#> ℹ [2026-11-12 03:51:25] Retrieving package list for environment: scop_env
#> ℹ [2026-11-12 03:51:25] Found 102 packages installed
#> ✔ [2026-11-12 03:51:25] doubletdetection == 4.3.0.post1
#> ℹ [2026-11-12 03:51:25] Checking 1 package in environment: scop_env
#> ℹ [2026-11-12 03:51:26] Retrieving package list for environment: scop_env
#> ℹ [2026-11-12 03:51:26] Found 102 packages installed
#> ✔ [2026-11-12 03:51:26] louvain == 0.8.2
#> ! [2026-11-12 03:51:45] Found 1 cells with unexpected DoubletDetection labels; setting to "doublet"
#> ✔ [2026-11-12 03:51:45] DoubletDetection doublet calling completed
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.DoubletDetection_class"
)


FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.DoubletDetection_score"
)
```
