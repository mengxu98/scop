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
#> ◌ [2026-11-21 11:07:44] Preparing python environment...
#> ℹ [2026-11-21 11:07:44] Environment name: scop_env and python version: 3.10-1
#> ℹ [2026-11-21 11:07:44] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-11-21 11:07:44] Using micromamba executable: /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#> ℹ [2026-11-21 11:07:44] Using existing environment: /home/runner/work/_temp/pkgdown-python/envs/scop_env
#> ℹ [2026-11-21 11:07:44] 14 requirements to check, 3 conda requirements and 11 pip requirements
#> ℹ [2026-11-21 11:07:44] Checking 3 conda requirements
#> ℹ [2026-11-21 11:07:45] Checking 13 requirements using uv
#> + /home/runner/work/_temp/pkgdown-python/envs/scop_env/bin/python3.10 -m pip check
#> ✔ [2026-11-21 11:07:47] Python environment ready
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
#> ℹ [2026-11-21 11:07:48] Start standard processing workflow...
#> ℹ [2026-11-21 11:07:48] Checking a list of <Seurat>...
#> ! [2026-11-21 11:07:48] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-11-21 11:07:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-11-21 11:07:48] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-11-21 11:07:48] Use the separate HVF from `srt_list`
#> ℹ [2026-11-21 11:07:48] Number of available HVF: 2000
#> ℹ [2026-11-21 11:07:48] Finished check
#> ℹ [2026-11-21 11:07:48] Perform `ScaleData()`
#> ℹ [2026-11-21 11:07:48] Perform pca linear dimension reduction
#> ℹ [2026-11-21 11:07:49] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-11-21 11:07:49] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-11-21 11:07:49] Reorder clusters...
#> ℹ [2026-11-21 11:07:49] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-21 11:07:49] Perform umap nonlinear dimension reduction
#> ✔ [2026-11-21 11:07:57] Standard processing workflow completed
pancreas_sub <- RunDoubletDetection(pancreas_sub)
#> ℹ [2026-11-21 11:07:57] Running DoubletDetection
#> ✔ [2026-11-21 11:07:58] Python environment completed
#> ✔                       Until next loading, the environment will be cached
#> ℹ [2026-11-21 11:07:58] Data type is raw counts
#> ℹ [2026-11-21 11:07:58] Checking 1 package in environment: scop_env
#> ℹ [2026-11-21 11:07:58] Retrieving package list for environment: scop_env
#> ℹ [2026-11-21 11:07:59] Found 102 packages installed
#> ✔ [2026-11-21 11:07:59] doubletdetection == 4.3.0.post1
#> ℹ [2026-11-21 11:07:59] Checking 1 package in environment: scop_env
#> ℹ [2026-11-21 11:07:59] Retrieving package list for environment: scop_env
#> ℹ [2026-11-21 11:08:00] Found 102 packages installed
#> ✔ [2026-11-21 11:08:00] louvain == 0.8.2
#> ! [2026-11-21 11:08:17] Found 1 cells with unexpected DoubletDetection labels; setting to "doublet"
#> ✔ [2026-11-21 11:08:17] DoubletDetection doublet calling completed
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
