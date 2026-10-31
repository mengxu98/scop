# Run doublet-calling with Scrublet

Run doublet-calling with Scrublet

## Usage

``` r
RunScrublet(
  object,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay used for doublet calling.

- db_rate:

  Expected doublet rate.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- ...:

  Additional arguments to be passed to
  [scrublet.Scrublet](https://github.com/swolock/scrublet).

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Examples

``` r
PrepareEnv(modules = "scrublet")
#> ◌ [2026-10-31 18:40:46] Preparing python environment...
#> ℹ [2026-10-31 18:40:46] Environment name: scop_env and python version: 3.10-1
#> ℹ [2026-10-31 18:40:46] Auto-detecting conda-compatible environment manager...
#> ℹ [2026-10-31 18:40:46] Using micromamba executable: /home/runner/.cache/R/scop/micromamba/linux-64/bin/micromamba
#> ℹ [2026-10-31 18:40:46] Using existing environment: /home/runner/work/_temp/pkgdown-python/envs/scop_env
#> ℹ [2026-10-31 18:40:46] 13 requirements to check, 3 conda requirements and 10 pip requirements
#> ℹ [2026-10-31 18:40:46] Checking 3 conda requirements
#> ℹ [2026-10-31 18:40:47] Checking 12 requirements using uv
#> + /home/runner/work/_temp/pkgdown-python/envs/scop_env/bin/python3.10 -m pip check
#> ✔ [2026-10-31 18:40:50] Python environment ready
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
#> ℹ [2026-10-31 18:40:50] Start standard processing workflow...
#> ℹ [2026-10-31 18:40:50] Checking a list of <Seurat>...
#> ! [2026-10-31 18:40:50] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:40:50] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-31 18:40:50] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-31 18:40:50] Use the separate HVF from `srt_list`
#> ℹ [2026-10-31 18:40:50] Number of available HVF: 2000
#> ℹ [2026-10-31 18:40:50] Finished check
#> ℹ [2026-10-31 18:40:50] Perform `ScaleData()`
#> ℹ [2026-10-31 18:40:50] Perform pca linear dimension reduction
#> ℹ [2026-10-31 18:40:51] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-31 18:40:51] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-31 18:40:51] Reorder clusters...
#> ℹ [2026-10-31 18:40:51] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-31 18:40:51] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-31 18:41:01] Standard processing workflow completed
pancreas_sub <- RunScrublet(pancreas_sub)
#> ℹ [2026-10-31 18:41:01] Running Scrublet
#> ✔ [2026-10-31 18:41:01] Python environment completed
#> ✔                       Until next loading, the environment will be cached
#> ℹ [2026-10-31 18:41:01] Data type is raw counts
#> ✔ [2026-10-31 18:41:03] Scrublet doublet calling completed
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.Scrublet_class"
)


FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.Scrublet_score"
)
```
