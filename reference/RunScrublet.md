# Run doublet-calling with Scrublet

Run doublet-calling with Scrublet

## Usage

``` r
RunScrublet(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE
)
```

## Arguments

- srt:

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

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:33:15] Start standard processing workflow...
#> ℹ [2026-09-06 22:33:16] Checking a list of <Seurat>...
#> ! [2026-09-06 22:33:16] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:33:16] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:33:16] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:33:16] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:33:16] Number of available HVF: 2000
#> ℹ [2026-09-06 22:33:16] Finished check
#> ℹ [2026-09-06 22:33:16] Perform `ScaleData()`
#> ℹ [2026-09-06 22:33:16] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:33:17] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:33:17] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:33:17] Reorder clusters...
#> ℹ [2026-09-06 22:33:17] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:33:17] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:33:25] Standard processing workflow completed
pancreas_sub <- RunScrublet(pancreas_sub)
#> ℹ [2026-09-06 22:33:25] Running Scrublet
#> Error in if (existing_minor %in% c("3.10", "3.11", "3.12")) {    version <- paste0(existing_minor, "-1")}: argument is of length zero
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.Scrublet_class"
)
#> Error in CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.Scrublet_class"): "db.Scrublet_class" is not in the meta.data of srt object

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.Scrublet_score"
)
#> ! [2026-09-06 22:33:28] "db.Scrublet_score" are not in the features of <Seurat>
#> Error in FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.Scrublet_score"): There are no valid features present.
```
