# Run doublet-calling with DoubletDetection

Run doublet-calling with DoubletDetection

## Usage

``` r
RunDoubletDetection(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  cores = 1,
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

- cores:

  The number of CPU cores to use for `doubletdetection`.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- ...:

  Additional arguments to be passed to
  [doubletdetection.BoostClassifier](https://github.com/JonathanShor/DoubletDetection).

- verbose:

  Whether to print messages.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:05:47] Start standard processing workflow...
#> ℹ [2026-09-06 22:05:48] Checking a list of <Seurat>...
#> ! [2026-09-06 22:05:48] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:05:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:05:48] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:05:48] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:05:48] Number of available HVF: 2000
#> ℹ [2026-09-06 22:05:48] Finished check
#> ℹ [2026-09-06 22:05:48] Perform `ScaleData()`
#> ℹ [2026-09-06 22:05:48] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:05:48] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:05:49] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:05:49] Reorder clusters...
#> ℹ [2026-09-06 22:05:49] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:05:49] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:05:56] Standard processing workflow completed
pancreas_sub <- RunDoubletDetection(pancreas_sub)
#> ℹ [2026-09-06 22:05:56] Running DoubletDetection
#> Error in if (existing_minor %in% c("3.10", "3.11", "3.12")) {    version <- paste0(existing_minor, "-1")}: argument is of length zero
CellDimPlot(
  pancreas_sub,
  reduction = "umap",
  group.by = "db.DoubletDetection_class"
)
#> Error in CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.DoubletDetection_class"): "db.DoubletDetection_class" is not in the meta.data of srt object

FeatureDimPlot(
  pancreas_sub,
  reduction = "umap",
  features = "db.DoubletDetection_score"
)
#> ! [2026-09-06 22:06:22] "db.DoubletDetection_score" are not in the features of <Seurat>
#> Error in FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.DoubletDetection_score"): There are no valid features present.
```
