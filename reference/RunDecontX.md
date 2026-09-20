# Ambient RNA decontamination with decontX

Ambient RNA decontamination with decontX

## Usage

``` r
RunDecontX(
  object,
  assay = "RNA",
  group.by = NULL,
  batch = NULL,
  background = NULL,
  background_assay = NULL,
  bg_batch = NULL,
  assay_name = "decontXcounts",
  store_assay = TRUE,
  round_counts = FALSE,
  data_type = NULL,
  seed = 11,
  ...,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay to decontaminate.

- group.by, batch:

  Cell cluster and batch labels passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).
  Column name, cell-aligned vector, or `NULL`.

- background:

  Background / empty-droplet input: a `Seurat` object,
  `SingleCellExperiment`, or count matrix.

- background_assay:

  Assay used when `background` is a `Seurat` or `SingleCellExperiment`.
  `NULL` uses `assay` (Seurat) or `"counts"` (SCE).

- bg_batch:

  Batch labels for `background`.

- assay_name, store_assay, round_counts:

  Store rounded decontaminated counts as a new assay.

- data_type:

  Optional
  [`CheckDataType()`](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- seed:

  Random seed.

- ...:

  Passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).

- verbose:

  Whether to print messages.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with decontX contamination in `meta.data` and optional
decontaminated counts in a new assay.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 22:13:30] Start standard processing workflow...
#> ℹ [2026-09-20 22:13:30] Checking a list of <Seurat>...
#> ! [2026-09-20 22:13:31] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 22:13:31] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:13:31] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:13:31] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:13:31] Number of available HVF: 2000
#> ℹ [2026-09-20 22:13:31] Finished check
#> ℹ [2026-09-20 22:13:31] Perform `ScaleData()`
#> ℹ [2026-09-20 22:13:31] Perform pca linear dimension reduction
#> ℹ [2026-09-20 22:13:31] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 22:13:31] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 22:13:32] Reorder clusters...
#> ℹ [2026-09-20 22:13:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:13:32] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 22:13:38] Standard processing workflow completed
pancreas_sub <- RunDecontX(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-09-20 22:13:38] Running decontX
#> ℹ [2026-09-20 22:13:39] Data type is raw counts
#> ℹ [2026-09-20 22:13:52] decontX contamination (median/mean/max): 0.0272 / 0.0875 / 0.6737
#> ℹ [2026-09-20 22:13:52] decontX assay stored as decontXcounts
#> ✔ [2026-09-20 22:13:52] decontX decontamination completed

FeatureStatPlot(
  pancreas_sub,
  stat.by = "decontX_contamination"
)


FeatureDimPlot(
  pancreas_sub,
  features = "decontX_contamination"
)
```
