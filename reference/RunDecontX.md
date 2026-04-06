# Run ambient RNA decontamination with decontX

Run ambient RNA decontamination with decontX

## Usage

``` r
RunDecontX(
  srt,
  assay = "RNA",
  group.by = NULL,
  batch = NULL,
  background = NULL,
  background_assay = NULL,
  bg_batch = NULL,
  assay_name = "decontXcounts",
  store_assay = TRUE,
  round_counts = FALSE,
  seed = 11,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for decontamination. Default is
  `"RNA"`.

- group.by:

  Cell cluster labels passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).
  Can be `NULL`, a meta.data column name, or a vector aligned to cells.
  Default is `NULL`.

- batch:

  Batch labels passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).
  Can be `NULL`, a meta.data column name, or a vector aligned to cells.
  Default is `NULL`.

- background:

  Optional background / empty-droplet input passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).
  Can be a `Seurat` object, `SingleCellExperiment`, or count matrix.
  Default is `NULL`.

- background_assay:

  Assay name used when `background` is a `Seurat` object or
  `SingleCellExperiment`. Default is `NULL`, which falls back to `assay`
  for `Seurat` background and `"counts"` for `SingleCellExperiment`
  background.

- bg_batch:

  Batch labels for `background` passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).
  Can be `NULL`, a metadata column name, or a vector aligned to the
  background droplets. Default is `NULL`.

- assay_name:

  Name of the assay used to store decontaminated counts. Default is
  `"decontXcounts"`.

- store_assay:

  Whether to store decontaminated counts as a new assay. Default is
  `TRUE`.

- round_counts:

  Whether to round decontaminated counts before creating the assay.
  Default is `FALSE`.

- seed:

  Random seed for reproducibility. Default is `11`.

- ...:

  Additional arguments passed to
  [`decontX::decontX()`](https://rdrr.io/pkg/decontX/man/decontX.html).

## Value

Returns a Seurat object with decontX contamination estimates stored in
the meta.data, and optional decontaminated counts stored in a new assay.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-06 11:15:31] Start standard processing workflow...
#> ℹ [2026-04-06 11:15:31] Checking a list of <Seurat>...
#> ! [2026-04-06 11:15:31] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 11:15:31] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 11:15:34] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 11:15:35] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 11:15:35] Number of available HVF: 2000
#> ℹ [2026-04-06 11:15:35] Finished check
#> ℹ [2026-04-06 11:15:35] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-06 11:15:36] Perform pca linear dimension reduction
#> ℹ [2026-04-06 11:15:36] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-06 11:15:37] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-06 11:15:37] Reorder clusters...
#> ℹ [2026-04-06 11:15:37] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 11:15:37] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-06 11:15:37] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-06 11:15:41] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-06 11:15:46] Standard processing workflow completed
pancreas_sub <- RunDecontX(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-04-06 11:15:47] Data type is raw counts
#> ℹ [2026-04-06 11:15:47] Running decontX
#> ℹ [2026-04-06 11:16:00] decontX contamination (median/mean/max): 0.0272 / 0.0875 / 0.6737
#> ℹ [2026-04-06 11:16:00] decontX assay stored as decontXcounts
#> ✔ [2026-04-06 11:16:00] decontX decontamination completed

FeatureStatPlot(
  pancreas_sub,
  stat.by = "decontX_contamination"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureDimPlot(
  pancreas_sub,
  features = "decontX_contamination"
)
```
