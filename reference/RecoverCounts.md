# Attempt to recover raw counts from the normalized matrix

Attempt to recover raw counts from the normalized matrix

## Usage

``` r
RecoverCounts(
  object,
  assay = NULL,
  trans = c("expm1", "exp", "none"),
  min_count = c(1, 2, 3),
  tolerance = 0.1,
  sf = NULL,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay to use. `NULL` uses the default assay.

- trans:

  The transformation function to applied when data is presumed to be
  log-normalized.

- min_count:

  Minimum UMI count of genes.

- tolerance:

  When recovering the raw counts, the nCount of each cell is
  theoretically calculated as an integer. However, due to decimal point
  preservation during normalization, the calculated nCount is usually a
  floating point number close to the integer. The tolerance is its
  difference from the integer.

- sf:

  Set the scaling factor manually.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## See also

[CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 22:00:15] Start standard processing workflow...
#> ℹ [2026-09-20 22:00:15] Checking a list of <Seurat>...
#> ! [2026-09-20 22:00:15] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 22:00:15] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:00:15] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:00:15] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:00:15] Number of available HVF: 2000
#> ℹ [2026-09-20 22:00:15] Finished check
#> ℹ [2026-09-20 22:00:15] Perform `ScaleData()`
#> ℹ [2026-09-20 22:00:15] Perform pca linear dimension reduction
#> ℹ [2026-09-20 22:00:15] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 22:00:16] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 22:00:16] Reorder clusters...
#> ℹ [2026-09-20 22:00:16] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:00:16] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 22:00:22] Standard processing workflow completed
raw_counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)

# Normalized the data
pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#> Normalizing layer: counts

# Now replace counts with the log-normalized data matrix
data <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "data"
)
new_pancreas_sub <- SeuratObject::SetAssayData(
  object = pancreas_sub,
  layer = "counts",
  new.data = data,
  assay = "RNA"
)
# Recover the counts and compare with the raw counts matrix
pancreas_sub <- RecoverCounts(new_pancreas_sub)
#> ℹ [2026-09-20 22:00:24] Data type is log-normalized
#> ℹ [2026-09-20 22:00:24] The data is presumed to be log-normalized
#> ℹ [2026-09-20 22:00:24] Perform "expm1" on the raw data
new_counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)
identical(raw_counts, new_counts)
#> [1] TRUE
```
