# Attempt to recover raw counts from the normalized matrix

Attempt to recover raw counts from the normalized matrix

## Usage

``` r
RecoverCounts(
  srt,
  assay = NULL,
  trans = c("expm1", "exp", "none"),
  min_count = c(1, 2, 3),
  tolerance = 0.1,
  sf = NULL,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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
  difference from the integer. Default is `0.1`

- sf:

  Set the scaling factor manually.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-11 15:34:50] Start standard processing workflow...
#> ℹ [2026-05-11 15:34:51] Checking a list of <Seurat>...
#> ! [2026-05-11 15:34:51] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 15:34:51] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 15:34:52] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 15:34:53] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 15:34:53] Number of available HVF: 2000
#> ℹ [2026-05-11 15:34:53] Finished check
#> ℹ [2026-05-11 15:34:53] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 15:34:53] Perform pca linear dimension reduction
#> ℹ [2026-05-11 15:34:54] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-11 15:34:54] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-11 15:34:54] Reorder clusters...
#> ℹ [2026-05-11 15:34:54] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-11 15:34:54] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-11 15:34:54] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-11 15:34:59] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-11 15:35:03] Standard processing workflow completed
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
#> ℹ [2026-05-11 15:35:05] Data type is log-normalized
#> ℹ [2026-05-11 15:35:05] The data is presumed to be log-normalized
#> ℹ [2026-05-11 15:35:05] Perform "expm1" on the raw data
new_counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)
identical(raw_counts, new_counts)
#> [1] TRUE
```
