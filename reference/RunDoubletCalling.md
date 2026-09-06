# Doublet calling

Doublet calling

## Usage

``` r
RunDoubletCalling(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  db_method = "scDblFinder",
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

- db_method:

  `"scDblFinder"`, `"Scrublet"`, `"DoubletDetection"`, `"scds_cxds"`,
  `"scds_bcds"`, or `"scds_hybrid"`.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- ...:

  Passed to the selected doublet method.

- verbose:

  Whether to print messages.

## Value

A `Seurat` object with doublet class and score columns in `meta.data`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:04:57] Start standard processing workflow...
#> ℹ [2026-09-06 22:04:57] Checking a list of <Seurat>...
#> ! [2026-09-06 22:04:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:04:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:04:57] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:04:58] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:04:58] Number of available HVF: 2000
#> ℹ [2026-09-06 22:04:58] Finished check
#> ℹ [2026-09-06 22:04:58] Perform `ScaleData()`
#> ℹ [2026-09-06 22:04:58] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:04:58] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:04:58] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:04:58] Reorder clusters...
#> ℹ [2026-09-06 22:04:58] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:04:58] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:05:05] Standard processing workflow completed
pancreas_sub <- RunDoubletCalling(
  pancreas_sub,
  db_method = "scDblFinder"
)
#> ℹ [2026-09-06 22:05:05] Data type is raw counts
#> ℹ [2026-09-06 22:05:05] Running scDblFinder
#> ℹ [2026-09-06 22:05:06] Data type is raw counts
table(pancreas_sub$db.scDblFinder_class)
#> 
#> singlet doublet 
#>     977      23 
head(pancreas_sub@meta.data[, c(
  "db.scDblFinder_class", "db.scDblFinder_score"
)])
#>                  db.scDblFinder_class db.scDblFinder_score
#> AAACCTGAGCCTTGAT              singlet         0.0011855627
#> AAACCTGGTAAGTGGC              singlet         0.0873451009
#> AAACGGGAGATATGGT              singlet         0.0003235179
#> AAACGGGCAAAGAATC              singlet         0.0389233232
#> AAACGGGGTACAGTTC              singlet         0.0113360332
#> AAACGGGTCAGCTCTC              singlet         0.0005065355
```
