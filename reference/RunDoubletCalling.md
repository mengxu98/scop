# Doublet calling

Doublet calling

## Usage

``` r
RunDoubletCalling(
  object,
  assay = "RNA",
  db_rate = ncol(object)/1000 * 0.01,
  db_method = "scDblFinder",
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

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with doublet class and score columns in `meta.data`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-13 22:16:58] Start standard processing workflow...
#> ℹ [2026-09-13 22:16:58] Checking a list of <Seurat>...
#> ! [2026-09-13 22:16:58] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-13 22:16:58] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:16:59] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:16:59] Use the separate HVF from `srt_list`
#> ℹ [2026-09-13 22:16:59] Number of available HVF: 2000
#> ℹ [2026-09-13 22:16:59] Finished check
#> ℹ [2026-09-13 22:16:59] Perform `ScaleData()`
#> ℹ [2026-09-13 22:16:59] Perform pca linear dimension reduction
#> ℹ [2026-09-13 22:16:59] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-13 22:16:59] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-13 22:17:00] Reorder clusters...
#> ℹ [2026-09-13 22:17:00] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-13 22:17:00] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-13 22:17:06] Standard processing workflow completed
pancreas_sub <- RunDoubletCalling(
  pancreas_sub,
  db_method = "scDblFinder"
)
#> ℹ [2026-09-13 22:17:06] Data type is raw counts
#> ℹ [2026-09-13 22:17:06] Running scDblFinder
#> ℹ [2026-09-13 22:17:07] Data type is raw counts
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
