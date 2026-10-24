# Export an ATAC object and explicit result/resource tables

Save the complete object, supplied result tables and external resource
inventory in a new directory with file checksums.

## Usage

``` r
ExportATACResults(
  object,
  outdir,
  assay = NULL,
  tables = list(),
  counts = FALSE,
  compress = TRUE
)
```

## Arguments

- object:

  A validated Seurat ChromatinAssay object.

- outdir:

  A new output directory. Existing paths are never overwritten.

- assay:

  Source chromatin assay; defaults to the active assay.

- tables:

  Uniquely named list of data frames, such as complete DAR tests,
  design/contrast tables, motif enrichment or peak annotation. Names
  must be safe filenames composed of letters, numbers, dots, underscores
  or hyphens, unique even when compared case-insensitively.

- counts:

  Whether to additionally write the sparse peak counts as Matrix Market.
  The full object RDS always retains its assays and analysis results.

- compress:

  Compression passed to saveRDS for the complete object.

## Value

Invisibly, a list with the normalized directory, export file manifest
and validation report. Tables are exported as supplied, without
filtering.

## Details

The output contains object.rds, cells.csv, peaks.bed (0-based
half-open), fragments.csv, fragment_cells.csv, check.rds,
sessionInfo.txt, supplied result tables and a file manifest with MD5
checksums. Fragments and reference genomes remain external; no automatic
copying, downloading or path repair occurs. The BED has three columns
and, with counts.mtx and barcodes.tsv when `counts = TRUE`, can be read
by
[`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md)
when the selected assay satisfies its native minimum of two peaks and
two cells. Smaller subsets can still be restored from object.rds.
cells.csv and barcodes.tsv contain the selected assay cells in matrix
column order. object.rds preserves all assays and cells. CSV files
preserve values, not R column classes. The `fragment` column joins
fragment_cells.csv to fragments.csv. The resource table retains inactive
resources, with `active` marking current assay use; `n_cells` is the
attached mapping size, not the exported cell count. Missing fragment
files do not prevent a valid matrix-only export. Stored Signac resource
hashes are labeled as such, and resource availability is only a
snapshot. After moving resources, explicitly update paths with native
Signac UpdatePath and rerun CheckATACData before fragment-dependent
computations. Files are staged in the destination parent and moved to
the requested new directory only after all writes succeed. This does not
modify the object.

## Examples

``` r
counts <- Matrix::sparseMatrix(
  i = c(1, 2, 1), j = c(1, 2, 2), x = c(2, 3, 1), dims = c(2, 2),
  dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
)
ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
destination <- tempfile("atac-export-")
result <- ExportATACResults(object, destination, counts = TRUE)
result$manifest$file
#>  [1] "README.txt"         "barcodes.tsv"       "cells.csv"         
#>  [4] "check.rds"          "counts.mtx"         "fragment_cells.csv"
#>  [7] "fragments.csv"      "object.rds"         "peaks.bed"         
#> [10] "sessionInfo.txt"   
restored <- readRDS(file.path(destination, "object.rds"))
CheckATACData(restored)$valid
#> [1] TRUE
matrix_only <- ReadATACData(file.path(destination, "counts.mtx"),
  peaks.file = file.path(destination, "peaks.bed"),
  barcodes.file = file.path(destination, "barcodes.tsv"), genome = "hg38"
)
CheckATACData(matrix_only)$valid
#> [1] TRUE
unlink(destination, recursive = TRUE)
```
