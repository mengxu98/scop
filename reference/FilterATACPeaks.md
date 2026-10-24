# Filter an ATAC peak universe with explicit criteria

Select peaks by raw-count support, interval width, chromosome, and
optional blacklist overlap, then create a clean counts-only assay.

## Usage

``` r
FilterATACPeaks(
  object,
  assay = NULL,
  new.assay = "filtered_peaks",
  min.cells = 0L,
  min.counts = 0,
  min.width = 1L,
  max.width = Inf,
  seqlevels = NULL,
  blacklist = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with a validated `ChromatinAssay`.

- assay:

  Source chromatin assay; defaults to the active assay.

- new.assay:

  An unused assay name, default `"filtered_peaks"`.

- min.cells:

  Minimum number of cells with a nonzero raw count, inclusive.

- min.counts:

  Minimum sum of raw counts across cells, inclusive.

- min.width, max.width:

  Inclusive integer interval-width limits in bases. The default retains
  every positive width; `max.width` may be `Inf`.

- seqlevels:

  Optional unique chromosome names to retain. Names must occur in the
  source reference; chromosome styles are never guessed.

- blacklist:

  Optional `GRanges` of excluded 1-based closed intervals. Any overlap
  excludes a peak regardless of strand. The blacklist must declare the
  same genome label as the source assay. Known chromosome lengths must
  agree and interval ends must not exceed them. Chromosomes absent from
  the source assay are retained as supplied; names are never converted.

- verbose:

  Whether to print progress messages.

## Value

The Seurat object with the filtered `ChromatinAssay`. All cells and
retained peak order/counts are preserved. The new assay retains
reference, gene annotation and fragment references, but not
normalized/scaled data, variable features, motifs, links or downstream
results. Native initialization puts raw counts in the data layer. Source
assays, active assay and existing reductions are unchanged.
`Misc(object[[new.assay]], "scop_peak_filter")` stores one row per
source peak with each criterion and the final keep flag. At least two
peaks and two cells must remain for native assay construction.

## See also

[`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md),
[`AnnotateATACPeaks()`](https://mengxu98.github.io/scop/reference/AnnotateATACPeaks.md)

## Examples

``` r
counts <- Matrix::sparseMatrix(
  i = c(1, 2, 1), j = c(1, 2, 2), x = c(2, 3, 1), dims = c(2, 2),
  dimnames = list(c("p1", "p2"), c("cell1", "cell2"))
)
ranges <- GenomicRanges::GRanges(c("chr1:1-50", "chr1:101-150"))
object <- ReadATACData(counts, ranges = ranges, genome = "hg38")
filtered <- FilterATACPeaks(object, min.cells = 1, verbose = FALSE)
rownames(filtered[["filtered_peaks"]])
#> [1] "chr1-1-50"    "chr1-101-150"
```
