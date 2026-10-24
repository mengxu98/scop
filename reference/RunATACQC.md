# Run scATAC quality control metrics

Calculate common scATAC QC metrics and optionally filter cells by
thresholds.

## Usage

``` r
RunATACQC(
  object,
  assay = NULL,
  tss.positions = NULL,
  blacklist = NULL,
  fast = TRUE,
  min_pct_reads_in_peaks = NULL,
  min_TSS_enrichment = NULL,
  max_nucleosome_signal = NULL,
  max_blacklist_ratio = NULL,
  verbose = TRUE,
  srt = NULL,
  frip.columns = NULL,
  frip = FALSE,
  frip.chunk.size = 100000L
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay to use. `NULL` uses the default assay.

- tss.positions:

  TSS positions passed to
  [`Signac::TSSEnrichment`](https://stuartlab.org/signac/reference/TSSEnrichment.html).
  An explicit value requests recalculation, even when an imported TSS
  score exists.

- blacklist:

  A `GRanges` blacklist used to compute `blacklist_ratio`.

- fast:

  Whether to use the fast mode in
  [`Signac::TSSEnrichment`](https://stuartlab.org/signac/reference/TSSEnrichment.html).

- min_pct_reads_in_peaks, min_TSS_enrichment, max_nucleosome_signal,
  max_blacklist_ratio:

  Optional thresholds used for filtering cells. Fraction thresholds use
  0–1, not percentages. A requested metric must be present and finite
  for every cell; missing or failed calculations never silently disable
  a requested filter.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- frip.columns:

  Optional named character vector with `in_peaks` and `total` metadata
  column names. Supplying these explicitly asserts that both count the
  same read/fragment unit over the same files, cells and filtering
  rules, and that the numerator refers to the current peak universe.
  Legacy column names are not guessed. Without either this argument or
  `frip = TRUE`, an existing `pct_reads_in_peaks` column is preserved,
  but a new FRiP is not inferred from matrix column sums.

- frip:

  Explicitly recalculate a fragment-record fraction overlapping the
  current peak union. Each record counts once, even across several
  peaks. Duplicate input rows remain separate records; column-five
  support is not used as a weight. Stores `frip` and the legacy 0–1
  alias.

- frip.chunk.size:

  Maximum fragment records read per sequential block.

## Value

A `Seurat` object, retaining input cell order or the ordered subset
passing all requested thresholds. Available metadata include fragment
counts, `nucleosome_signal`, `nucleosome_percentile`, `TSS.enrichment`,
`pct_reads_in_peaks` and `blacklist_ratio`. Direct fragment FRiP also
writes `frip` and a `scop_frip` receipt in the assay's `misc` slot,
recording the peak union, mapped cells, source hashes and per-cell
numerator/denominator. Explicit source-column recalculation also
synchronizes `frip` and clears any previous direct-fragment receipt.
Metrics are calculated only for cells in the selected assay; other cells
and their existing metadata are preserved. Filtering removes only failed
cells from the selected assay; at least one cell in that assay must pass
the requested thresholds. Zero denominators produce `NA`; requested
filters reject unavailable or non-finite metrics. Fraction metrics use
0–1, not percentages.

## Details

`total.fragments` counts fragment-file records (`frequency_count`), not
column-five read support (`reads_count`). Each file's explicit Signac
barcode mapping is applied separately, so captures may reuse the same
raw barcode. Resources without selected-assay cells are ignored during
calculation but remain attached to the returned assay. `blacklist_ratio`
is the fraction of peak-matrix counts assigned to peaks overlapping any
blacklist interval. Each peak contributes once; it is not a
fragment-level blacklist fraction. Existing QC columns are imported
values; callers must verify their source, units and relevance after
changing peaks. Explicit FRiP source columns must not use names
overwritten by fragment counting (`frequency_count`, `reads_count`,
`mononucleosomal`, `nucleosome_free`, `total.fragments`); retain
upstream values under distinct names when attaching fragments. TSS
enrichment uses Signac's R implementation and does not select an
external backend merely because a newer Signac release exports it.
Missing optional inputs leave metrics unavailable. This function does
not perform doublet detection or establish that all fragment QC is
complete.

## See also

[`ReadATACData()`](https://mengxu98.github.io/scop/reference/ReadATACData.md),
[`CheckATACData()`](https://mengxu98.github.io/scop/reference/CheckATACData.md)

## Examples

``` r
counts <- Matrix::sparseMatrix(
  i = c(1, 2, 1, 2), j = c(1, 1, 2, 2), x = c(2, 1, 1, 3),
  dimnames = list(c("chr1-1-10", "chr1-21-30"), c("cell1", "cell2"))
)
ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
  seqnames = "chr1", start = c(1L, 21L), end = c(10L, 30L)
))
object <- ReadATACData(counts, ranges = ranges, genome = "toy", verbose = FALSE)
object$source_in_peaks <- c(3, 4)
object$source_total <- c(5, 8)
object <- RunATACQC(object,
  frip.columns = c(in_peaks = "source_in_peaks", total = "source_total"),
  min_pct_reads_in_peaks = 0.5, verbose = FALSE
)
object$pct_reads_in_peaks
#> cell1 cell2 
#>   0.6   0.5 
```
