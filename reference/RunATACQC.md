# Run scATAC quality control metrics

Calculate common scATAC QC metrics and optionally filter cells by
thresholds.

## Usage

``` r
RunATACQC(
  srt,
  assay = NULL,
  tss.positions = NULL,
  blacklist = NULL,
  fast = TRUE,
  min_pct_reads_in_peaks = NULL,
  min_TSS_enrichment = NULL,
  max_nucleosome_signal = NULL,
  max_blacklist_ratio = NULL,
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

- tss.positions:

  TSS positions passed to
  [`Signac::TSSEnrichment`](https://stuartlab.org/signac/reference/TSSEnrichment.html).

- blacklist:

  A `GRanges` blacklist used to compute `blacklist_ratio`.

- fast:

  Whether to use the fast mode in
  [`Signac::TSSEnrichment`](https://stuartlab.org/signac/reference/TSSEnrichment.html).

- min_pct_reads_in_peaks, min_TSS_enrichment, max_nucleosome_signal,
  max_blacklist_ratio:

  Optional thresholds used for filtering cells.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with QC metadata added.

## Examples

``` r
# \donttest{
data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub <- RunATACQC(
  pbmcmultiome_sub,
  assay = "peaks",
  fast = TRUE
)
#> ℹ [2026-06-26 11:30:09] Calculating ATAC QC metrics...
#> ! [2026-06-26 11:30:09] Skip nucleosome signal: "No fragment files present in assay"
#> ! [2026-06-26 11:30:09] Skip FRiP calculation: no total fragment count column or local fragments available
#> ✔ [2026-06-26 11:30:09] ATAC QC completed
# }
```
