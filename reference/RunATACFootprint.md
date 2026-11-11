# Bias-corrected ATAC footprinting

Compute Signac insertion bias and footprints with an explicit sequence
reference. No reference is downloaded. Existing cached footprints are
never overwritten.

## Usage

``` r
RunATACFootprint(
  object,
  genome,
  genome.id,
  bias.region,
  motifs,
  assay = NULL,
  upstream = 250L,
  downstream = 250L,
  in.peaks = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing a ChromatinAssay.

- genome:

  A loaded BSgenome or indexed Rsamtools FaFile.

- genome.id:

  Verified reference label matching the assay.

- bias.region:

  One explicit chromosome interval string, for example
  \`chr1-4-1000000\`. It must fit the reference. Bias estimation
  requires adequate coverage of all hexamers and can fail on small
  datasets.

- motifs:

  Unique motif IDs present in the selected assay's positions.

- assay:

  Source assay; defaults to the active assay.

- upstream, downstream:

  Nonnegative flanking widths (at least 50 each).

- in.peaks:

  Restrict motif instances to peaks.

- verbose:

  Print progress.

## Value

A Seurat object with positionEnrichment and a \`scop_footprints\`
receipt in assay Misc. Input counts, cells and default assay are
preserved. At least two assay cells are required. A receipt is a
snapshot; recalculate after changing cells, peaks, fragments or
reference.

## Details

Requires indexed local fragment files with complete cell/barcode
mappings. Bias is recalculated using the explicit region, never silently
reused from another reference. Estimation errors are preserved.
Footprints summarize Tn5 insertion patterns, not proven TF binding or
causality. Factors sharing a motif cannot be distinguished by this
analysis.

## See also

\[ATACFootprintPlot()\]
