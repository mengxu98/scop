# Bias-corrected ATAC footprinting

Bias-corrected ATAC footprinting.

## Usage

``` r
RunATACFootprint(object, genome, genome.id, bias.region, motifs,
  assay = NULL, upstream = 250L, downstream = 250L,
  in.peaks = FALSE, verbose = TRUE)
```

## Arguments

- object:

  Seurat object with a ChromatinAssay.

- genome:

  Loaded BSgenome or indexed FaFile; no reference is downloaded.

- genome.id:

  Verified reference label matching the assay.

- bias.region:

  Explicit chromosome-start-end interval with three-base boundary
  margins.

- motifs:

  Unique motif IDs with stored positions.

- assay:

  Chromatin assay.

- upstream, downstream:

  Integer flanking widths, at least 50 bases.

- in.peaks:

  Restrict motif positions to peaks.

- verbose:

  Print progress.

## Value

Seurat object with native bias and footprints, plus assay Misc
scop_footprints receipts.

## Details

Native Signac bias is recalculated on an explicit region. Existing
footprint keys are not overwritten. Adequate hexamer coverage is
required; native failures are preserved. Footprints are insertion
patterns, not proof of TF binding. Recompute after changing reference,
cells, peaks or fragments.
