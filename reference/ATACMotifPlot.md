# Plot motif sequence logos from a selected ATAC assay

Plot motif sequence logos from a selected ATAC assay

## Usage

``` r
ATACMotifPlot(object, motifs, assay = NULL, use.names = TRUE, ...)
```

## Arguments

- object:

  A Seurat object with native motifs.

- motifs:

  Unique motif IDs (not ambiguous common names).

- assay:

  Chromatin assay.

- use.names:

  Use native motif names in plot labels.

- ...:

  Arguments passed to Signac::MotifPlot and ggseqlogo.

## Value

Native motif-logo plot. Motif similarity is not TF binding evidence.

## Details

Requires the optional ggseqlogo package.
