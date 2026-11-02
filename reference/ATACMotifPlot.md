# Plot motif sequence logos

Plot motif sequence logos.

## Usage

``` r
ATACMotifPlot(object, motifs, assay = NULL, use.names = TRUE, ...)
```

## Arguments

- object:

  Seurat object with motifs.

- motifs:

  Unique motif IDs in the selected assay PWM list.

- assay:

  Chromatin assay.

- use.names:

  Use native motif display names.

- ...:

  Arguments passed to Signac::MotifPlot and ggseqlogo.

## Value

Native motif sequence logo plot.

## Details

Uses exact assay-specific IDs rather than ambiguous common names.
Requires ggseqlogo.
