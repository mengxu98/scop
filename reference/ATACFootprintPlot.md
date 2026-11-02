# Plot bias-corrected footprints

Plot bias-corrected footprints.

## Usage

``` r
ATACFootprintPlot(object, motifs, assay = NULL, group.by = NULL,
  normalization = "subtract", show.expected = TRUE, ...)
```

## Arguments

- object:

  Seurat object returned by RunATACFootprint.

- motifs:

  Unique stored footprint IDs.

- assay:

  Chromatin assay.

- group.by:

  Complete cell metadata column, or identities when NULL.

- normalization:

  Native subtract, divide, or NULL normalization.

- show.expected:

  Show expected insertion profile.

- ...:

  Additional Signac::PlotFootprint arguments.

## Value

Native ggplot or patchwork object.

## Details

Missing or stale cell/peak receipts are rejected. This plot does not
establish TF occupancy.
