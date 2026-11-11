# Plot bias-corrected ATAC footprints

Plot bias-corrected ATAC footprints

## Usage

``` r
ATACFootprintPlot(
  object,
  motifs,
  assay = NULL,
  group.by = NULL,
  normalization = "subtract",
  show.expected = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object returned by RunATACFootprint.

- motifs:

  Unique footprint keys.

- assay:

  Chromatin assay.

- group.by:

  Metadata column, or NULL for identities.

- normalization:

  Signac normalization: \`subtract\`, \`divide\`, or NULL.

- show.expected:

  Show expected insertion profile.

- ...:

  Arguments passed to other functions.

## Value

A ggplot/patchwork object from Signac.
