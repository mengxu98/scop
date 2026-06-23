# Run semla radial distance analysis

Use
[`semla::RadialDistance()`](https://spatial-research.github.io/semla/reference/radial-distance.html)
to calculate distances from selected spatial regions and write the
returned columns to Seurat metadata.

## Usage

``` r
RunSemlaRadialDistance(
  srt,
  column_name,
  selected_groups = NULL,
  column_suffix = NULL,
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object with spatial image data.

- column_name:

  Metadata column containing region labels.

- selected_groups:

  Region labels used by semla. If `NULL`, semla uses all labels in
  `column_name`.

- column_suffix:

  Optional suffix for metadata columns returned by semla.

- image_type:

  Image scale used by
  [`semla::UpdateSeuratForSemla()`](https://spatial-research.github.io/semla/reference/UpdateSeuratForSemla.html)
  when the object does not already contain a Staffli object.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to semla.

## Value

A `Seurat` object.
