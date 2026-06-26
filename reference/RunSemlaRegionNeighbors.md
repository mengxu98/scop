# Run semla region neighbor detection

Use
[`semla::RegionNeighbors()`](https://spatial-research.github.io/semla/reference/region-neighbors.html)
to identify neighboring spots for selected metadata labels and write the
returned columns to Seurat metadata.

## Usage

``` r
RunSemlaRegionNeighbors(
  srt,
  column_name,
  column_labels = NULL,
  mode = "outer",
  column_key = NULL,
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object with spatial image data.

- column_name:

  Metadata column containing labels.

- column_labels:

  Labels to find neighbors for. If `NULL`, semla uses all labels in
  `column_name`.

- mode:

  Neighbor selection mode passed to semla.

- column_key:

  Prefix for metadata columns returned by semla.

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

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
spatial$region <- ifelse(
  spatial$col > stats::median(spatial$col),
  "right",
  "left"
)

spatial <- RunSemlaRegionNeighbors(
  spatial,
  column_name = "region",
  column_labels = "right",
  mode = "outer",
  column_key = "right_border"
)

grep("right_border", colnames(spatial[[]]), value = TRUE)
SpatialSpotPlot(spatial, group.by = "region")
} # }
```
