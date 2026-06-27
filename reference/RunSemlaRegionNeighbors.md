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
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:120],
  features = rownames(visium_human_pancreas_sub)[1:400]
)
#> Warning: Not validating Centroids objects
#> Warning: Not validating Centroids objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating Seurat objects
spatial$region <- ifelse(
  spatial$x > stats::median(spatial$x),
  "right",
  "left"
)
spatial$right_border <- spatial$region == "right" &
  abs(spatial$x - stats::median(spatial$x)) < stats::sd(spatial$x) * 0.25

SpatialSpotPlot(
  spatial,
  group.by = c("region", "right_border"),
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


if (
  requireNamespace("semla", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
spatial <- RunSemlaRegionNeighbors(
  spatial,
  column_name = "region",
  column_labels = "right",
  mode = "outer",
  column_key = "right_border",
  verbose = FALSE
)
}
```
