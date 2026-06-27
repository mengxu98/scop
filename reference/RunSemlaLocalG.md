# Run semla local G spatial autocorrelation

Use
[`semla::RunLocalG()`](https://spatial-research.github.io/semla/reference/local-G.html)
on a Staffli-enabled Seurat object. Results are written by semla to
metadata or to an assay, depending on `store_in_metadata`.

## Usage

``` r
RunSemlaLocalG(
  srt,
  features,
  alternative = NULL,
  store_in_metadata = TRUE,
  assay_name = "GiScores",
  image_type = "tissue_lowres",
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object with spatial image data.

- features:

  Features passed to
  [`semla::RunLocalG()`](https://spatial-research.github.io/semla/reference/local-G.html).

- alternative:

  Alternative hypothesis passed to semla. Use `NULL` to keep semla's
  default behavior.

- store_in_metadata:

  Whether semla should store results in metadata.

- assay_name:

  Assay name used by semla when `store_in_metadata = FALSE`.

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
spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
features <- rownames(spatial)[1:3]
spatial[[paste0(features[1], "_localG")]] <- as.numeric(scale(spatial$x))

SpatialSpotPlot(
  spatial,
  group.by = paste0(features[1], "_localG"),
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


if (
  requireNamespace("semla", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
spatial <- RunSemlaLocalG(
  spatial,
  features = features,
  store_in_metadata = TRUE,
  verbose = FALSE
)
}
```
