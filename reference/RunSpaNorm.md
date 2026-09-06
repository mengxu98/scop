# Run SpaNorm spatial normalization

Normalize spatial transcriptomics counts with the optional Bioconductor
`SpaNorm` backend and store the normalized expression in a new Seurat
assay. The example is a non-executing template because the optional
backend and its platform-specific numerical requirements are not part of
a standard SCOP installation.

## Usage

``` r
RunSpaNorm(
  srt,
  assay = NULL,
  layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  new_assay = "SpaNorm",
  tool_name = "SpaNorm",
  store_results = TRUE,
  store_spe = FALSE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer used for expression values.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- new_assay:

  Name of the assay used to store SpaNorm-normalized data.

- tool_name:

  Name used to store detailed SpaNorm results in `srt@tools`.

- store_results:

  Whether to store the full result in `srt@tools`.

- store_spe:

  Whether to store the backend `SpatialExperiment` returned by
  `SpaNorm`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- coordinate_space:

  Coordinate system supplied to SpaNorm. The default is raw acquisition
  coordinates; use `"legacy_display"` explicitly to reproduce the
  display-scaled coordinates used before scop 0.9.0.

- ...:

  Additional arguments passed to
  [`SpaNorm::SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.html),
  such as `sample.p`.

## Value

A `Seurat` object with SpaNorm-normalized expression stored in
`new_assay`. When `store_results = TRUE`, parameters, coordinates,
features, cells, and optional backend output are stored in
`srt@tools[[tool_name]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
thisutils::check_r("SpaNorm", verbose = FALSE)
data(visium_human_pancreas_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]
spatial <- RunSpaNorm(
  spatial,
  assay = "Spatial",
  layer = "counts",
  coord.cols = c("x", "y"),
  new_assay = "SpaNorm",
  store_spe = FALSE,
  sample.p = 0.25,
  verbose = FALSE
)
SpatialSpotPlot(
  spatial,
  features = rownames(spatial[["SpaNorm"]])[1],
  assay = "SpaNorm",
  layer = "data",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
} # }
```
