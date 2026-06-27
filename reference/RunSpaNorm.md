# Run SpaNorm spatial normalization

Normalize spatial transcriptomics counts with the optional Bioconductor
`SpaNorm` backend and store the normalized expression in a new Seurat
assay.

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
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used for expression values.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

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
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:80],
  features = rownames(visium_human_pancreas_sub)[1:300]
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

SpatialSpotPlot(
  spatial,
  features = rownames(spatial)[1:2],
  assay = "Spatial",
  layer = "data",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


if (
  requireNamespace("SpaNorm", quietly = TRUE) &&
    requireNamespace("SpatialExperiment", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
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
    features = rownames(spatial[["SpaNorm"]])[1:2],
    assay = "SpaNorm",
    layer = "data",
    overlay_image = FALSE,
    coord.cols = c("x", "y")
  )
}
```
