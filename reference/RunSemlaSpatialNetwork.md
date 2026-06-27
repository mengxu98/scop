# Run semla spatial network construction

Use the optional `semla` package as a backend to prepare a
Staffli-enabled Seurat object and compute spot-level spatial networks.
The network is stored in `srt@tools[[tool_name]]` when
`store_results = TRUE`.

## Usage

``` r
RunSemlaSpatialNetwork(
  srt,
  image_type = "tissue_lowres",
  nNeighbors = 6,
  maxDist = NULL,
  minK = 0,
  coords = "pixels",
  tool_name = "SemlaSpatialNetwork",
  store_results = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object with spatial image data.

- image_type:

  Image scale used by
  [`semla::UpdateSeuratForSemla()`](https://spatial-research.github.io/semla/reference/UpdateSeuratForSemla.html)
  when the object does not already contain a Staffli object.

- nNeighbors:

  Number of nearest spatial neighbors.

- maxDist:

  Optional maximum neighbor distance.

- minK:

  Minimum number of retained neighbors per spot.

- coords:

  Coordinate system passed to
  [`semla::GetSpatialNetwork()`](https://spatial-research.github.io/semla/reference/get-network.html).

- tool_name:

  Name used to store results in `srt@tools`.

- store_results:

  Whether to store the semla spatial network in `srt@tools`.

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
spatial@tools$SemlaSpatialNetwork <- list(
  network = data.frame(
    from = colnames(spatial)[1:6],
    to = colnames(spatial)[2:7],
    distance = sqrt(diff(spatial$x[1:7])^2 + diff(spatial$y[1:7])^2)
  )
)

head(spatial@tools$SemlaSpatialNetwork$network)
#>                                  from                 to distance
#> ATTATCTCGACAGATC-1 TGGTATCGGTCTGTAT-1 ATTATCTCGACAGATC-1       55
#> TGAGATCAAATACTCA-1 ATTATCTCGACAGATC-1 TGAGATCAAATACTCA-1       55
#> CTGGTCCTAACTTGGC-1 TGAGATCAAATACTCA-1 CTGGTCCTAACTTGGC-1       55
#> ATAGTCTTTGACGTGC-1 CTGGTCCTAACTTGGC-1 ATAGTCTTTGACGTGC-1       55
#> GGGTGGTCCAGCCTGT-1 ATAGTCTTTGACGTGC-1 GGGTGGTCCAGCCTGT-1       55
#> ACACGGCACTATGCAT-1 GGGTGGTCCAGCCTGT-1 ACACGGCACTATGCAT-1       55
SpatialSpotPlot(
  spatial,
  group.by = "coda_label",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


if (
  requireNamespace("semla", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
spatial <- RunSemlaSpatialNetwork(
  spatial,
  nNeighbors = 6,
  coords = "pixels",
  verbose = FALSE
)
}
```
