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
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)

spatial <- RunSemlaSpatialNetwork(
  visium_human_pancreas_sub,
  nNeighbors = 6,
  coords = "array"
)

head(spatial@tools$SemlaSpatialNetwork$network)
SpatialSpotPlot(spatial, group.by = "orig.ident")
} # }
```
