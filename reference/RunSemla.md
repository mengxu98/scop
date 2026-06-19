# Run semla spatial analysis backends

Lightweight wrappers around optional `semla` spatial analysis methods.
`RunSemlaSpatialNetwork()` stores the network in `srt@tools`;
`RunSemlaLocalG()`, `RunSemlaRegionNeighbors()`, and
`RunSemlaRadialDistance()` delegate result placement to semla.

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

  Whether to print progress messages.

- ...:

  Additional arguments passed to semla.

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

- column_name:

  Metadata column containing labels or region labels.

- column_labels:

  Labels to find neighbors for. If `NULL`, semla uses all labels in
  `column_name`.

- mode:

  Neighbor selection mode passed to semla.

- column_key:

  Prefix for metadata columns returned by semla.

- selected_groups:

  Region labels used by semla. If `NULL`, semla uses all labels in
  `column_name`.

- column_suffix:

  Optional suffix for metadata columns returned by semla.

## Value

A `Seurat` object.
