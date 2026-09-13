# Build a native spatial network

Build a k-nearest-neighbor or radius spatial network from raw Seurat
spatial coordinates. Results are stored as named graphs in
\`srt@tools\$SpatialNetwork\`.

## Usage

``` r
RunSpatialNetwork(
  object = NULL,
  method = c("knn", "radius"),
  image = NULL,
  coord.cols = c("col", "row"),
  k = 6,
  radius = NULL,
  graph.name = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A \`Seurat\` object.

- method:

  Network method, either \`"knn"\` or \`"radius"\`.

- image:

  Seurat image name. A single image is selected automatically;
  multi-image objects require an explicit value.

- coord.cols:

  Metadata columns used when the object has no image.

- k:

  Number of neighbors for \`method = "knn"\`.

- radius:

  Positive distance threshold for \`method = "radius"\`, expressed in
  the raw coordinate units.

- graph.name:

  Optional graph name. If \`NULL\`, a deterministic name is generated
  from the image, method, and method parameter.

- overwrite:

  Whether an existing graph with the same name may be replaced.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for \`object\`; supply exactly one of the two. It
  will be removed in scop 1.0.0.

## Value

The input \`Seurat\` object with a \`SpatialNetwork\` result in
\`srt@tools\`.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
spatial <- RunSpatialNetwork(object = spatial, k = 6, verbose = FALSE)
SpatialNetworkPlot(object = spatial, group.by = "coda_label")

```
