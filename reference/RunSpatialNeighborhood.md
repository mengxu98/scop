# Run spatial neighborhood statistics

Build a scop-style spatial neighborhood result bundle and optionally
dispatch to a supported backend for colocalization or local-effect
statistics.

## Usage

``` r
RunSpatialNeighborhood(
  srt,
  group.by,
  method = c("spicyR", "mistyR", "Statial", "HoodscanR"),
  assay = NULL,
  layer = "data",
  coord.cols = c("col", "row"),
  image = NULL,
  sample.by = NULL,
  split.by = NULL,
  subject.by = NULL,
  radius = NULL,
  k = NULL,
  features = NULL,
  from = NULL,
  to = NULL,
  tool_name = "SpatialNeighborhood",
  store_results = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- group.by:

  Metadata column containing spatial cell or spot labels.

- method:

  Neighborhood backend. Currently `"spicyR"` is supported.

- assay:

  Assay used when `features` are requested.

- layer:

  Assay layer used when `features` are requested.

- coord.cols:

  Metadata coordinate columns used when no Seurat image coordinates are
  available.

- image:

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- sample.by:

  Metadata column identifying images or samples. If `NULL`, all spots
  are treated as one sample.

- split.by:

  Optional metadata column identifying conditions for differential
  neighborhood statistics.

- subject.by:

  Optional metadata column identifying subjects for backends that
  support paired or repeated designs.

- radius:

  Optional spatial radius used for scop-native neighborhood summaries.

- k:

  Optional number of nearest neighbors used for scop-native neighborhood
  summaries. When both `radius` and `k` are `NULL`, `k = 6` is used.

- features:

  Optional features to extract into the backend input table.

- from, to:

  Optional cell or spot label filters.

- tool_name:

  Name used to store results in `srt@tools`.

- store_results:

  Whether to store results in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to the selected backend.

## Value

A `Seurat` object with results stored in `srt@tools[[tool_name]]`.
