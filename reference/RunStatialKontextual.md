# Run Statial Kontextual spatial relationships

Compute spatial relationships between cell or spot labels relative to a
parent context using Statial Kontextual.

## Usage

``` r
RunStatialKontextual(
  object,
  group.by,
  r,
  from = NULL,
  to = NULL,
  parent = NULL,
  parent_df = NULL,
  image = NULL,
  sample.by = NULL,
  images = NULL,
  coord.cols = c("col", "row"),
  inhom = FALSE,
  edge_correct = TRUE,
  window = c("convex", "square", "concave"),
  window.length = NA_real_,
  include_original = TRUE,
  cores = 1,
  tool_name = "StatialKontextual",
  store_results = TRUE,
  store_input = FALSE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- group.by:

  Metadata column containing cell or spot labels.

- r:

  Numeric radius or radii used by `Statial::Kontextual()`, expressed in
  the selected coordinate units.

- from, to, parent:

  Cell or spot labels passed to `Statial::Kontextual()`. Ignored when
  `parent_df` is supplied.

- parent_df:

  Optional data frame from `Statial::parentCombinations()`.

- image:

  Name of the Seurat spatial image. Required when multiple images are
  present; a single image is selected automatically when `NULL`.

- sample.by:

  Optional metadata column used as Statial `imageID`. If `NULL`, all
  cells or spots are treated as one image.

- images:

  Optional Statial image filter passed to `Kontextual(image = )`.

- coord.cols:

  Metadata coordinate columns used when no Seurat image coordinates are
  available.

- inhom:

  Whether Statial should account for inhomogeneity.

- edge_correct:

  Whether Statial should perform edge correction.

- window, window.length:

  Window arguments passed to `Statial::Kontextual()`. Numeric window
  lengths use the selected coordinate units.

- include_original:

  Whether to include original L-function values.

- cores:

  Number of cores passed to `Statial::Kontextual()`.

- tool_name:

  Name used to store results in `srt@tools`.

- store_results:

  Whether to store results in `srt@tools`.

- store_input:

  Whether to store the backend input cell table in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- coordinate_space:

  Coordinate system used for spatial relationships. The default is raw
  acquisition coordinates; `"legacy_display"` remains an explicit
  compatibility option.

- ...:

  Additional named arguments passed to `Statial::Kontextual()`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with Statial results stored in
`srt@tools[[tool_name]]` when `store_results = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
# Compare tissue-labelled spots at two radii (full-resolution image pixels).
spatial <- RunStatialKontextual(
  visium_human_pancreas_sub,
  group.by = "coda_label", image = "slice1",
  r = c(500, 1000), from = "collagen", to = "acini", parent = c("collagen", "acini"),
  verbose = FALSE
)
StatialKontextualPlot(spatial)
} # }
```
