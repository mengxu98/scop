# Run Statial Kontextual spatial relationships

Run `Statial::Kontextual()` on a spatial `Seurat` object to quantify
pairwise cell or spot label relationships relative to a parent context.
Results are stored as a compact SCOP bundle with raw Statial output,
standardized summary, and parameters. `Statial` is an optional
Bioconductor dependency installable with
`BiocManager::install("Statial")`.

## Usage

``` r
RunStatialKontextual(
  srt,
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
  ...
)
```

## Arguments

- srt:

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

## Value

A `Seurat` object with Statial results stored in
`srt@tools[[tool_name]]` when `store_results = TRUE`.

## Examples

``` r
data(visium_human_pancreas_results_sub)
statial <- visium_human_pancreas_results_sub@tools$StatialKontextual
statial$summary
#> $n_records
#> [1] 1
#> 
#> $n_images
#> [1] 1
#> 
#> $n_tests
#> [1] 1
#> 
#> $radii
#> [1] 50
#> 
#> $n_localized
#> [1] 0
#> 
#> $n_dispersed
#> [1] 1
#> 
#> $top_relationships
#>   imageID            test original kontextual  r inhomL
#> 1 sample1 collagen__acini      -50        -50 50  FALSE
#> 
StatialKontextualPlot(res = statial)
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?


if (FALSE) { # \dontrun{
check_r("sydney-informatics-hub/Statial", verbose = FALSE)
spatial <- RunStatialKontextual(
  visium_human_pancreas_results_sub, group.by = "coda_label", r = 50,
  from = "collagen", to = "acini", parent = c("collagen", "acini"),
  coord.cols = c("x", "y"), verbose = FALSE
)
} # }
```
