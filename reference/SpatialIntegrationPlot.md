# Plot spatial integration results

Visualize standardized results produced by
[`RunSpatialIntegration()`](https://mengxu98.github.io/scop/reference/RunSpatialIntegration.md).

## Usage

``` r
SpatialIntegrationPlot(
  srt,
  method = NULL,
  plot_type = c("spatial", "embedding", "alignment", "composition"),
  group.by = NULL,
  sample.by = NULL,
  reduction = NULL,
  cluster_colname = NULL,
  coord.cols = c("col", "row"),
  use_aligned = FALSE,
  tool_name = "SpatialIntegration",
  combine = TRUE,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
)
```

## Arguments

- srt:

  A `Seurat` object containing spatial integration results.

- method:

  Stored integration method. If `NULL`, the active method stored in
  `srt@tools[[tool_name]]` is used.

- plot_type:

  Plot type: `"spatial"`, `"embedding"`, `"alignment"`, or
  `"composition"`.

- group.by:

  Metadata column used for coloring. Defaults to the stored spatial
  domain column.

- sample.by:

  Metadata column used for facets or composition grouping. Defaults to
  the stored sample column.

- reduction:

  Reduction used for embedding plots. Defaults to the stored integration
  reduction.

- cluster_colname:

  Backward-compatible alias for `group.by`.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- use_aligned:

  Whether spatial plots should use aligned coordinates when available.

- tool_name:

  Name of the `srt@tools` entry created by
  [`RunSpatialIntegration()`](https://mengxu98.github.io/scop/reference/RunSpatialIntegration.md).

- combine:

  Whether to combine plots when delegated plotting returns a list.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md)
  or
  [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md).

## Value

A `ggplot`, patchwork object, or list of plots.
