# Plot cell2location spatial results

Plot cell2location spatial results

## Usage

``` r
Cell2locationPlot(
  srt,
  plot_type = c("proportion", "abundance", "dominant", "pie"),
  cell_types = NULL,
  prefix = "Cell2location",
  tool_name = "Cell2location",
  image = NULL,
  overlay_image = TRUE,
  coord.cols = c("col", "row"),
  ...,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A `Seurat` object returned by
  [`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md).

- plot_type:

  Result to draw: normalized proportion, q05 absolute abundance,
  dominant cell type, or a spot-level proportion pie.

- cell_types:

  Optional cell types to display for abundance, proportion, or pie
  plots.

- prefix:

  Metadata prefix used by
  [`RunCell2location()`](https://mengxu98.github.io/scop/reference/RunCell2location.md).

- tool_name:

  Name of the `srt@tools` result entry.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- overlay_image:

  Whether to draw the selected spatial image.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

- image.scale:

  Image scale factor matching the raster stored in the selected image.
  Use `"hires"` for a hires raster; do not modify Seurat scale-factor
  slots.

## Value

A `ggplot`, `patchwork`, or list of plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Result from the official Human Lymph Node example in RunCell2location().
# The tutorial result is an external file and is not downloaded by pkgdown.
spatial <- readRDS("human_lymph_node_cell2location/official_human_lymph_node.rds")
selected <- names(sort(
  colMeans(spatial@tools$Cell2location$proportions),
  decreasing = TRUE
))[1:6]
Cell2locationPlot(
  spatial,
  plot_type = "proportion",
  cell_types = selected,
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  ncol = 3
)
Cell2locationPlot(
  spatial,
  plot_type = "dominant",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
} # }
```
