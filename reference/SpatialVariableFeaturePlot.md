# Plot spatial variable feature results

Visualize normalized results produced by
[`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md).
The summary view shows feature ranks and, when finite p- or q-values are
stored, significance. Results without permutation statistics are shown
without a significance size mapping or legend. The surface view reuses
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md)
to draw spatial expression for selected features.

## Usage

``` r
SpatialVariableFeaturePlot(
  srt = NULL,
  plot_type = c("summary", "surface", "combined"),
  features = NULL,
  nfeatures = 10,
  score_col = "score",
  assay = NULL,
  layer = NULL,
  image = NULL,
  overlay_image = TRUE,
  image.alpha = 1,
  coord.cols = c("col", "row"),
  flip.y = TRUE,
  pt.size = NULL,
  pt.alpha = 0.9,
  stroke = 0.1,
  palette = "Spectral",
  palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  object = NULL,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A `Seurat` object. The same object may be supplied as `object =` for
  consistency with spatial plotting APIs.

- object:

  Optional alias for `srt`. Supply exactly one of `srt` or `object`.

- plot_type:

  Plot type: `"summary"`, `"surface"`, or `"combined"`.

- features:

  Features to plot. If `NULL`, top features from the stored spatial
  variable feature result are used.

- nfeatures:

  Number of top features used when `features = NULL`.

- score_col:

  Result column used for the summary x-axis.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- overlay_image, image.alpha:

  Draw the spatial image beneath spots.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- flip.y:

  Reverse the y axis for metadata coordinates.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- stroke:

  Point border width.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- legend.position:

  Legend position for surface plots.

- theme_use:

  Theme name or function.

- theme_args:

  Theme name or function, plus extra theme arguments.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- image.scale:

  Image scale factor matching the raster stored in the selected image.
  Use `"hires"` for a hires raster; do not modify Seurat scale-factor
  slots.

## Value

A `ggplot` or `patchwork` object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
spatial <- RunSpatialVariableFeatures(
  spatial,
  assay = "Spatial",
  nfeatures = 10,
  verbose = FALSE
)
SpatialVariableFeaturePlot(spatial, plot_type = "summary")
```
