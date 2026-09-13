# Plot spatial cell boundaries

Plot cell segmentation polygons from a boundary table or Seurat spatial
image.

## Usage

``` r
SpatialCellPlot(
  object = NULL,
  res = NULL,
  boundaries = NULL,
  cells = NULL,
  image = NULL,
  crop = TRUE,
  group.by = NULL,
  features = NULL,
  palette = "Paired",
  palcolor = NULL,
  fill.alpha = 0.7,
  boundary.color = "grey30",
  boundary.linewidth = 0.1,
  theme_use = "theme_spatial",
  theme_args = list(),
  assay = NULL,
  layer = "data",
  ...
)
```

## Arguments

- object:

  Optional \`Seurat\` object used to extract boundaries or values.

- res:

  Optional result list containing a \`boundaries\` data frame.

- boundaries:

  Optional boundary data frame.

- cells:

  Optional cell or spot identifiers to retain.

- image:

  Seurat image name. Multi-image objects require an explicit name.

- crop:

  Whether to crop the plot to the selected boundaries.

- group.by:

  Boundary column or Seurat metadata column used for filling.

- features:

  Features to display. Multiple features return a patchwork.

- palette, palcolor:

  Palette name or explicit colors.

- fill.alpha:

  Polygon fill opacity.

- boundary.color, boundary.linewidth:

  Boundary appearance.

- theme_use, theme_args:

  Theme used to style the plot. Default is \`"theme_spatial"\`.

- assay:

  Assay used when \`features\` are fetched from \`object\`.

- layer:

  Layer used when \`features\` are fetched from \`object\`.

- ...:

  Additional arguments passed to \`ggplot2::geom_polygon()\`.

## Value

A \`ggplot\` or patchwork object.

## Details

Boundary tables require \`cell_id\`, \`x\`, and \`y\` columns in polygon
vertex order. Seurat images must contain segmentation boundaries.

## Examples

``` r
data(visium_human_pancreas_sub)
SpatialCellPlot(visium_human_pancreas_sub, group.by = "CellType")
#> Error: The selected image does not contain segmentation boundaries
```
