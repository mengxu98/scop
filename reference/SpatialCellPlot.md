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
  palette = NULL,
  palcolor = NULL,
  fill.alpha = 0.7,
  boundary.color = "grey30",
  boundary.linewidth = 0.1,
  theme_use = "theme_spatial",
  theme_args = list(),
  assay = NULL,
  layer = "data",
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  boundary = NULL,
  xlim = NULL,
  ylim = NULL,
  image.scale = c("lowres", "hires"),
  overlay_image = TRUE,
  image.alpha = 1,
  flip.y = NULL,
  pt.size = 1,
  molecules = NULL,
  mols.size = 0.1,
  mols.alpha = 1,
  mols.colors = NULL,
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

  Palette name or explicit colors. The default uses \`"Paired"\` for
  categories and sequential \`"YlGnBu"\` for numeric values. Named
  colors are matched to category names.

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

- combine:

  Return combined panels, or a named list when \`FALSE\`.

- nrow, ncol, byrow:

  Layout controls for multiple feature panels.

- legend.position, legend.direction, legend.title:

  Legend controls.

- boundary:

  Boundary name in the image. \`NULL\` selects segmentation when
  available and otherwise plots centroids. Explicit names must exist.

- xlim, ylim:

  Optional increasing ROI limits in raw coordinate units. Cells are
  selected by centroids; boundary vertices are clipped at the limits.

- image.scale, overlay_image, image.alpha:

  Raster scale, visibility, and opacity.

- flip.y:

  Reverse the y axis. \`NULL\` uses the official FOV orientation; use
  \`FALSE\` for the historical FOV layout. Visium rasters include the
  flip.

- pt.size:

  Centroid point size.

- molecules:

  Optional names of already loaded molecules to overlay. No files are
  read. Molecules are filtered by ROI, not assigned to cells.

- mols.size, mols.alpha, mols.colors:

  Molecule size, opacity, and optional named colors matched to molecule
  names.

- ...:

  Arguments passed to other functions.

## Value

A \`ggplot\`, patchwork, or named list of plots.

## Details

Boundary tables require \`cell_id\`, \`x\`, and \`y\` columns in polygon
vertex order, with optional \`polygon_id\` and \`ring_id\` columns for
multipart cells and holes. Object-derived boundaries use the same
display transform as \[SpatialSpotPlot()\]. ROI filtering does not
subset the input object.

## Examples

``` r
# Constructed polygons demonstrate plotting, not a segmentation algorithm.
boundaries <- data.frame(
  cell_id = rep(c("cell1", "cell2"), each = 4),
  x = c(0, 2, 2, 0, 3, 5, 5, 3), y = rep(c(0, 0, 1, 1), 2),
  cell_type = rep(c("A", "B"), each = 4)
)
SpatialCellPlot(boundaries = boundaries, group.by = "cell_type")
```
