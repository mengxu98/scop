# 3D-Dimensional reduction plot for cell classification visualization.

Plotting cell points on a reduced 3D space and coloring according to the
groups of the cells.

## Usage

``` r
CellDimPlot3D(
  srt,
  group.by,
  reduction = NULL,
  dims = c(1, 2, 3),
  axis_labs = NULL,
  palette = "Chinese",
  palcolor = NULL,
  bg_color = "grey80",
  pt.size = 1.5,
  cells.highlight = NULL,
  cols.highlight = "black",
  shape.highlight = "circle-open",
  sizes.highlight = 2,
  lineages = NULL,
  lineages_palette = "Dark2",
  span = 0.75,
  width = NULL,
  height = NULL,
  save = NULL,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a three-length numeric vector specifying
  x-, y- and z-dimensions

- axis_labs:

  A character vector of length 3 indicating the labels for the axes.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- bg_color:

  Color value for background(NA) points.

- pt.size:

  The size of the points in the plot.

- cells.highlight:

  A logical or character vector specifying the cells to highlight in the
  plot. If `TRUE`, all cells are highlighted. If `FALSE`, no cells are
  highlighted. Default is `NULL`.

- cols.highlight:

  Color used to highlight the cells.

- shape.highlight:

  Shape of the cell to highlight. See
  [scattergl-marker-symbol](https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol)

- sizes.highlight:

  Size of highlighted cell points.

- lineages:

  Lineages/pseudotime to add to the plot. If specified, curves will be
  fitted using [stats::loess](https://rdrr.io/r/stats/loess.html)
  method.

- lineages_palette:

  Color palette used for lineages.

- span:

  The span of the loess smoother for lineages line.

- width:

  Width in pixels, defaults to automatic sizing.

- height:

  Height in pixels, defaults to automatic sizing.

- save:

  The name of the file to save the plot to. Must end in ".html".

- force:

  Whether to force drawing regardless of maximum levels in any cell
  group is greater than 100. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[FeatureDimPlot3D](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(
  pancreas_sub,
  nonlinear_reduction_dims = 3,
)
CellDimPlot3D(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D"
)

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D",
  show_plot = FALSE
)
CellDimPlot3D(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "StandardpcaUMAP3D",
  lineages = "Lineage1"
)
```
