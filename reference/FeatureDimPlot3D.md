# 3D-Dimensional reduction plot for gene expression visualization.

3D-Dimensional reduction plot for gene expression visualization.

## Usage

``` r
FeatureDimPlot3D(
  srt,
  features,
  reduction = NULL,
  dims = c(1, 2, 3),
  axis_labs = NULL,
  split.by = NULL,
  layer = "data",
  assay = NULL,
  calculate_coexp = FALSE,
  pt.size = 1.5,
  cells.highlight = NULL,
  cols.highlight = "black",
  shape.highlight = "circle-open",
  sizes.highlight = 2,
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

- features:

  A character vector or a named list of features to plot. Features can
  be gene names in Assay or names of numeric columns in meta.data.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a two-length numeric vector specifying x-
  and y-dimensions

- axis_labs:

  A character vector of length 3 indicating the labels for the axes.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- calculate_coexp:

  Whether to calculate the co-expression value (geometric mean) of the
  features.

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

- width:

  Width in pixels, defaults to automatic sizing.

- height:

  Height in pixels, defaults to automatic sizing.

- save:

  The name of the file to save the plot to. Must end in ".html".

- force:

  Whether to force drawing regardless of the number of features greater
  than 100. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[CellDimPlot3D](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
FeatureDimPlot3D(
  pancreas_sub,
  features = c("Ghrl", "Ins1", "Gcg", "Ins2"),
  reduction = "StandardpcaUMAP3D"
)
```
