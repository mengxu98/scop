# Persist the coordinate axis convention of a VisiumV2 image

Record whether the centroid column named `x` represents the horizontal
image-column axis or the vertical image-row axis. This does not move
coordinates or change the raster. It replaces metadata-based guessing
and survives Seurat subsetting and cell renaming. Stored spatial
analyses must be rerun after changing the convention.

## Usage

``` r
SetSpatialImageAxes(
  object,
  image = NULL,
  x_orientation = c("horizontal", "vertical")
)
```

## Arguments

- object:

  A Seurat object.

- image:

  Image name. Required when several images are present.

- x_orientation:

  `"horizontal"` for custom images whose centroid x is image-column, or
  `"vertical"` for Seurat Read10X_Image's positional imagerow/imagecol
  convention.

## Value

A Seurat object with an explicit image coordinate convention.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- SetSpatialImageAxes(
  visium_human_pancreas_sub, image = "slice1", x_orientation = "horizontal"
)
SpatialCoordinates(spatial, image = "slice1")$source$coord.cols
#> [1] "x" "y"
```
