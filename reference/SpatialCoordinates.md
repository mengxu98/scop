# Read spatial coordinates with an explicit coordinate contract

Return raw analysis coordinates or display coordinates together with
their source and reversible transform. Image-backed coordinates are
normalized to horizontal `x` and vertical `y`. Seurat's standard
`VisiumV2` loader stores source `imagerow`/`imagecol` values
positionally as centroid `x`/`y`; SCOP restores horizontal image-column
and vertical image-row order. Generic FOV centroid tables with explicit
`x`/`y` retain their order. SCOP then applies only the selected scale
and display-axis flip. Custom VisiumV2 loaders storing horizontal
centroid `x` must call
[`SetSpatialImageAxes()`](https://mengxu98.github.io/scop/reference/SetSpatialImageAxes.md)
with `x_orientation = "horizontal"` once. This persists the convention
in the Seurat misc slot through subsetting and renaming; `"vertical"`
(or no marker) uses the standard Seurat loader convention. Legacy
`coords_x_orientation` image attributes are readable but should be
migrated with
[`SetSpatialImageAxes()`](https://mengxu98.github.io/scop/reference/SetSpatialImageAxes.md)
before subsetting. Cell metadata never determines image axis order.
Metadata-only coordinates retain their numerical orientation in both raw
and display space. This function does not modify the object.

## Usage

``` r
SpatialCoordinates(
  object,
  image = NULL,
  coord.cols = c("col", "row"),
  space = c("raw", "display"),
  image_policy = "strict",
  image.scale = c("lowres", "hires")
)
```

## Arguments

- object:

  A `Seurat` object.

- image:

  Optional Seurat image name.

- coord.cols:

  Metadata columns used when no image is available.

- space:

  Coordinate space to return.

- image_policy:

  Multi-image selection policy. The default requires an explicit image
  when more than one image is available.

- image.scale:

  Image scale factor used only when `space = "display"`. Choose
  `"hires"` when the image slot contains a hires raster. Raw coordinates
  are always returned in full-resolution units.

## Value

A plain list with `data`, `source`, and `transform` entries.

## Examples

``` r
data(visium_human_pancreas_sub)
coords_raw <- SpatialCoordinates(
  visium_human_pancreas_sub,
  image = "slice1",
  coord.cols = c("x", "y"),
  space = "raw"
)
coords_display <- SpatialCoordinates(
  visium_human_pancreas_sub,
  image = "slice1",
  coord.cols = c("x", "y"),
  space = "display",
  image.scale = "lowres"
)
head(coords_raw$data[, c("cell_id", "x", "y")])
#>                               cell_id    x   y
#> TGGTATCGGTCTGTAT-1 TGGTATCGGTCTGTAT-1 3859 662
#> ATTATCTCGACAGATC-1 ATTATCTCGACAGATC-1 3914 662
#> TGAGATCAAATACTCA-1 TGAGATCAAATACTCA-1 3969 662
#> CTGGTCCTAACTTGGC-1 CTGGTCCTAACTTGGC-1 4024 662
#> ATAGTCTTTGACGTGC-1 ATAGTCTTTGACGTGC-1 4079 662
#> GGGTGGTCCAGCCTGT-1 GGGTGGTCCAGCCTGT-1 4134 662
head(coords_display$data[, c("cell_id", "x", "y")])
#>                               cell_id        x        y
#> TGGTATCGGTCTGTAT-1 TGGTATCGGTCTGTAT-1 301.4844 288.2812
#> ATTATCTCGACAGATC-1 ATTATCTCGACAGATC-1 305.7812 288.2812
#> TGAGATCAAATACTCA-1 TGAGATCAAATACTCA-1 310.0781 288.2812
#> CTGGTCCTAACTTGGC-1 CTGGTCCTAACTTGGC-1 314.3750 288.2812
#> ATAGTCTTTGACGTGC-1 ATAGTCTTTGACGTGC-1 318.6719 288.2812
#> GGGTGGTCCAGCCTGT-1 GGGTGGTCCAGCCTGT-1 322.9688 288.2812
```
