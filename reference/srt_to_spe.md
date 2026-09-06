# Convert Seurat to SpatialExperiment

Create a lightweight `SpatialExperiment` from a spatial Seurat object
using one assay layer, metadata, and resolved spatial coordinates.

## Usage

``` r
srt_to_spe(
  srt,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  image = NULL,
  include_meta = TRUE,
  coordinate_space = c("raw", "legacy_display")
)
```

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay to export. If `NULL`, the default assay is used.

- layer:

  Assay layer to export.

- coord.cols:

  Metadata coordinate columns. By default, SCOP resolves `x/y` first and
  then `col/row`.

- image:

  Optional Seurat image name. When present, image-derived coordinates
  are used.

- include_meta:

  Whether to include Seurat metadata as `colData`.

- coordinate_space:

  Coordinate space exported to `spatialCoords`. `"legacy_display"`
  preserves the historical scaled/y-flipped behavior; `"raw"` (the
  default) preserves analysis distances. The source and transform are
  stored in `metadata(spe)$scop_spatial_coordinates` so that an explicit
  display export can be inverted by
  [`spe_to_srt()`](https://mengxu98.github.io/scop/reference/spe_to_srt.md).

## Value

A `SpatialExperiment`.
