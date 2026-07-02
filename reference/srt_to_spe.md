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
  include_meta = TRUE
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

## Value

A `SpatialExperiment`.
