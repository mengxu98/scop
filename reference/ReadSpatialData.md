# Load and validate Visium, Visium HD, or Xenium data

Delegate file parsing to Seurat, validate each image against its assay,
and record observation type, resolution and coordinate units. This
function does not download data, normalize counts, segment cells, or
install packages.

## Usage

``` r
ReadSpatialData(
  data.dir,
  technology = c("visium", "visium_hd", "xenium"),
  bin.size = c(8L, 16L),
  sample_id = "sample1",
  ...
)
```

## Arguments

- data.dir:

  Existing vendor output directory.

- technology:

  One of visium, visium_hd, or xenium.

- bin.size:

  Positive integer HD bin sizes in microns. Ignored for other
  technologies; Seurat stores each resolution in a separate assay/image.

- sample_id:

  Nonempty sample identifier stored in spatial_sample.

- ...:

  Extra arguments to Seurat's Load10X_Spatial or LoadXenium. Xenium
  defaults to cell segmentation and does not load molecule coordinates
  unless requested; cell segmentation must already exist in the vendor
  output.

## Value

A validated Seurat object with misc\$scop_spatial_input provenance.

## See also

SpatialCellPlot
