# Load and validate Visium, Visium HD, Xenium, or Atera data

Delegate file parsing to Seurat, validate each image against its assay,
and record observation type, resolution and coordinate units. This
function does not download data, normalize counts, segment cells, or
install packages.

## Usage

``` r
ReadSpatialData(
  data.dir,
  technology = c("visium", "visium_hd", "xenium", "atera"),
  bin.size = c(8L, 16L),
  sample_id = "sample1",
  ...
)
```

## Arguments

- data.dir:

  Existing vendor output directory.

- technology:

  One of visium, visium_hd, xenium, or atera.

- bin.size:

  Positive integer HD bin sizes in microns, \`"polygons"\`, or a mixture
  such as \`c(8, "polygons")\`. Ignored for other technologies; Seurat
  stores each resolution in a separate assay/image.

- sample_id:

  Nonempty sample identifier stored in spatial_sample.

- ...:

  Arguments to the official Seurat loader. For HD, use
  \`segmentation.type = "cell"\` or \`"nucleus"\` and \`compact =
  TRUE\`. Atera retains official defaults, including deferred
  segmentations and molecules; use \`bpcells.dir\` for its counts cache.
  \`genes\` selects molecules, not expression-matrix rows. Additional
  Atera data can be loaded later with Seurat's molecule, segmentation,
  and codeword loaders.

## Value

A validated Seurat object with misc\$scop_spatial_input provenance.

## See also

SpatialCellPlot
