# Convert SpatialExperiment to Seurat

Create a Seurat object from a `SpatialExperiment`, preserving `colData`
and spatial coordinates as metadata columns.

## Usage

``` r
spe_to_srt(
  spe,
  assay = "Spatial",
  layer = NULL,
  coord.cols = c("x", "y"),
  project = "SpatialExperiment"
)
```

## Arguments

- spe:

  A `SpatialExperiment` or `SummarizedExperiment`.

- assay:

  Assay name for the created Seurat assay.

- layer:

  Assay from `spe` to use as counts. If `NULL`, the first assay is used.

- coord.cols:

  Metadata names used for spatial coordinates in Seurat.

- project:

  Project name passed to
  [`Seurat::CreateSeuratObject()`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html).

## Value

A `Seurat` object.
