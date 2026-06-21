# Run multi-sample spatial integration

Integrate multi-slice or multi-sample spatial transcriptomics data with
an optional spatial backend and store standardized embeddings, domains,
and aligned coordinates in a `Seurat` object.

## Usage

``` r
RunSpatialIntegration(
  object,
  method = c("PRECAST", "BASS", "SpatialMNN"),
  sample.by = NULL,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  features = NULL,
  image = NULL,
  reduction.name = NULL,
  cluster_colname = NULL,
  tool_name = "SpatialIntegration",
  store_results = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A merged spatial `Seurat` object or a list of spatial `Seurat`
  objects.

- method:

  Spatial integration backend.

- sample.by:

  Metadata column identifying samples for a merged `Seurat` object. For
  list input, list names are copied into this column.

- assay:

  Assay used for expression extraction. If `NULL`, the default assay is
  used.

- layer:

  Assay layer used as backend input.

- coord.cols:

  Metadata coordinate columns used when no image coordinate source is
  available.

- features:

  Optional features to use. If `NULL`, shared assay features are used.

- image:

  Name of the Seurat spatial image used to recover coordinates when
  available.

- reduction.name:

  Name of the integrated embedding reduction. If `NULL`, a
  method-specific name is used.

- cluster_colname:

  Metadata column used for spatial domain labels. If `NULL`, a
  method-specific name is used.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store detailed results in `srt@tools`.

- verbose:

  Whether to print progress messages.

- ...:

  Additional backend-specific arguments.

## Value

A `Seurat` object with spatial integration results stored in metadata,
reductions, and `srt@tools[[tool_name]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
srt <- RunSpatialIntegration(
  object = spatial,
  method = "PRECAST",
  sample.by = "sample",
  assay = "Spatial",
  coord.cols = c("col", "row")
)

SpatialIntegrationPlot(srt, plot_type = "spatial")
SpatialIntegrationPlot(srt, plot_type = "embedding")
} # }
```
