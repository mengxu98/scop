# Run lightweight SpatialDWLS-style deconvolution

Estimate spot-level cell-type proportions by fitting spatial expression
to reference cell-type signatures with non-negative least squares-style
coefficients. The output follows SCOP's deconvolution metadata
convention so it can be plotted directly with
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

## Usage

``` r
RunSpatialDWLS(
  srt,
  reference,
  reference_label,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  min_cells = 2,
  prefix = "SpatialDWLS",
  tool_name = "SpatialDWLS",
  normalize = TRUE,
  store_results = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  Spatial `Seurat` object used as the RCTD query.

- reference:

  Reference `Seurat` object containing annotated single cells.

- reference_label:

  Metadata column in `reference` with cell type labels.

- assay:

  Assay used in `srt`. If `NULL`, the default assay is used.

- reference_assay:

  Assay used in `reference`.

- layer, reference_layer:

  Assay layers used for spatial and reference raw counts.

- features:

  Features used for RCTD. If `NULL`, shared features are used.

- image:

  Name of the Seurat spatial image used to recover coordinates when
  `coord.cols` are not available.

- coord.cols:

  Metadata coordinate columns used when no image coordinate source is
  requested or available.

- min_cells:

  Minimum reference cells required per cell type.

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- normalize:

  Whether to library-size normalize and `log1p` transform spatial and
  reference matrices before fitting.

- store_results:

  Whether to store detailed RCTD results in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with `"<prefix>_prop_*"`, `"<prefix>_dominant_type"`,
and `"<prefix>_max_prop"` metadata columns. Detailed results are stored
in `srt@tools[[tool_name]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
spatial <- RunSpatialDWLS(
  spatial,
  reference = reference,
  reference_label = "celltype",
  coord.cols = c("x", "y")
)
SpatialSpotPlot(spatial, group.by = "SpatialDWLS_dominant_type")
SpatialSpotPlot(spatial, group.by = "SpatialDWLS_dominant_type", plot_type = "pie")
} # }
```
