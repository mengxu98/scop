# Run CARD spatial deconvolution

Estimate spot-level cell-type proportions from a spatial `Seurat` object
using a single-cell `Seurat` reference and the optional `CARD`/`CARDspa`
backend. The example is a non-executing template because CARD
installation and execution are optional.

## Usage

``` r
RunCARD(
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
  sample_varname = NULL,
  minCountGene = 100,
  minCountSpot = 5,
  ct_select = NULL,
  prefix = "CARD",
  tool_name = "CARD",
  store_results = TRUE,
  round_counts = TRUE,
  create_card_params = list(),
  card_deconvolution_params = list(),
  verbose = TRUE,
  ...,
  coordinate_space = c("raw", "legacy_display")
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

- sample_varname:

  Optional metadata column in `reference` containing sample labels. When
  `NULL`, all reference cells are assigned to one sample.

- minCountGene, minCountSpot:

  Filtering parameters passed to `CARD::createCARDObject()` or
  [`CARDspa::createCARDObject()`](https://rdrr.io/pkg/CARDspa/man/createCARDObject.html)
  when supported.

- ct_select:

  Optional cell types to keep in CARD.

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store detailed RCTD results in `srt@tools`.

- round_counts:

  Whether to round non-integer counts to the nearest integer before
  passing data to `spacexr`. RCTD requires integer count matrices; this
  defaults to `TRUE` so bundled example data with scaled non-integer
  reference counts can run directly.

- create_card_params:

  Additional parameters passed to `createCARDObject()`.

- card_deconvolution_params:

  Additional parameters passed to `CARD_deconvolution()`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to the RCTD run step.

- coordinate_space:

  Coordinate space used for distance-sensitive input. The default is raw
  acquisition coordinates, so backend distances use raw coordinate
  units. Use `"legacy_display"` explicitly to reproduce the
  display-scaled coordinates used before scop 0.9.0.

## Value

A `Seurat` object with CARD proportion columns in metadata and dominant
cell type summaries. When `store_results = TRUE`, detailed results are
stored in `srt@tools[[tool_name]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
thisutils::check_r("YingMa0107/CARD", verbose = FALSE)
data(visium_human_pancreas_sub)
data(panc8_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]
reference <- panc8_sub[, panc8_sub@meta.data[["celltype"]] %in%
  c("ductal", "alpha", "beta")]
features_use <- head(intersect(rownames(spatial), rownames(reference)), 300)
spatial <- RunCARD(
  srt = spatial,
  reference = reference,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  features = features_use,
  coord.cols = c("x", "y"),
  verbose = FALSE
)
SpatialDeconvolutionPlot(
  spatial,
  tool_name = "CARD",
  plot_type = "dominant",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
} # }
```
