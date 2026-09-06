# Run RCTD spatial deconvolution

Estimate spot-level cell type proportions from a spatial `Seurat` object
using a single-cell `Seurat` reference and the optional `spacexr` RCTD
backend. The example is a non-executing template so pkgdown does not
need to install or run this heavy optional dependency.

## Usage

``` r
RunRCTD(
  srt,
  reference,
  reference_label = "celltype",
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  rctd_mode = c("full", "multi", "doublet"),
  max_cores = 1,
  min_cells = 25,
  prefix = "RCTD",
  store_results = TRUE,
  round_counts = TRUE,
  create_rctd_params = list(),
  run_rctd_params = list(),
  verbose = TRUE,
  ...,
  coordinate_space = c("raw", "legacy_display"),
  tool_name = "RCTD"
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

- rctd_mode:

  RCTD mode passed to `spacexr`. `"full"` is the default for Visium spot
  deconvolution.

- max_cores:

  Number of cores passed to `spacexr`.

- min_cells:

  Minimum number of reference cells required for each cell type. Old
  `spacexr` RCTD requires at least 25 cells per type.

- prefix:

  Prefix for metadata columns.

- store_results:

  Whether to store detailed RCTD results in `srt@tools`.

- round_counts:

  Whether to round non-integer counts to the nearest integer before
  passing data to `spacexr`. RCTD requires integer count matrices; this
  defaults to `TRUE` so bundled example data with scaled non-integer
  reference counts can run directly.

- create_rctd_params:

  Additional parameters passed to
  [`spacexr::createRctd()`](https://rdrr.io/pkg/spacexr/man/createRCTD.html)
  or `spacexr::create.RCTD()`.

- run_rctd_params:

  Additional parameters passed to
  [`spacexr::runRctd()`](https://rdrr.io/pkg/spacexr/man/runRCTD.html)
  or `spacexr::run.RCTD()`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to the RCTD run step.

- coordinate_space:

  Coordinate space used for distance-sensitive input. The default is raw
  acquisition coordinates, so backend distances use raw coordinate
  units. Use `"legacy_display"` explicitly to reproduce the
  display-scaled coordinates used before scop 0.9.0.

- tool_name:

  Name used to store the plain result bundle in `srt@tools`.

## Value

A `Seurat` object with RCTD proportion columns in metadata and dominant
cell type summaries. When `store_results = TRUE`, detailed results are
also stored in `srt@tools[[tool_name]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
thisutils::check_r("dmcable/spacexr", verbose = FALSE)
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
reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
features_use <- head(intersect(
  SeuratObject::VariableFeatures(reference),
  rownames(spatial)
), 300)
spatial <- RunRCTD(
  srt = spatial,
  reference = reference,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  features = features_use,
  coord.cols = c("x", "y"),
  rctd_mode = "full",
  max_cores = 1,
  min_cells = 25,
  prefix = "RCTD",
  verbose = FALSE
)
SpatialDeconvolutionPlot(
  spatial,
  tool_name = "RCTD",
  plot_type = "dominant",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
} # }
```
