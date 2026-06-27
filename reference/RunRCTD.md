# Run RCTD spatial deconvolution

Estimate spot-level cell type proportions from a spatial `Seurat` object
using a single-cell `Seurat` reference and `spacexr` RCTD.

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
  ...
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

  Additional parameters passed to `spacexr::createRctd()` or
  [`spacexr::create.RCTD()`](https://rdrr.io/pkg/spacexr/man/create.RCTD.html).

- run_rctd_params:

  Additional parameters passed to `spacexr::runRctd()` or
  [`spacexr::run.RCTD()`](https://rdrr.io/pkg/spacexr/man/run.RCTD.html).

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to the RCTD run step.

## Value

A `Seurat` object with RCTD proportion columns in metadata and dominant
cell type summaries. When `store_results = TRUE`, detailed results are
also stored in `srt@tools[["RCTD"]]`.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:120],
  features = rownames(visium_human_pancreas_sub)[1:400]
)
#> Warning: Not validating Centroids objects
#> Warning: Not validating Centroids objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating Seurat objects
rctd_weights <- data.frame(
  RCTD_prop_Ductal = seq(0.75, 0.15, length.out = ncol(spatial)),
  RCTD_prop_Endocrine = seq(0.15, 0.65, length.out = ncol(spatial)),
  RCTD_prop_Stromal = 0.10,
  row.names = colnames(spatial)
)
rctd_weights <- rctd_weights / rowSums(rctd_weights)
spatial <- Seurat::AddMetaData(spatial, rctd_weights)
spatial$RCTD_dominant_type <- sub(
  "^RCTD_prop_",
  "",
  colnames(rctd_weights)[max.col(rctd_weights)]
)
spatial$RCTD_max_prop <- apply(rctd_weights, 1, max)

SpatialSpotPlot(
  spatial,
  group.by = "RCTD_dominant_type",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)

if (requireNamespace("scatterpie", quietly = TRUE)) {
  SpatialSpotPlot(
    spatial,
    group.by = "RCTD_dominant_type",
    plot_type = "pie",
    overlay_image = FALSE,
    coord.cols = c("x", "y")
  )
}


if (
  requireNamespace("spacexr", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
data(pancreas_sub)
features_use <- head(intersect(rownames(spatial), rownames(pancreas_sub)), 300)

spatial <- RunRCTD(
  srt = spatial,
  reference = pancreas_sub,
  reference_label = "CellType",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  features = features_use,
  rctd_mode = "full",
  max_cores = 1,
  min_cells = 5,
  prefix = "RCTD"
)

rctd_cols <- grep("^RCTD_prop_", colnames(spatial@meta.data), value = TRUE)
SpatialSpotPlot(
  spatial,
  group.by = "RCTD_dominant_type",
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  theme_use = "theme_scop"
)
SpatialSpotPlot(
  spatial,
  group.by = rctd_cols[1:min(3, length(rctd_cols))],
  palette = "Spectral",
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  theme_use = "theme_scop"
)
}
```
