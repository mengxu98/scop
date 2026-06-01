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

## Value

A `Seurat` object with RCTD proportion columns in metadata and dominant
cell type summaries. When `store_results = TRUE`, detailed results are
also stored in `srt@tools[["RCTD"]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scop)

data(visium_human_pancreas_sub)
data(panc8_sub)

drop_celltypes <- c(
  "epsilon",
  "macrophage",
  "mast",
  "quiescent-stellate",
  "schwann"
)

panc8_rctd <- subset(
  panc8_sub,
  subset = !celltype %in% drop_celltypes
)

table(panc8_sub$celltype)
table(panc8_rctd$celltype)

spatial <- RunRCTD(
  srt = visium_human_pancreas_sub,
  reference = panc8_rctd,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  rctd_mode = "full",
  max_cores = 1,
  prefix = "RCTD"
)

rctd_cols <- grep("^RCTD_prop_", colnames(spatial@meta.data), value = TRUE)
head(spatial@meta.data[, c("RCTD_dominant_type", "RCTD_max_prop", rctd_cols[1:3])])

SpatialSpotPlot(
  spatial,
  group.by = "RCTD_dominant_type",
  theme_use = "theme_scop"
)

SpatialSpotPlot(
  spatial,
  group.by = rctd_cols[1:min(4, length(rctd_cols))],
  palette = "Spectral",
  theme_use = "theme_scop"
)

SpatialSpotPlot(
  spatial,
  group.by = "RCTD_dominant_type",
  plot_type = "pie",
  pie.radius.scale = 0.45,
  theme_use = "theme_scop"
)
} # }
```
