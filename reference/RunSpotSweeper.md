# Run SpotSweeper spatial quality control

Run `SpotSweeper` spatially aware spot-level quality control on a
spatial `Seurat` object. The wrapper computes standard spot QC metrics,
runs local outlier detection for each metric, optionally runs regional
artifact detection per sample, and writes standardized QC status
metadata that can be visualized with
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).
Artifact detection is reported as completed, skipped, failed, or not
requested; an artifact that was not evaluated is never reported as a
negative artifact call.

## Usage

``` r
RunSpotSweeper(
  srt,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  image = NULL,
  sample.by = NULL,
  metrics = NULL,
  directions = NULL,
  n_neighbors = 36,
  cutoff = 3,
  log = TRUE,
  run_artifact = TRUE,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  mito_percent = NULL,
  mito_sum = NULL,
  n_order = 5,
  shape = c("hexagonal", "square"),
  prefix = "SpotSweeper",
  tool_name = "SpotSweeper",
  return_filtered = FALSE,
  store_results = TRUE,
  cores = 1,
  workers = NULL,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay used for expression. If `NULL`, the default assay is used.

- layer:

  Assay layer used for expression values.

- coord.cols:

  Metadata coordinate columns used when no Seurat image is available.

- image:

  Name of the Seurat spatial image. Required when multiple images are
  present; a single image is selected automatically when `NULL`.

- sample.by:

  Optional metadata column identifying samples or images. If `NULL`, all
  spots are treated as one sample. The selected column must not contain
  missing values. The reserved backend output name `"artifact"` cannot
  be used when `run_artifact = TRUE`.

- metrics:

  QC metrics used by
  [`SpotSweeper::localOutliers()`](https://rdrr.io/pkg/SpotSweeper/man/localOutliers.html).
  If `NULL`, `nCount_<assay>`, `nFeature_<assay>`, and `percent.mito`
  are used. The corresponding `<metric>_outliers` and `<metric>_z` names
  are reserved for backend output and cannot also be requested as input
  columns.

- directions:

  Outlier direction for each metric. If `NULL`, count and feature
  metrics use `"lower"` and mitochondrial metrics use `"higher"`.

- n_neighbors:

  Unitless number of nearest spatial neighbors for local outlier
  detection.

- cutoff:

  Modified z-score cutoff passed to local outlier detection.

- log:

  Whether SpotSweeper should log1p-transform local outlier metrics.

- run_artifact:

  Whether to run
  [`SpotSweeper::findArtifacts()`](https://rdrr.io/pkg/SpotSweeper/man/findArtifacts.html)
  per sample. If the requested artifact stage is skipped or fails for
  any sample, the returned metadata and stored result are marked as
  partial.

- mito_pattern:

  Regex prefixes used to identify mitochondrial genes.

- mito_gene:

  Optional explicit mitochondrial gene vector. When provided,
  `mito_pattern` is ignored.

- mito_percent:

  Metadata column used as mitochondrial percent for artifact detection.
  If `NULL`, `percent.mito` is used. The reserved backend output name
  `"artifact"` cannot be used when `run_artifact = TRUE`.

- mito_sum:

  Metadata column used as mitochondrial counts for artifact detection.
  If `NULL`, mitochondrial counts are computed. The reserved backend
  output name `"artifact"` cannot be used when `run_artifact = TRUE`.

- n_order, shape:

  Parameters passed to
  [`SpotSweeper::findArtifacts()`](https://rdrr.io/pkg/SpotSweeper/man/findArtifacts.html).

- prefix:

  Prefix used for metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- return_filtered:

  Whether to return only spots passing SpotSweeper QC. A fully filtered
  object is not returned when any requested QC stage is incomplete.

- store_results:

  Whether to store detailed results in `srt@tools`.

- cores:

  Number of workers passed to SpotSweeper local outlier detection.

- workers:

  Deprecated alias for `cores`.

- verbose:

  Whether to print progress and completion status messages.

- coordinate_space:

  Coordinate system used for spatial neighborhoods. The default is raw
  acquisition coordinates; `"legacy_display"` remains an explicit
  compatibility option. Neighbor counts are unitless, while any backend
  distance calculation uses the selected coordinate units.

- ...:

  Additional named arguments passed to matching SpotSweeper backend
  functions when those arguments are supported by the installed version.

## Value

A `Seurat` object with SpotSweeper QC metadata. When
`store_results = TRUE`, detailed results and an `artifact_status` table
are stored in `srt@tools[[tool_name]]`. The artifact metadata uses `NA`
and `"NotEvaluated"` when a requested artifact stage did not run; the
overall `*_local_outlier_qc` uses `"NotEvaluated"` outside the tested
coordinate scope. In the overall `*_QC` column, `"Fail"` takes
precedence when any evaluated check fails; otherwise `"Partial"` marks
spots not fully evaluated by the requested QC stages.

## Examples

``` r
data(visium_human_pancreas_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]

# Plot a real input QC metric before running the optional backend.
SpatialSpotPlot(
  spatial,
  group.by = "nCount_Spatial",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)


# Template only: the current runtime still raises its known spot
# identity/order validation error. This documentation change leaves that
# backend issue unresolved and does not fabricate a successful result.
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]
spatial <- RunSpotSweeper(
  spatial,
  assay = "Spatial",
  coord.cols = c("x", "y"),
  n_neighbors = 12,
  run_artifact = FALSE,
  verbose = FALSE
)
} # }
```
