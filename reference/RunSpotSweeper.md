# Run SpotSweeper spatial quality control

Run `SpotSweeper` spatially aware spot-level quality control on a
spatial `Seurat` object. The wrapper computes standard spot QC metrics,
runs local outlier detection for each metric, optionally runs regional
artifact detection per sample, and writes scop-style pass/fail metadata
that can be visualized with
[`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

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
  workers = 1,
  verbose = TRUE,
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

  Name of the Seurat spatial image. If `NULL`, the first image is used
  when present.

- sample.by:

  Optional metadata column identifying samples or images. If `NULL`, all
  spots are treated as one sample.

- metrics:

  QC metrics used by
  [`SpotSweeper::localOutliers()`](https://rdrr.io/pkg/SpotSweeper/man/localOutliers.html).
  If `NULL`, `nCount_<assay>`, `nFeature_<assay>`, and `percent.mito`
  are used.

- directions:

  Outlier direction for each metric. If `NULL`, count and feature
  metrics use `"lower"` and mitochondrial metrics use `"higher"`.

- n_neighbors:

  Number of nearest spatial neighbors for local outlier detection.

- cutoff:

  Modified z-score cutoff passed to local outlier detection.

- log:

  Whether SpotSweeper should log1p-transform local outlier metrics.

- run_artifact:

  Whether to run
  [`SpotSweeper::findArtifacts()`](https://rdrr.io/pkg/SpotSweeper/man/findArtifacts.html)
  per sample.

- mito_pattern:

  Regex prefixes used to identify mitochondrial genes.

- mito_gene:

  Optional explicit mitochondrial gene vector. When provided,
  `mito_pattern` is ignored.

- mito_percent:

  Metadata column used as mitochondrial percent for artifact detection.
  If `NULL`, `percent.mito` is used.

- mito_sum:

  Metadata column used as mitochondrial counts for artifact detection.
  If `NULL`, mitochondrial counts are computed.

- n_order, shape:

  Parameters passed to
  [`SpotSweeper::findArtifacts()`](https://rdrr.io/pkg/SpotSweeper/man/findArtifacts.html).

- prefix:

  Prefix used for metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- return_filtered:

  Whether to return only spots passing SpotSweeper QC.

- store_results:

  Whether to store detailed results in `srt@tools`.

- workers:

  Number of workers passed to SpotSweeper local outlier detection.

- verbose:

  Whether to print progress messages.

- ...:

  Additional named arguments passed to matching SpotSweeper backend
  functions when those arguments are supported by the installed version.

## Value

A `Seurat` object with SpotSweeper QC metadata. When
`store_results = TRUE`, detailed results are stored in
`srt@tools[[tool_name]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)

spatial <- RunSpotSweeper(
  visium_human_pancreas_sub,
  assay = "Spatial",
  n_neighbors = 18,
  n_order = 2
)

SpatialSpotPlot(spatial, group.by = "SpotSweeper_QC")
SpatialSpotPlot(spatial, group.by = "SpotSweeper_nCount_Spatial_z")
} # }
```
