# Run MERINGUE spatial autocorrelation analysis

Run `MERINGUE` spatial autocorrelation, spatial cross-correlation, and
spatial module analysis for a spatial `Seurat` object.

## Usage

``` r
RunMERINGUE(
  srt,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("col", "row"),
  features = NULL,
  mode = c("autocorrelation", "cross_correlation", "modules"),
  nfeatures = 2000,
  min_spots = 5,
  filterDist = NA_real_,
  binary = TRUE,
  alternative = "greater",
  nperm = 0,
  ncores = 1,
  pairwise_features = NULL,
  set_variable_features = FALSE,
  store_results = TRUE,
  verbose = TRUE,
  seed = 11,
  neighbor_params = list(),
  moran_params = list(),
  cross_cor_params = list(),
  module_params = list()
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used for expression values.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

- features:

  Features to score. If `NULL`, current variable features are used; if
  no variable features are present, all assay features are used.

- mode:

  MERINGUE analysis modes to run. `"autocorrelation"` computes spatial
  autocorrelation, `"cross_correlation"` computes pairwise spatial
  cross-correlation, and `"modules"` detects spatial gene modules.

- nfeatures:

  Number of top spatial features stored in
  `srt@misc[["SpatialVariableFeatures"]]`.

- min_spots:

  Minimum number of spots with non-zero expression required for a
  feature to be tested.

- filterDist:

  Euclidean distance cutoff passed to
  [`MERINGUE::getSpatialNeighbors()`](https://rdrr.io/pkg/MERINGUE/man/getSpatialNeighbors.html).

- binary:

  Whether to binarize the MERINGUE spatial neighbor matrix.

- alternative:

  Alternative hypothesis passed to MERINGUE Moran tests.

- nperm:

  Number of label permutations used for empirical p values. The default
  `0` skips p-value calculation.

- ncores:

  Number of cores passed to MERINGUE permutation tests.

- pairwise_features:

  Features used for spatial cross-correlation. If `NULL`, top spatially
  autocorrelated features are used.

- set_variable_features:

  Whether to set the top spatial features as variable features for
  `assay`.

- store_results:

  Whether to store the full result in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed used for permutation tests.

- neighbor_params, moran_params, cross_cor_params, module_params:

  Named lists of additional arguments passed to the corresponding
  MERINGUE steps.

## Value

A `Seurat` object with MERINGUE results stored in
`srt@tools[["MERINGUE"]]` and top autocorrelated features stored in
`srt@misc[["MERINGUEFeatures"]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)

spatial <- RunMERINGUE(
  spatial,
  assay = "Spatial",
  mode = c("autocorrelation", "cross_correlation"),
  nfeatures = 50
)

head(spatial@tools[["MERINGUE"]]$autocorrelation)
SpatialSpotPlot(
  spatial,
  features = spatial@misc[["MERINGUEFeatures"]][1:2]
)
} # }
```
