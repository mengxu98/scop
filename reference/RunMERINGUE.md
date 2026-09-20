# Run MERINGUE spatial autocorrelation analysis

Run `MERINGUE` spatial autocorrelation, spatial cross-correlation, and
spatial module analysis for a spatial `Seurat` object.

## Usage

``` r
RunMERINGUE(
  object,
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
  cores = 1,
  ncores = NULL,
  pairwise_features = NULL,
  set_variable_features = FALSE,
  store_results = TRUE,
  verbose = TRUE,
  seed = 11,
  neighbor_params = list(),
  moran_params = list(),
  cross_cor_params = list(),
  module_params = list(),
  coordinate_space = c("raw", "legacy_display"),
  backend = c("cpp", "r"),
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer used for expression values.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- features:

  Features to score. If `NULL`, current variable features are used; if
  no variable features are present, all assay features are used.

- mode:

  MERINGUE analysis modes to run. `"autocorrelation"` computes spatial
  autocorrelation, `"cross_correlation"` computes pairwise spatial
  cross-correlation, and `"modules"` detects spatial gene modules.

- nfeatures:

  Number of top spatial features stored in
  `srt@tools[["SpatialVariableFeatures"]]`.

- min_spots:

  Minimum number of spots with non-zero expression required for a
  feature to be tested.

- filterDist:

  Euclidean distance cutoff passed to `MERINGUE::getSpatialNeighbors()`,
  expressed in the selected coordinate units.

- binary:

  Whether to binarize the MERINGUE spatial neighbor matrix.

- alternative:

  Alternative hypothesis passed to MERINGUE Moran tests.

- nperm:

  Number of label permutations used for empirical p values. The default
  `0` skips p-value calculation.

- cores:

  Number of cores passed to MERINGUE permutation tests. Mapped to
  `ncores` for the R backend and `n_threads` for the native backend.

- ncores:

  Deprecated alias for `cores`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

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

- coordinate_space:

  Coordinate system used for MERINGUE distances. The default is raw
  acquisition coordinates; `"legacy_display"` remains an explicit
  compatibility option.

- backend:

  Implementation used for the Moran autocorrelation tests. `"cpp"` uses
  a scop-compiled kernel that is numerically equivalent to
  `MERINGUE::moranTest()` and `MERINGUE::moranPermutationTest()` and is
  substantially faster for permutation tests; `"r"` calls the original
  MERINGUE functions directly. When `"cpp"` is selected, `moran_params`
  are ignored with a warning because the compiled kernel has no extra
  arguments.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with MERINGUE results stored in
`srt@tools[["MERINGUE"]]`. Top autocorrelated features are available at
`srt@tools[["MERINGUE"]]$summary$top_features`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
visium_human_pancreas_sub <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial", verbose = FALSE
)
spatial <- RunMERINGUE(
  visium_human_pancreas_sub,
  assay = "Spatial", coord.cols = c("x", "y"),
  mode = "autocorrelation", nfeatures = 20, verbose = FALSE
)
spatial@tools$MERINGUE$summary
} # }
```
