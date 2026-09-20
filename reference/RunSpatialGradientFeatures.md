# Run spatial gradient feature screening

Screen features for expression trends along a spatial trajectory or with
distance from annotated regions.

## Usage

``` r
RunSpatialGradientFeatures(
  object,
  reference = c("trajectory", "annotation"),
  backend = "cpp",
  result_name = NULL,
  assay = NULL,
  layer = "data",
  variables = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  start = NULL,
  end = NULL,
  traj_df = NULL,
  annotation.by = NULL,
  annotation.groups = NULL,
  annotation.variable = NULL,
  annotation.threshold = NULL,
  unit = NULL,
  sign_var = "fdr",
  sign_threshold = 0.05,
  n_random = 10000,
  seed = 123,
  n_bins = 50,
  min_spots = 3,
  nfeatures = 2000,
  set_variable_features = FALSE,
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- reference:

  Spatial reference type: `"trajectory"` for STS or `"annotation"` for
  SAS.

- backend:

  Computation backend. The supported backend is `"cpp"`, which uses the
  compiled native spatial implementation.

- result_name:

  Name used to store this result. If `NULL`, a name is generated from
  `reference`.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer used for expression values.

- variables:

  Numeric variables or genes tested by the native backend. If `NULL`,
  `srt@tools[["SpatialVariableFeatures"]]` is used first, then variable
  features, then all assay features.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used by the `"cpp"` backend when no image
  coordinates are available.

- start, end, traj_df:

  Trajectory geometry used by the native backend when
  `reference = "trajectory"`.

- annotation.by, annotation.groups:

  Metadata grouping used to create the native annotation mask.

- annotation.variable, annotation.threshold:

  Numeric variable and threshold used to create a native numeric
  annotation mask. Numeric thresholds are interpreted as
  `">{threshold}"`.

- unit:

  Optional coordinate-unit label stored in the result source. It does
  not transform coordinates.

- sign_var, sign_threshold:

  Column and cutoff used to select top gradient variables from the
  native significance table.

- n_random, seed:

  Permutation count and random seed used by the native screening
  backend.

- n_bins:

  Number of distance bins used for the `"cpp"` backend screening curve.

- min_spots:

  Minimum number of non-zero spots required for a variable in the
  `"cpp"` backend.

- nfeatures:

  Number of top gradient variables retained in `top_variables` and
  optionally set as Seurat variable features.

- set_variable_features:

  Whether to set top gradient variables as Seurat variable features
  independently of `store_results`. If explicitly `TRUE`, an empty
  selection clears the target assay's variable features.

- store_results:

  Whether to store the normalized result in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- coordinate_space:

  Coordinate system used by the C++ distance calculations. The default
  is raw acquisition coordinates, so `start`, `end`, trajectory
  positions, and C++ distances share raw coordinate units. Use
  `"legacy_display"` explicitly for pre-0.9.0 display coordinates.
  Results retain the raw coordinate units used by the native backend.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object. When `store_results = TRUE`, spatial gradient
screening results are stored in
`srt@tools[["SpatialGradientFeatures"]]`. Otherwise existing stored
results are unchanged. Variable features are independently updated only
when `set_variable_features = TRUE`.

## Details

The current native interface is intentionally small. Older calls that
pass `sample_name`, `platform`, `img_scale_fct`, `assay_modality`,
`trajectory_id`, `width`, `annotation_id`, `core`, `distance`,
`angle_span`, `resolution`, `model_add`, `model_subset`, `model_remove`,
`control`, or `annotation_ids` now fail as unused arguments. Use
`image`, `coord.cols`, and `coordinate_space` for coordinate selection;
use `annotation.by`/`annotation.groups` or
`annotation.variable`/`annotation.threshold` for annotation references;
and use `n_random`, `n_bins`, `min_spots`, `sign_threshold`, and
`nfeatures` for native screening controls.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
counts <- GetAssayData5(spatial, assay = "Spatial", layer = "counts")
gradient_features <- names(sort(Matrix::rowSums(counts), decreasing = TRUE))[
  seq_len(min(8L, nrow(counts)))
]
# Use a diagonal example axis without permutation testing.
spatial <- RunSpatialGradientFeatures(
  spatial,
  reference = "trajectory",
  backend = "cpp",
  result_name = "example_axis",
  variables = gradient_features,
  start = c(min(spatial$x), min(spatial$y)),
  end = c(max(spatial$x), max(spatial$y)),
  layer = "counts",
  coord.cols = c("x", "y"),
  n_random = 0,
  n_bins = 5,
  min_spots = 3,
  sign_threshold = 1,
  nfeatures = 4,
  verbose = FALSE
)

SpatialGradientPlot(spatial, plot_type = "line")

SpatialGradientPlot(spatial, plot_type = "model")
```
