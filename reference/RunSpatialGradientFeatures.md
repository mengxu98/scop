# Run spatial gradient feature screening

Run native spatial trajectory or annotation gradient screening for
Seurat objects. The compiled C++ backend computes distance-based
screening and stores validated result tables in `srt@tools` when
requested.

## Usage

``` r
RunSpatialGradientFeatures(
  srt,
  reference = c("trajectory", "annotation"),
  backend = "cpp",
  result_name = NULL,
  assay = NULL,
  layer = "data",
  variables = NULL,
  sample_name = NULL,
  platform = "Undefined",
  image = NULL,
  coord.cols = c("x", "y"),
  img_scale_fct = "lowres",
  assay_modality = "gene",
  trajectory_id = "scop_gradient",
  start = NULL,
  end = NULL,
  traj_df = NULL,
  width = NULL,
  annotation_ids = NULL,
  annotation.by = NULL,
  annotation.groups = NULL,
  annotation.variable = NULL,
  annotation.threshold = NULL,
  annotation_id = "scop_gradient",
  core = FALSE,
  distance = "dte",
  angle_span = c(0, 360),
  resolution = NULL,
  unit = NULL,
  sign_var = "fdr",
  sign_threshold = 0.05,
  model_add = NULL,
  model_subset = NULL,
  model_remove = NULL,
  n_random = 10000,
  seed = 123,
  control = NULL,
  n_bins = 50,
  min_spots = 3,
  nfeatures = 2000,
  set_variable_features = FALSE,
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
)
```

## Arguments

- srt:

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

- sample_name, platform, img_scale_fct, assay_modality, trajectory_id,
  width, annotation_id, core, distance, angle_span, resolution,
  model_add, model_subset, model_remove, control:

  Legacy compatibility inputs. The native backend ignores these values
  and does not include them in the stored effective-parameter summary.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used by the `"cpp"` backend when no image
  coordinates are available.

- start, end, traj_df:

  Trajectory geometry used by the native backend when
  `reference = "trajectory"`.

- annotation_ids:

  Reserved compatibility input. It is not supported by the native
  backend; use `annotation.by` with `annotation.groups`, or
  `annotation.variable` with `annotation.threshold`.

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

- ...:

  Additional named arguments accepted for compatibility. The native
  backend ignores them and does not store them as effective parameters.

## Value

A `Seurat` object. When `store_results = TRUE`, spatial gradient
screening results are stored in
`srt@tools[["SpatialGradientFeatures"]]`. Otherwise existing stored
results are unchanged. Variable features are independently updated only
when `set_variable_features = TRUE`.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
spatial <- RunSpatialGradientFeatures(
  spatial,
  reference = "trajectory",
  backend = "cpp",
  result_name = "ductal_axis",
  variables = rownames(spatial)[1:8],
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

SpatialGradientPlot(
  spatial,
  plot_type = "surface",
  features = rownames(spatial)[1:2],
  overlay_image = FALSE,
  coord.cols = c("x", "y"),
  pt.size = 1.2
)
```
