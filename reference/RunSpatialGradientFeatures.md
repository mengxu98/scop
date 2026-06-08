# Run spatial gradient feature screening

Run spatial trajectory or annotation gradient screening for Seurat
objects. The native `"cpp"` backend avoids SPATA2 object construction
for fast distance-based screening, while the `"spata2"` backend keeps
full upstream SPATA2 SAS/STS behavior. Results are normalized into plain
data.frames and stored in `srt@tools[["SpatialGradientFeatures"]]`; the
SPATA2 object itself is never stored.

## Usage

``` r
RunSpatialGradientFeatures(
  srt,
  reference = c("trajectory", "annotation"),
  backend = c("cpp", "spata2"),
  result_name = NULL,
  spata_object = NULL,
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
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- reference:

  Spatial reference type: `"trajectory"` for STS or `"annotation"` for
  SAS.

- backend:

  Computation backend. `"cpp"` uses SCOP's native fast spatial gradient
  implementation and avoids SPATA2 object construction. `"spata2"` uses
  SPATA2 directly for full upstream SAS/STS behavior.

- result_name:

  Name used to store this result. If `NULL`, a name is generated from
  `reference`.

- spata_object:

  Optional pre-built SPATA2 object. If `NULL`, `srt` is converted with
  `SPATA2::asSPATA2()`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used for expression values.

- variables:

  Numeric variables or genes passed to SPATA2. If `NULL`,
  `srt@misc[["SpatialVariableFeatures"]]` is used first, then variable
  features, then all assay features.

- sample_name, platform, img_scale_fct, assay_modality:

  Arguments forwarded to `SPATA2::asSPATA2()` when `spata_object` is
  `NULL`.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the native `"cpp"` backend when no
  image coordinates are available.

- trajectory_id, start, end, traj_df, width:

  Trajectory setup passed to `SPATA2::addSpatialTrajectory()` and
  `SPATA2::spatialTrajectoryScreening()`.

- annotation_ids:

  Existing SPATA2 spatial annotation ids. If `NULL`, annotations are
  created from `annotation.by` and `annotation.groups`, or from
  `annotation.variable` and `annotation.threshold`.

- annotation.by, annotation.groups:

  Metadata grouping used to create SPATA2 group annotations.

- annotation.variable, annotation.threshold:

  Numeric variable and threshold used to create SPATA2 numeric
  annotations. Numeric thresholds are interpreted as `">{threshold}"`.

- annotation_id:

  Base id used when creating annotations.

- core, distance, angle_span:

  SAS parameters forwarded to SPATA2.

- resolution, unit, sign_var, sign_threshold, model_add, model_subset,
  model_remove, n_random, seed, control:

  SPATA2 screening parameters.

- n_bins:

  Number of distance bins used for the native `"cpp"` backend screening
  curve.

- min_spots:

  Minimum number of non-zero spots required for a variable in the native
  `"cpp"` backend.

- nfeatures:

  Number of top gradient variables retained in `top_variables` and
  optionally set as Seurat variable features.

- set_variable_features:

  Whether to set top gradient variables as Seurat variable features.

- store_results:

  Whether to store the normalized result in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments forwarded to the SPATA2 screening function.

## Value

A `Seurat` object with spatial gradient screening results stored in
`srt@tools[["SpatialGradientFeatures"]]`.

## Examples

``` r
if (FALSE) { # \dontrun{
spatial <- RunSpatialGradientFeatures(
  spatial,
  reference = "annotation",
  annotation.by = "BayesSpace_cluster",
  annotation.groups = "1",
  variables = spatial@misc[["SpatialVariableFeatures"]]
)
SpatialGradientPlot(spatial, plot_type = "combined")
} # }
```
