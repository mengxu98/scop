# Infer gene regulatory networks with GRNBoost2

Run GRNBoost2-like regulatory network inference and return a
standardized adjacency table with columns \`TF\`, \`target\`, and
\`importance\`.

## Usage

``` r
RunGRNBoost2(object, ...)

# S3 method for class 'Seurat'
RunGRNBoost2(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  backend = c("cpp", "python"),
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  exclude_self = TRUE,
  importance_norm_power = 0,
  correlation_fill = FALSE,
  boost_weight = 0.5,
  correlation_weight = 0.8,
  covariance_weight = 0.4,
  correlation_method = c("pearson", "spearman"),
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grnboost2",
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
  cores = 1,
  seed = 1234,
  force = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'matrix'
RunGRNBoost2(object, ...)

# Default S3 method
RunGRNBoost2(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  backend = c("cpp", "python"),
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  exclude_self = TRUE,
  importance_norm_power = 0,
  correlation_fill = FALSE,
  boost_weight = 0.5,
  correlation_weight = 0.8,
  covariance_weight = 0.4,
  correlation_method = c("pearson", "spearman"),
  seed = 1234,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grnboost2",
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object or expression matrix.

- ...:

  Additional backend-specific arguments.

- assay:

  Assay used when \`object\` is a Seurat object.

- layer:

  Assay layer used when \`object\` is a Seurat object.

- regulators:

  Candidate transcription factor genes.

- targets:

  Optional target genes. If \`NULL\`, all genes in the GRN matrix are
  considered as candidate targets.

- backend:

  Runtime backend. Supports \`"cpp"\` and \`"python"\`.

- max_edges_per_target:

  Maximum incoming regulator edges retained per target. The default
  \`Inf\` keeps all positive-importance links, matching arboreto
  GRNBoost2 output.

- n_rounds:

  Number of boosting rounds for GRNBoost2-like tree ensemble inference.
  The default follows arboreto \`SGBM_KWARGS\`.

- learning_rate:

  GRNBoost2-like tree ensemble learning rate.

- max_depth:

  Maximum depth of each regression tree.

- max_features:

  Fraction of candidate regulators sampled at each tree split.

- subsample:

  Fraction of cells sampled for each boosting round.

- early_stop_window_length:

  Out-of-bag improvement window used for GRNBoost2 early stopping.

- exclude_self:

  Whether native GRNBoost2-like inference excludes a target gene from
  its own regulator feature set.

- importance_norm_power:

  Power used to normalize native GRNBoost2-like edge importance by the
  total importance for each target. Set to \`0\` to disable
  normalization.

- correlation_fill:

  Whether native GRNBoost2-like inference should fill missing TF-target
  candidates with expression-correlation support before applying
  \`max_edges_per_target\`.

- boost_weight:

  Weight of the native boosting score used when correlation fill is
  enabled.

- correlation_weight:

  Weight of the correlation-fill score relative to native boosting
  importance.

- covariance_weight:

  Weight of the covariance score used when correlation fill is enabled.

- correlation_method:

  Correlation method used for native GRNBoost2-like candidate filling.

- output_file:

  Optional path where the adjacency table is written.

- work_dir:

  Working directory used by Python backends.

- prefix:

  Prefix for temporary Python backend files.

- envname:

  Python environment used by Python backends.

- conda:

  Conda-compatible executable used by Python backends.

- prepare_env:

  Whether to prepare the Python environment before running Python
  backends.

- cores:

  Number of workers used by native/Python GRNBoost2.

- seed:

  Random seed passed to supported backends.

- force:

  Whether to rebuild existing \`output_file\`.

- verbose:

  Whether to print progress messages.

- genes_in:

  Matrix orientation for matrix inputs. \`"rows"\` means genes x cells;
  \`"columns"\` means cells x genes.

## Value

A data frame with columns \`TF\`, \`target\`, and \`importance\`.
