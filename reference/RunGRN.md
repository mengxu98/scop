# Infer gene regulatory networks with a selected backend method

Unified GRN inference entry point for GRNBoost2, GENIE3, RegDiffusion,
and GNIPLR. It returns a standardized adjacency table with at least
\`TF\`, \`target\`, and \`importance\` columns.

## Usage

``` r
RunGRN(object, ...)

# S3 method for class 'Seurat'
RunGRN(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
  backend = c("cpp", "python"),
  ...
)

# S3 method for class 'matrix'
RunGRN(object, ...)

# Default S3 method
RunGRN(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
  backend = c("cpp", "python"),
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grn",
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  exclude_self = TRUE,
  correlation_threshold = 0.3,
  lasso_degree = 30,
  lasso_alpha = 0.1,
  max_lag = 3,
  envname = NULL,
  conda = "auto",
  cores = 1,
  seed = 1234,
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

  Optional target genes. If \`NULL\`, all genes are considered.

- grn_method:

  GRN inference method.

- backend:

  Runtime backend. \`"cpp"\` is available for GRNBoost2 and GNIPLR;
  \`"python"\` is required for RegDiffusion.

- genes_in:

  Matrix orientation for matrix inputs. \`"rows"\` means genes x cells;
  \`"columns"\` means cells x genes.

- output_file:

  Optional path where the adjacency table is written.

- work_dir:

  Working directory used by Python backends.

- prefix:

  Prefix for temporary backend files.

- max_edges_per_target:

  Maximum incoming regulator edges retained per target.

- n_rounds:

  Number of boosting rounds for GRNBoost2-like inference.

- learning_rate:

  GRNBoost2-like tree ensemble learning rate.

- max_depth:

  Maximum depth of each regression tree.

- max_features:

  Fraction of candidate regulators sampled at each split.

- subsample:

  Fraction of cells sampled for each boosting round.

- early_stop_window_length:

  Out-of-bag improvement window used for GRNBoost2 early stopping.

- exclude_self:

  Whether GRNBoost2-like inference excludes a target gene from its own
  regulator feature set.

- correlation_threshold:

  Relative correlation filter used by GNIPLR.

- lasso_degree:

  Polynomial degree used by GNIPLR.

- lasso_alpha:

  LASSO regularization strength used by GNIPLR.

- max_lag:

  Maximum lag used by GNIPLR.

- envname:

  Python environment used by Python backends. If \`NULL\`, GNIPLR uses
  the default \`"scop_env"\` through \[RunGNIPLR()\], while pySCENIC
  GRNBoost2 and RegDiffusion use the isolated \`"scenic_env"\`
  environment.

- conda:

  Conda-compatible executable used by Python backends.

- cores:

  Number of workers used by supported methods.

- seed:

  Random seed passed to supported backends.

- force:

  Whether to rebuild existing \`output_file\`.

- verbose:

  Whether to print progress messages.

## Value

A data frame with standardized GRN edges.

## Examples

``` r
data(pancreas_sub)
expr <- GetAssayData5(
  pancreas_sub,
  assay = SeuratObject::DefaultAssay(pancreas_sub),
  layer = "counts"
)
expr <- as.matrix(expr[, seq_len(8)])
expr <- expr[
  names(sort(apply(expr, 1, stats::var), decreasing = TRUE))[seq_len(5)],
]

gniplr_grn <- RunGRN(
  expr,
  genes_in = "rows",
  grn_method = "gniplr",
  backend = "cpp",
  correlation_threshold = 0,
  lasso_degree = 1,
  max_lag = 1,
  max_edges_per_target = 2,
  verbose = FALSE
)
```
