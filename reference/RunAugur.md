# Prioritize perturbed cell types using Augur

Prioritize perturbed cell types using Augur

## Usage

``` r
RunAugur(
  srt,
  celltype.by,
  label.by,
  assay = NULL,
  layer = "counts",
  backend = c("cpp", "r"),
  features = NULL,
  n_subsamples = 50,
  subsample_size = 20,
  folds = 3,
  min_cells = NULL,
  var_quantile = 0.5,
  feature_perc = 0.5,
  cores = 1,
  select_var = TRUE,
  augur_mode = c("default", "velocity", "permute"),
  classifier = c("rf", "lr"),
  rf_params = list(trees = 100, mtry = 2, min_n = NULL, importance = "accuracy"),
  lr_params = list(mixture = 1, penalty = "auto"),
  prefix = "augur",
  tool_name = "Augur",
  add_meta = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- celltype.by:

  Metadata column or vector defining the cell types used by Augur.

- label.by:

  Metadata column or vector defining the labels to predict, such as
  condition, stimulation, sample type, or technology.

- assay:

  Assay used to extract the feature matrix. If `NULL`, the default assay
  is used.

- layer:

  Assay layer used as the feature matrix.

- backend:

  Backend used to run Augur. `"cpp"` uses a parity-preserving `scop`
  backend that keeps Augur's sampling, random forest, and metric
  semantics but removes the heavier recipe/bake loop and uses C++ sparse
  matrix preparation where possible. `"r"` calls the native
  `Augur::calculate_auc` implementation.

- features:

  Features used by Augur. If `NULL`, all features in `assay` are used.

- n_subsamples, subsample_size, folds, min_cells, var_quantile,
  feature_perc, select_var, augur_mode, classifier, rf_params,
  lr_params:

  Arguments passed to `Augur::calculate_auc`.

- cores:

  Number of cores used by Augur.

- prefix:

  Prefix for metadata columns written to `srt@meta.data`.

- tool_name:

  Name of the `srt@tools` entry used to store Augur results.

- add_meta:

  Whether to write `prefix_auc` and `prefix_rank` metadata columns back
  to each cell by matching `celltype.by`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to `Augur::calculate_auc`.

## Value

A `Seurat` object with native Augur results stored in
`srt@tools[[tool_name]]`. When `add_meta = TRUE`, cell-level metadata
columns `prefix_auc` and `prefix_rank` are added for use with existing
`scop` plotting functions such as
[`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md).

## References

Skinnider, M.A., Squair, J.W., Kathe, C., et al. (2021). Cell type
prioritization in single-cell data. *Nature Biotechnology*, 39, 30-34.
[doi:10.1038/s41587-020-0605-1](https://doi.org/10.1038/s41587-020-0605-1)

## Examples

``` r
data(panc8_sub)
panc8_sub <- subset(panc8_sub, subset = tech %in% c("celseq", "celseq2"))
panc8_sub <- standard_scop(panc8_sub, verbose = FALSE)
panc8_sub <- RunAugur(
  panc8_sub,
  celltype.by = "celltype",
  label.by = "tech",
  n_subsamples = 5,
  subsample_size = 20,
  min_cells = 20,
  cores = 1,
  verbose = FALSE,
  rf_params = list(
    trees = 20,
    mtry = 2,
    min_n = NULL,
    importance = "accuracy"
  )
)

panc8_sub@tools$Augur$AUC
FeatureDimPlot(
  panc8_sub,
  features = "augur_auc",
  reduction = "StandardUMAP2D",
  bg_cutoff = -Inf
)
```
