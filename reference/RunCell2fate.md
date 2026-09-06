# Run Cell2fate RNA velocity analysis

Run the official Python Cell2fate model on raw spliced and unspliced
counts from a `Seurat` object. Cell2fate uses an isolated Python 3.9
environment because its upstream dependency stack is not compatible with
the default scvi-tools environment. Inputs, model files, posterior
output, logs, and a reproducible manifest are persisted under
`result_dir`.

A returned object can be passed back to the function for a matching
resumed run. Output previously recorded by the same `tool_name` and
`prefix` is replaced in the returned copy; unrelated metadata or tool
entries are never overwritten.

## Usage

``` r
RunCell2fate(
  srt,
  result_dir,
  spliced_assay = "spliced",
  unspliced_assay = "unspliced",
  spliced_layer = "counts",
  unspliced_layer = "counts",
  cluster.by,
  features = NULL,
  remove_clusters = NULL,
  cells_per_cluster = 100L,
  min_shared_counts = 10L,
  n_var_genes = 2000L,
  n_modules = NULL,
  model_params = list(),
  train_params = list(max_epochs = 500L, batch_size = 1000L, train_size = 1, lr = 0.01,
    accelerator = "auto"),
  posterior_params = list(num_samples = 30L, batch_size = NULL, use_gpu = FALSE,
    return_samples = FALSE),
  seed = 1L,
  envname = NULL,
  resume = TRUE,
  overwrite = FALSE,
  prefix = "Cell2fate",
  tool_name = "Cell2fate",
  store_velocity = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object containing raw spliced and unspliced counts.

- result_dir:

  Empty directory, or a directory owned by an earlier `RunCell2fate()`
  run, used to persist inputs, model files, posterior output,
  per-attempt logs, and the run manifest.

- spliced_assay, unspliced_assay:

  Assays containing raw spliced and unspliced counts.

- spliced_layer, unspliced_layer:

  Raw-count layers in the corresponding assays.

- cluster.by:

  Metadata column containing cell-state labels used for Cell2fate
  training-data selection.

- features:

  Optional features to consider before Cell2fate filtering and
  variable-gene selection.

- remove_clusters:

  Optional cluster labels to remove before training.

- cells_per_cluster:

  Maximum cells retained per cluster. Use `NULL` to retain every cell.
  Cells excluded by this sampling receive `NA` posterior values and
  `FALSE` in the generated `<prefix>_selected` metadata column.

- min_shared_counts:

  Minimum total shared spliced and unspliced counts required for a gene.

- n_var_genes:

  Number of variable genes retained for model fitting.

- n_modules:

  Number of Cell2fate modules. If `NULL`, the upstream
  `get_max_modules()` heuristic is used.

- model_params:

  Named arguments passed to `Cell2fate_DynamicalModel()`.

- train_params:

  Named arguments passed to the model `train()` method.

- posterior_params:

  Named arguments passed to `export_posterior()` through its
  `sample_kwargs` argument.

- seed:

  Random seed used by Python, NumPy, PyTorch, and scvi-tools.

- envname:

  Name of the isolated Cell2fate environment. If `NULL`,
  `"cell2fate_env"` is used.

- resume:

  Reuse a completed run only when its input fingerprint, parameters, and
  artifact hashes match.

- overwrite:

  Permit replacement of incompatible artifacts in an owned `result_dir`.

- prefix:

  Prefix used for Cell2fate metadata columns. Existing columns without
  matching Cell2fate provenance are rejected.

- tool_name:

  Name of the `srt@tools` result entry. An existing unrelated entry is
  rejected.

- store_velocity:

  Whether to write the dense posterior velocity matrix to CSV and read
  it into `srt@tools`. The posterior `.h5ad`, including its velocity
  layer, is always retained on disk.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with Cell2fate time, uncertainty, module activation,
module-state, and training-cell-selection metadata. Cells not selected
for training have `NA` posterior values. Detailed provenance and
optional velocity values are stored in `srt@tools[[tool_name]]`.

## See also

[RunStandardWorkflow](https://mengxu98.github.io/scop/reference/RunStandardWorkflow.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunCell2fate(
  pancreas_sub,
  result_dir = "pancreas_cell2fate",
  cluster.by = "SubCellType",
  n_modules = 10
)
FeatureDimPlot(
  pancreas_sub,
  c("Cell2fate_time", "Cell2fate_time_uncertainty")
)
FeatureDimPlot(
  pancreas_sub,
  grep(
    "^Cell2fate_module_.*_activation$",
    colnames(pancreas_sub@meta.data),
    value = TRUE
  )
)
CellDimPlot(pancreas_sub, group.by = "Cell2fate_module_0_state")
} # }
```
