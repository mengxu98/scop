# Run cell2location spatial deconvolution

Estimate cell-type abundances and proportions in spatial transcriptomics
data using cell2location.

## Usage

``` r
RunCell2location(
  object,
  result_dir,
  reference = NULL,
  reference_label = "celltype",
  reference_signatures = NULL,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  spatial_batch = NULL,
  reference_batch = NULL,
  reference_covariates = NULL,
  min_cells = 5L,
  N_cells_per_location = 30,
  detection_alpha = 20,
  gene_filter_params = list(cell_count_cutoff = 5, cell_percentage_cutoff2 = 0.03,
    nonz_mean_cutoff = 1.12),
  reference_train_params = list(max_epochs = 250L, batch_size = 2500L, train_size = 1, lr
    = 0.002, accelerator = "auto", device = "auto"),
  spatial_train_params = list(max_epochs = 30000L, batch_size = NULL, train_size = 1,
    accelerator = "auto", device = "auto"),
  reference_posterior_params = list(num_samples = 1000L, batch_size = 2500L),
  spatial_posterior_params = list(batch_size = 2500L),
  envname = NULL,
  resume = TRUE,
  overwrite = FALSE,
  prefix = "Cell2location",
  tool_name = "Cell2location",
  store_results = TRUE,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  Spatial `Seurat` object containing raw counts.

- result_dir:

  Directory used to persist inputs, models, posterior results, tables,
  logs, and the run manifest.

- reference:

  Optional single-cell `Seurat` reference. Required when
  `reference_signatures` is `NULL`.

- reference_label:

  Metadata column in `reference` containing cell types.

- reference_signatures:

  Optional gene-by-cell-type numeric matrix, data.frame, or CSV file.
  When supplied, reference-model training is skipped.

- assay, reference_assay:

  Assays containing spatial and reference counts.

- layer, reference_layer:

  Raw-count layers.

- features:

  Optional features to consider before shared non-zero genes are
  selected.

- spatial_batch, reference_batch:

  Optional metadata columns describing spatial and reference batches.

- reference_covariates:

  Optional categorical reference metadata columns used to model
  technical effects.

- min_cells:

  Minimum number of reference cells retained per cell type.

- N_cells_per_location:

  Expected average number of cells per spatial location.

- detection_alpha:

  Prior controlling within-batch variation in RNA detection sensitivity.

- gene_filter_params:

  Named arguments passed to
  `cell2location.utils.filtering.filter_genes()`.

- reference_train_params, spatial_train_params:

  Named arguments passed to the reference and spatial model `train()`
  methods.

- reference_posterior_params:

  Named sampling arguments used when exporting reference signatures.

- spatial_posterior_params:

  Named sampling arguments used when exporting the spatial q05
  posterior.

- envname:

  Name of the shared SCOP Python environment.

- resume:

  Reuse completed stages only when their manifest matches the current
  inputs and parameters.

- overwrite:

  Permit replacement of incompatible existing artifacts.

- prefix:

  Prefix used for abundance/proportion metadata columns.

- tool_name:

  Name of the `srt@tools` result entry.

- store_results:

  Whether to store detailed result matrices and paths in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with cell2location abundance, proportion, dominant
cell type, and maximum-proportion metadata.

## Details

Requires raw spatial counts and an annotated single-cell reference or
gene-by-cell-type signatures with matching gene identifiers.

## See also

[`Cell2locationPlot()`](https://mengxu98.github.io/scop/reference/Cell2locationPlot.md),
[cell2location
tutorial](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html)
