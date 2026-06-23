# Run tAge transcriptomic aging-clock prediction

Run tAge transcriptomic aging-clock prediction

## Usage

``` r
RuntAge(
  object,
  model_paths = NULL,
  backend = c("python", "cpp"),
  clock = c("Chronoage", "NormalizedAge", "Mortality"),
  model_species = c("auto", "Multispecies", "Mouse", "Rodents"),
  model_tissue = "Multitissue",
  model_preprocessing = "scaled_diff",
  model_cache_dir = NULL,
  datasets_base_url = "https://raw.githubusercontent.com/mengxu98/datasets/main/tAge/EN",
  zenodo_record = "18763485",
  max_model_size = 1024^3,
  metadata = NULL,
  group.by = NULL,
  split.by = NULL,
  assay = NULL,
  layer = "counts",
  species = c("mouse", "human", "rat", "monkey"),
  mode = c("EN", "BR"),
  gene_mapping_type = c("Gene.Symbol", "Ensembl"),
  coverage_threshold = 1e+06,
  count_threshold = 10,
  percent_threshold = 20,
  control_group_column = NULL,
  control_group_label = NULL,
  remove_outliers = FALSE,
  outlier.by = split.by,
  outlier_n_components = 10,
  outlier_threshold_quantile = 0.99,
  outlier_min_samples = 10,
  min_samples = 5,
  shuffle = FALSE,
  seed = NULL,
  check_python = TRUE,
  tool_name = "tAge",
  store_eset = FALSE,
  store_processed = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object, `ExpressionSet`, or expression matrix-like object
  with genes in rows and pseudobulk/bulk samples in columns.

- model_paths:

  Named list or character vector of local tAge model files (`.rds` for
  `backend = "cpp"`, `.pkl` for `backend = "python"`). Names must match
  one or more of `"scaled"`, `"scaled_diff"`, `"yugene"`, and
  `"yugene_diff"`. If `NULL`, model files are downloaded from
  `mengxu98/datasets` for the C++ backend or Zenodo for the Python
  backend.

- backend:

  Model prediction backend. `"cpp"` uses converted Elastic Net models
  and C++ prediction. `"python"` uses `tAge::predict_tAge()`.

- clock:

  tAge clock family used when `model_paths = NULL`.

- model_species:

  Model species scope used when `model_paths = NULL`. `"auto"` chooses
  `"Mouse"` for `species = "mouse"`, `"Rodents"` for `species = "rat"`,
  and `"Multispecies"` otherwise.

- model_tissue:

  Model tissue scope used when `model_paths = NULL`.

- model_preprocessing:

  Preprocessing-specific model file(s) used when `model_paths = NULL`.
  Values map to tAge model filename suffixes.

- model_cache_dir:

  Directory used to cache downloaded tAge models. If `NULL`, uses
  `tools::R_user_dir("scop", "data")`.

- datasets_base_url:

  Base URL or local directory containing converted tAge EN models from
  `mengxu98/datasets`.

- zenodo_record:

  Zenodo record ID for tAge model files.

- max_model_size:

  Maximum model file size, in bytes, allowed for automatic download.
  Increase this or use `Inf` for large Bayesian Ridge models.

- metadata:

  Sample metadata for matrix input. Row names must match columns of
  `object`. If `NULL`, minimal `sample_id` metadata is created.

- group.by:

  Metadata columns used to form pseudobulk groups for `Seurat` input. If
  `NULL`, cells are aggregated sequentially across the whole object.
  `split.by` and `control_group_column` are automatically included in
  the aggregation metadata when provided.

- split.by:

  Optional pseudobulk metadata column used to run tAge separately by
  group, such as tissue.

- assay, layer:

  Assay and layer used as raw counts for `Seurat` input.

- species:

  Species passed to tAge. Supported values are `"mouse"`, `"human"`,
  `"rat"`, and `"monkey"`.

- mode:

  tAge model mode: `"EN"` for Elastic Net or `"BR"` for Bayesian Ridge.

- gene_mapping_type:

  Gene identifier type passed to tAge preprocessing: `"Gene.Symbol"` or
  `"Ensembl"`.

- coverage_threshold:

  Minimum cumulative read count per pseudobulk sample for `Seurat`
  input.

- count_threshold, percent_threshold:

  Gene filtering thresholds passed to `tAge::tAge_preprocessing()`.

- control_group_column, control_group_label:

  Optional control group used by tAge control subtraction.

- remove_outliers:

  Whether to run `tAge::remove_outliers()` on the pseudobulk
  `ExpressionSet`.

- outlier.by:

  Optional metadata column used to split samples before outlier
  detection. Defaults to `split.by`.

- outlier_n_components, outlier_threshold_quantile, outlier_min_samples:

  Parameters passed to `tAge::remove_outliers()`.

- min_samples:

  Minimum pseudobulk samples per `split.by` group for
  `tAge::tAge_by_group()`.

- shuffle, seed:

  Whether to shuffle cells during pseudobulk aggregation and the random
  seed used for that shuffle.

- check_python:

  Whether to verify that the active `reticulate` Python can import
  `joblib`, `pandas`, and `sklearn` before tAge prediction.

- tool_name:

  Name of the Seurat tool entry used to store results.

- store_eset, store_processed:

  Whether to store the pseudobulk `ExpressionSet` and processed tAge
  `ExpressionSet` list in the returned Seurat tool entry or list result.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with tAge results stored in `object@tools`, or a list
with `predictions`, metadata, and parameters for non-Seurat input.

## References

Tyshkovskiy, A., Glubokov, D., Moliere, A., et al. (2026). Universal
transcriptomic hallmarks of mammalian ageing and mortality. *Nature*.
[doi:10.1038/s41586-026-10542-3](https://doi.org/10.1038/s41586-026-10542-3)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RuntAge(
  pancreas_sub,
  group.by = "CellType",
  species = "mouse",
  mode = "EN"
)
head(pancreas_sub@tools$tAge$predictions)
} # }
```
