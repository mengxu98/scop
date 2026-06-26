# Run scMalignantFinder malignant cell identification

Run the Python package `scMalignantFinder` on a Seurat or AnnData object
and append malignant-cell predictions to Seurat metadata. The pretrained
model files are not bundled with `scop`; provide a directory containing
`model.joblib` and `ordered_feature.tsv` through `pretrain_dir`.

## Usage

``` r
RunscMalignantFinder(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "data",
  cells = NULL,
  pretrain_dir = NULL,
  train_h5ad_path = NULL,
  feature_path = NULL,
  model_method = c("LogisticRegression", "RandomForest", "XGBoost"),
  norm_type = FALSE,
  use_raw = FALSE,
  n_thread = 1,
  prefix = "",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- adata:

  Optional Python AnnData object.

- h5ad:

  Optional path to an `.h5ad` file.

- assay:

  Assay used when `srt` is supplied. Default is `"RNA"`.

- layer:

  Layer used when `srt` is supplied. Default is `"data"`.

- cells:

  Optional cells to run. If supplied with `srt`, results are appended to
  these cells and other cells receive `NA`.

- pretrain_dir:

  Directory containing pretrained `scMalignantFinder` model files:
  `model.joblib` and `ordered_feature.tsv`.

- train_h5ad_path:

  Optional training `.h5ad` file used when training a model from
  scratch.

- feature_path:

  Optional feature file used when training from scratch.

- model_method:

  Model used when training from scratch. One of `"LogisticRegression"`,
  `"RandomForest"`, or `"XGBoost"`.

- norm_type:

  Passed to `scMalignantFinder`. Use `TRUE` for raw counts that should
  be library-size normalized; use `FALSE` for already normalized input.
  Default is `FALSE`.

- use_raw:

  Whether to use `adata.raw.X` when available.

- n_thread:

  Number of threads used by `scMalignantFinder`.

- prefix:

  Optional prefix for output metadata columns. Default preserves the
  original `scMalignantFinder` column names.

- return_seurat:

  Whether to return a Seurat object when `srt` is supplied. If `FALSE`,
  returns a data frame of predictions.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with `scMalignantFinder_prediction` and
`malignancy_probability` added, or a data frame when
`return_seurat = FALSE`.

## References

Yu Q, Li YY, Chen Y. scMalignantFinder distinguishes malignant cells in
single-cell and spatial transcriptomics by leveraging cancer signatures.
Communications Biology. 2025.

## See also

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md),
[RunscMalignantRegion](https://mengxu98.github.io/scop/reference/RunscMalignantRegion.md),
[RunscMalignantStates](https://mengxu98.github.io/scop/reference/RunscMalignantStates.md)

## Examples

``` r
if (FALSE) {
data(pancreas_sub)
pancreas_sub <- RunscMalignantFinder(
  pancreas_sub,
  assay = "RNA",
  layer = "data",
  pretrain_dir = "path/to/pretrained_model",
  norm_type = FALSE
)
CellDimPlot(pancreas_sub, group.by = "malignancy_probability")
}
```
