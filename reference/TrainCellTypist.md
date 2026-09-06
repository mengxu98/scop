# Train a CellTypist model

Train a CellTypist model from a Seurat object, AnnData object, or h5ad
file.

## Usage

``` r
TrainCellTypist(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "data",
  labels,
  genes = NULL,
  transpose_input = FALSE,
  with_mean = TRUE,
  check_expression = TRUE,
  C = 1,
  solver = NULL,
  max_iter = NULL,
  n_jobs = 1,
  use_SGD = FALSE,
  alpha = 1e-04,
  use_GPU = FALSE,
  mini_batch = FALSE,
  batch_number = 100,
  batch_size = 1000,
  epochs = 10,
  balance_cell_type = FALSE,
  feature_selection = FALSE,
  top_genes = 300,
  date = "",
  details = "",
  url = "",
  source = "",
  version = "",
  model_path = NULL,
  return = c("summary", "path", "model"),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object. If provided, `adata` will be ignored.

- adata:

  Optional AnnData object used as input.

- h5ad:

  Optional path to an input `.h5ad` file.

- assay:

  Which assay to use.

- layer:

  Assay layer to use.

- labels:

  Cell labels used for training. Can be a metadata column name when
  `srt`/`adata` is supplied, or a vector aligned with cells.

- genes:

  Optional gene names. Usually inferred from the input object.

- transpose_input:

  Whether to transpose the input matrix before training.

- with_mean:

  Whether to center features during scaling.

- check_expression:

  Whether to validate expected CellTypist input format.

- C:

  Inverse regularization strength for logistic regression.

- solver:

  Optional solver passed to CellTypist.

- max_iter:

  Optional maximum iterations.

- n_jobs:

  Number of CPUs used by CellTypist training.

- use_SGD:

  Whether to use SGD training.

- alpha:

  Regularization strength for SGD training.

- use_GPU:

  Whether to use GPU for over-clustering.

- mini_batch:

  Whether to enable mini-batch training.

- batch_number:

  Number of batches per epoch for mini-batch training.

- batch_size:

  Batch size for mini-batch training.

- epochs:

  Number of epochs for mini-batch training.

- balance_cell_type:

  Whether to balance cell types during mini-batch training.

- feature_selection:

  Whether to run CellTypist feature selection.

- top_genes:

  Number of top genes used during feature selection.

- date, details, url, source, version:

  Free-text metadata stored in the model.

- model_path:

  Optional output path for the trained model.

- return:

  Return mode. One of `"summary"`, `"path"`, or `"model"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Depends on `return`: a summary list, the model path, or a Python model
object.

## See also

[RunCellTypist](https://mengxu98.github.io/scop/reference/RunCellTypist.md),
[CellTypistModels](https://mengxu98.github.io/scop/reference/CellTypistModels.md)

## Examples

``` r
if (FALSE) { # \dontrun{
check_python(c("celltypist", "anndata"))
data(pbmcmultiome_sub)
pbmcmultiome_sub <- RunStandardWorkflow(
  pbmcmultiome_sub, assay = "RNA", linear_reduction_dims = 10
)
model_info <- TrainCellTypist(
  srt = pbmcmultiome_sub, labels = "CellType",
  model_path = tempfile(fileext = ".pkl")
)
} # }
```
