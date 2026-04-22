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

  A Seurat object. Default is `NULL`. If provided, `adata` will be
  ignored.

- adata:

  An anndata object. Default is `NULL`.

- assay:

  Which assay to use. Default is `"RNA"`.

- layer:

  Which layer to use. Default is `data`.

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

  Whether to use GPU for over-clustering. Default is `FALSE`.

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

  Return mode. One of `"summary"`, `"path"`, or `"model"`. Default is
  `"summary"`.

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
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)

model_info <- TrainCellTypist(
  srt = pancreas_sub,
  labels = "SubCellType",
  model_path = tempfile(fileext = ".pkl")
)

data(panc8_sub)
genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames
)
panc8_sub <- standard_scop(panc8_sub)

pancreas_sub <- RunCellTypist(
  srt = pancreas_sub,
  model = model_info
)
CellDimPlot(
  pancreas_sub,
  group.by = c("SubCellType", "celltypist_predicted_labels")
)

panc8_sub <- RunCellTypist(
  srt = panc8_sub,
  model = model_info
)
CellDimPlot(
  panc8_sub,
  group.by = c("celltype", "celltypist_predicted_labels")
)

ht <- CellCorHeatmap(
  srt_query = panc8_sub,
  srt_ref = panc8_sub,
  query_group = "celltypist_predicted_labels",
  ref_group = "celltype",
  width = 4,
  height = 3
)
ht$plot
} # }
```
