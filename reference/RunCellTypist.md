# Run CellTypist cell type annotation

CellTypist is an automated cell type annotation tool for scRNA-seq
datasets based on logistic regression classifiers. This function runs
CellTypist annotation on a Seurat object or AnnData object.

## Usage

``` r
RunCellTypist(
  object = NULL,
  adata = NULL,
  assay = "RNA",
  layer = "data",
  model = "Immune_All_Low.pkl",
  mode = "best match",
  p_thres = 0.5,
  majority_voting = FALSE,
  over_clustering = NULL,
  min_prop = 0,
  use_GPU = FALSE,
  insert_labels = TRUE,
  insert_conf = TRUE,
  insert_conf_by = "predicted_labels",
  insert_prob = FALSE,
  insert_decision = FALSE,
  prefix = "celltypist_",
  return_seurat = !is.null(srt),
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A Seurat object. If provided, `adata` will be ignored.

- adata:

  Optional AnnData object used as input.

- assay:

  Which assay to use.

- layer:

  Assay layer to use.

- model:

  Model name or path. Supports three formats:

  1.  Model name (e.g., `"Immune_All_Low.pkl"`): automatically searched
      in `~/.celltypist/data/models/`

  2.  Full path (contains `/`): use the provided path directly

  3.  `NULL`: use default model

  4.  A summary list returned by
      [`TrainCellTypist()`](https://mengxu98.github.io/scop/reference/TrainCellTypist.md):
      use its `model_path`

- mode:

  Prediction mode: `"best match"` or `"prob match"`.

- p_thres:

  Probability threshold for `"prob match"` mode.

- majority_voting:

  Whether to use majority voting.

- over_clustering:

  Over-clustering result. Can be:

  - String: column name in Seurat metadata or AnnData obs

  - Vector: over-clustering labels

  - `NULL`: use heuristic over-clustering

- min_prop:

  Minimum proportion for majority voting.

- use_GPU:

  Whether to use GPU for over-clustering.

- insert_labels:

  Whether to insert predicted labels.

- insert_conf:

  Whether to insert confidence scores.

- insert_conf_by:

  Which prediction type to base confidence on.

- insert_prob:

  Whether to insert probability matrix.

- insert_decision:

  Whether to insert decision matrix.

- prefix:

  Prefix for inserted columns.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

An AnnData object or a Seurat object depending on the `return_seurat`
argument.

## See also

[CellTypistModels](https://mengxu98.github.io/scop/reference/CellTypistModels.md),
[TrainCellTypist](https://mengxu98.github.io/scop/reference/TrainCellTypist.md),
[RunSingleR](https://mengxu98.github.io/scop/reference/RunSingleR.md),
[RunScmap](https://mengxu98.github.io/scop/reference/RunScmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
check_python(c("celltypist", "anndata"))
data(pbmcmultiome_sub)
pbmcmultiome_sub <- RunStandardWorkflow(
  pbmcmultiome_sub,
  assay = "RNA", linear_reduction_dims = 10
)
pbmcmultiome_sub <- RunCellTypist(
  pbmcmultiome_sub,
  model = "Immune_All_Low.pkl", verbose = FALSE
)
} # }
```
