# Map query cells into a reference space

Reference mapping workflow modeled after Seurat's `MapQuery` pattern.
The current implementation computes gene activity for ATAC query cells,
finds anchors to an RNA reference, integrates the query into the
reference embedding space, transfers labels, and projects query cells to
the reference UMAP.

## Usage

``` r
RunReferenceMapping(
  srt,
  reference,
  assay = NULL,
  prefix = "ATAC",
  reference_assay = NULL,
  reference_reduction = "pca",
  ref_umap = NULL,
  reference_dims = 1:30,
  reference_label = NULL,
  label_method = c("Seurat", "scOMM"),
  add_gene_activity = TRUE,
  gene_activity_assay = "ACTIVITY",
  dims = 2:30,
  features = NULL,
  reduction_project_method = "pcaproject",
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  k.weight = 100,
  projection_method = c("model", "knn"),
  nn_method = NULL,
  k = 30,
  distance_metric = "cosine",
  vote_fun = "mean",
  evaluate = FALSE,
  truth_col = NULL,
  tool_name = NULL,
  rare_threshold = 0.05,
  scomm_python = NULL,
  scomm_hidden_nodes = c(128, 64),
  scomm_epochs = 10,
  scomm_batch_size = 32,
  scomm_threshold = 0.5,
  scomm_seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- reference:

  RNA reference `Seurat` object used for mapping.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- prefix:

  Prefix used to resolve ATAC reductions. Default is `"ATAC"`.

- reference_assay:

  Assay used in the reference object.

- reference_reduction:

  Reduction used in the reference object.

- ref_umap:

  UMAP reduction in the RNA reference used for query projection.

- reference_dims:

  Dimensions used from the reference reduction.

- reference_label:

  Metadata column in the reference used as transfer labels.

- label_method:

  Label-transfer backend used after anchor mapping. One of `"Seurat"` or
  `"scOMM"`.

- add_gene_activity:

  Whether to calculate a gene activity assay for the query.

- gene_activity_assay:

  Name of the gene activity assay used for mapping.

- dims:

  Query reduction dimensions used by `TransferData`.

- features:

  Features used by `FindTransferAnchors`. If `NULL`, reference variable
  features are used.

- reduction_project_method:

  Anchor projection method. Default is `"pcaproject"` for RNA reference
  mapping.

- k.anchor, k.filter, k.score, k.weight:

  Parameters passed to Seurat anchor finding/integration.

- projection_method, nn_method, k, distance_metric, vote_fun:

  Passed to
  [RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md)
  for UMAP projection and neighbor voting.

- evaluate:

  Whether to compute mapping metrics against a truth label.

- truth_col:

  Metadata column in `srt` used as the truth label when
  `evaluate = TRUE`.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- rare_threshold:

  Maximum class proportion used to define rare classes when calculating
  `rare_recall`.

- scomm_python:

  Optional Python binary used by the `scOMM` backend. If `NULL`,
  `SCOP_SCOMM_PYTHON` is consulted and reticulate defaults are used
  otherwise.

- scomm_hidden_nodes, scomm_epochs, scomm_batch_size, scomm_threshold,
  scomm_seed:

  Parameters passed to the optional `scOMM` backend.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A query `Seurat` object mapped into the RNA reference space.

## Examples

``` r
data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub <- standard_scop(
  pbmcmultiome_sub,
  assay = "RNA",
  linear_reduction_dims = 20
)
reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:250])
query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[251:350])
query <- standard_scop(
  query,
  assay = "peaks",
  normalization_method = "TFIDF",
  linear_reduction_dims = 20
)
query <- RunReferenceMapping(
  srt = query,
  reference = reference,
  assay = "peaks",
  reference_assay = "RNA",
  reference_reduction = "Standardpca",
  ref_umap = "StandardUMAP2D",
  reference_label = "CellType",
  reference_dims = 1:10,
  dims = 2:10
)
```
