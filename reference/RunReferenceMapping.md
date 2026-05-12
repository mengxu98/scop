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
#> ℹ [2026-05-12 05:12:43] Start standard processing workflow...
#> ℹ [2026-05-12 05:12:44] Checking a list of <Seurat>...
#> ! [2026-05-12 05:12:44] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-12 05:12:44] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 05:12:46] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> Warning: pseudoinverse used at -2.3979
#> Warning: neighborhood radius 0.30103
#> Warning: reciprocal condition number  9.9917e-16
#> ℹ [2026-05-12 05:12:46] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 05:12:47] Number of available HVF: 2000
#> ℹ [2026-05-12 05:12:47] Finished check
#> ℹ [2026-05-12 05:12:47] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-12 05:12:47] Perform pca linear dimension reduction
#> ℹ [2026-05-12 05:12:47] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-12 05:12:48] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-12 05:12:48] Reorder clusters...
#> ℹ [2026-05-12 05:12:48] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-12 05:12:48] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-12 05:12:48] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-12 05:12:53] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-12 05:12:57] Standard processing workflow completed
reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:250])
query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[251:350])
query <- standard_scop(
  query,
  assay = "peaks",
  normalization_method = "TFIDF",
  linear_reduction_dims = 20
)
#> ℹ [2026-05-12 05:12:58] Start standard processing workflow...
#> ℹ [2026-05-12 05:12:58] Checking a list of <Seurat>...
#> ! [2026-05-12 05:12:58] Data 1/1 of the `srt_list` is "raw_counts"
#> ℹ [2026-05-12 05:12:58] Perform `RunTFIDF()` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 05:12:58] Perform `FindTopFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-12 05:12:58] Use the separate HVF from `srt_list`
#> ℹ [2026-05-12 05:12:58] Number of available HVF: 11426
#> ℹ [2026-05-12 05:12:58] Finished check
#> ℹ [2026-05-12 05:12:58] `normalization_method` is TFIDF. Use lsi workflow
#> ℹ [2026-05-12 05:12:58] Perform svd linear dimension reduction
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-05-12 05:12:59] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-12 05:12:59] Reorder clusters...
#> ℹ [2026-05-12 05:12:59] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-12 05:12:59] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-12 05:12:59] Perform umap nonlinear dimension reduction using ATACsvd (2:30)
#> ℹ [2026-05-12 05:13:03] Perform umap nonlinear dimension reduction using ATACsvd (2:30)
#> ✔ [2026-05-12 05:13:08] Standard processing workflow completed
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
#> ℹ [2026-05-12 05:13:08] Calculating gene activity assay...
#> Extracting gene coordinates
#> Error in atac_add_activity(srt = srt, assay = assay, gene_activity_assay = gene_activity_assay,     verbose = verbose): Unable to calculate gene activity assay: "No fragment information found
#> for requested assay". Please provide a valid ATAC annotation.
```
