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
#> ℹ [2026-06-29 04:34:30] Start standard processing workflow...
#> ℹ [2026-06-29 04:34:31] Checking a list of <Seurat>...
#> ! [2026-06-29 04:34:31] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-29 04:34:31] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-29 04:34:31] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> Warning: pseudoinverse used at -2.3979
#> Warning: neighborhood radius 0.30103
#> Warning: reciprocal condition number  1.2589e-15
#> ℹ [2026-06-29 04:34:31] Use the separate HVF from `srt_list`
#> ℹ [2026-06-29 04:34:31] Number of available HVF: 2000
#> ℹ [2026-06-29 04:34:31] Finished check
#> ℹ [2026-06-29 04:34:31] Perform `ScaleData()`
#> ℹ [2026-06-29 04:34:31] Perform pca linear dimension reduction
#> ℹ [2026-06-29 04:34:32] Use stored estimated dimensions 1:9 for Standardpca
#> ℹ [2026-06-29 04:34:32] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-29 04:34:32] Reorder clusters...
#> ℹ [2026-06-29 04:34:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-29 04:34:32] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-29 04:34:39] Standard processing workflow completed
reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:250])
query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[251:350])
query <- standard_scop(
  query,
  assay = "peaks",
  normalization_method = "TFIDF",
  linear_reduction_dims = 20
)
#> ℹ [2026-06-29 04:34:39] Start standard processing workflow...
#> ℹ [2026-06-29 04:34:39] Checking a list of <Seurat>...
#> ! [2026-06-29 04:34:40] Data 1/1 of the `srt_list` is "raw_counts"
#> ℹ [2026-06-29 04:34:40] Perform `RunTFIDF()` on 1/1 of `srt_list`...
#> ℹ [2026-06-29 04:34:40] Perform `FindTopFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-29 04:34:40] Use the separate HVF from `srt_list`
#> ℹ [2026-06-29 04:34:40] Number of available HVF: 11426
#> ℹ [2026-06-29 04:34:40] Finished check
#> ℹ [2026-06-29 04:34:40] `normalization_method` is TFIDF. Use lsi workflow
#> ℹ [2026-06-29 04:34:40] Perform svd linear dimension reduction
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-06-29 04:34:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-29 04:34:40] Reorder clusters...
#> ℹ [2026-06-29 04:34:40] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-29 04:34:40] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-29 04:34:47] Standard processing workflow completed
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
#> ℹ [2026-06-29 04:34:47] Use existing query assay "RNA" as `gene_activity_assay`
#> ℹ [2026-06-29 04:34:47] Adjust `k.filter` from 200 to 99 for small-sample ATAC mapping
#> ℹ [2026-06-29 04:34:47] Finding RNA-to-ATAC anchors for query mapping...
#> ℹ [2026-06-29 04:34:50] Adjust `k.weight` from 100 to 95 for small-sample ATAC mapping
#> Warning: Max dims.to.integrate is larger than the max dims for at least one of the reductions specified. Setting dims.to.integrate to 2,3,4,5,6,7,8,9 and continuing.
#> Requested to reuse weights matrix, but no weights found. Computing new weights.
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> Warning: Layer counts isn't present in the assay object; returning NULL
#> 
#> Integrating dataset 2 with reference dataset
#> Finding integration vectors
#> Finding integration vector weights
#> Integrating data
#> ℹ [2026-06-29 04:34:51] Adjust `k.filter` from 200 to 99 for small-sample ATAC mapping
#> ℹ [2026-06-29 04:34:51] Running RNA reference label transfer for ATAC cells...
#> ℹ [2026-06-29 04:34:54] Adjust `k.weight` from 95 to 94 for small-sample ATAC mapping
#> ℹ [2026-06-29 04:34:54] No UMAP model detected. Set the `projection_method` to "knn"
#> ℹ [2026-06-29 04:34:54] Use the reduction to calculate distance metric
#> ℹ [2026-06-29 04:34:54] Use raw method to find neighbors
#> ℹ [2026-06-29 04:34:54] Predicting cell types based on ref_group
```
