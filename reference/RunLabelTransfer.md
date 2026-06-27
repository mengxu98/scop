# Transfer reference labels to query cells

Standalone label-transfer workflow for query cells using a reference
object. The current implementation is optimized for scATAC query objects
mapped to a scRNA-seq reference via gene activity.

## Usage

``` r
RunLabelTransfer(
  srt,
  reference,
  assay = NULL,
  method = c("Seurat", "scOMM"),
  prefix = "ATAC",
  reference_assay = NULL,
  reference_reduction = "pca",
  reference_dims = 1:30,
  reference_label = NULL,
  add_gene_activity = TRUE,
  gene_activity_assay = "ACTIVITY",
  weight_reduction = NULL,
  dims = 2:30,
  features = NULL,
  prediction_prefix = NULL,
  k.weight = 100,
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

  RNA reference `Seurat` object used for label transfer.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- method:

  Label-transfer backend. One of `"Seurat"` or `"scOMM"`.

- prefix:

  Prefix used to resolve ATAC reductions. Default is `"ATAC"`.

- reference_assay:

  Assay used in the reference object.

- reference_reduction:

  Reduction used in the reference object.

- reference_dims:

  Dimensions used from the reference reduction.

- reference_label:

  Metadata column in the reference used as transfer labels.

- add_gene_activity:

  Whether to calculate a gene activity assay for the query.

- gene_activity_assay:

  Name of the gene activity assay used for mapping.

- weight_reduction:

  Reduction in `srt` used to weight transferred labels. If `NULL`, an
  ATAC linear reduction is resolved automatically from
  `ATAC_default_linear_reduction`, `{prefix}lsi`, `{prefix}svd`, or the
  current default reduction.

- dims:

  Query reduction dimensions used by `TransferData`.

- features:

  Features used by `FindTransferAnchors`. If `NULL`, reference variable
  features are used.

- prediction_prefix:

  Prefix added to prediction metadata columns. If `NULL`, `"predicted_"`
  is used for `method = "Seurat"` and `"scomm_"` is used for
  `method = "scOMM"`.

- k.weight:

  Number of neighbors used when weighting transfer anchors.

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

A `Seurat` object with prediction metadata added.

## Examples

``` r
data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub <- standard_scop(
  pbmcmultiome_sub,
  assay = "RNA",
  linear_reduction_dims = 20
)
#> ℹ [2026-06-27 18:09:11] Start standard processing workflow...
#> ℹ [2026-06-27 18:09:12] Checking a list of <Seurat>...
#> ! [2026-06-27 18:09:12] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-27 18:09:12] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 18:09:12] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> Warning: pseudoinverse used at -2.3979
#> Warning: neighborhood radius 0.30103
#> Warning: reciprocal condition number  1.2589e-15
#> ℹ [2026-06-27 18:09:12] Use the separate HVF from `srt_list`
#> ℹ [2026-06-27 18:09:12] Number of available HVF: 2000
#> ℹ [2026-06-27 18:09:12] Finished check
#> ℹ [2026-06-27 18:09:12] Perform `ScaleData()`
#> ℹ [2026-06-27 18:09:12] Perform pca linear dimension reduction
#> ℹ [2026-06-27 18:09:13] Use stored estimated dimensions 1:9 for Standardpca
#> ℹ [2026-06-27 18:09:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-27 18:09:13] Reorder clusters...
#> ℹ [2026-06-27 18:09:13] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-27 18:09:13] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-27 18:09:19] Standard processing workflow completed
reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:250])
query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[251:350])
query <- standard_scop(
  query,
  assay = "peaks",
  normalization_method = "TFIDF",
  linear_reduction_dims = 20
)
#> ℹ [2026-06-27 18:09:20] Start standard processing workflow...
#> ℹ [2026-06-27 18:09:20] Checking a list of <Seurat>...
#> ! [2026-06-27 18:09:20] Data 1/1 of the `srt_list` is "raw_counts"
#> ℹ [2026-06-27 18:09:20] Perform `RunTFIDF()` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 18:09:20] Perform `FindTopFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 18:09:20] Use the separate HVF from `srt_list`
#> ℹ [2026-06-27 18:09:20] Number of available HVF: 11426
#> ℹ [2026-06-27 18:09:20] Finished check
#> ℹ [2026-06-27 18:09:20] `normalization_method` is TFIDF. Use lsi workflow
#> ℹ [2026-06-27 18:09:20] Perform svd linear dimension reduction
#> Running SVD
#> Scaling cell embeddings
#> ℹ [2026-06-27 18:09:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-27 18:09:21] Reorder clusters...
#> ℹ [2026-06-27 18:09:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-27 18:09:21] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-27 18:09:27] Standard processing workflow completed
query <- RunLabelTransfer(
  srt = query,
  reference = reference,
  assay = "peaks",
  reference_assay = "RNA",
  reference_reduction = "Standardpca",
  reference_label = "CellType",
  reference_dims = 1:10,
  dims = 2:10
)
#> ℹ [2026-06-27 18:09:27] Use existing query assay "RNA" as `gene_activity_assay`
#> ℹ [2026-06-27 18:09:27] Use "ATAClsi" as the ATAC weight reduction
#> ℹ [2026-06-27 18:09:27] Adjust `k.filter` from 200 to 99 for small-sample ATAC mapping
#> ℹ [2026-06-27 18:09:27] Running RNA reference label transfer for ATAC cells...
#> ℹ [2026-06-27 18:09:30] Adjust `k.weight` from 100 to 96 for small-sample ATAC mapping
```
