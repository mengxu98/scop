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

  A `Seurat` object.

- reference:

  RNA reference `Seurat` object used for label transfer.

- assay:

  Assay to use. `NULL` uses the default assay.

- method:

  Label-transfer backend. One of `"Seurat"` or `"scOMM"`.

- prefix:

  Prefix used to resolve ATAC reductions.

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
pbmcmultiome_sub <- pbmcmultiome_sub[, seq_len(200)]
pbmcmultiome_sub <- RunStandardWorkflow(
  pbmcmultiome_sub,
  assay = "RNA",
  linear_reduction_dims = 10,
  verbose = FALSE
)
#> ℹ [2026-09-06 22:27:43] Skip `log1p()` because `layer = data` is not "counts"
reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:120])
query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[121:200])
query <- RunStandardWorkflow(
  query,
  assay = "peaks",
  normalization_method = "TFIDF",
  linear_reduction_dims = 10,
  verbose = FALSE
)
#> ! [2026-09-06 22:27:51] Only one cluster found
query <- RunLabelTransfer(
  srt = query,
  reference = reference,
  assay = "peaks",
  reference_assay = "RNA",
  reference_reduction = "Standardpca",
  reference_label = "CellType",
  reference_dims = 1:10,
  dims = 2:5,
  verbose = FALSE
)
head(query[[]][, grep("^predicted_", colnames(query[[]]), value = TRUE), drop = FALSE])
#>                    predicted_predicted.id predicted_prediction.score.NK
#> CAAGGCCTCTAGCGTG-1                      B                             0
#> CAAGGCTGTGTTTGCT-1                  CD4 T                             0
#> CAAGGGAGTGGTTCCC-1                      B                             0
#> CAAGGTAAGCCACATG-1                  CD4 T                             0
#> CAAGTAACAACAGGAT-1                   Mono                             0
#> CAAGTAACACAACAGG-1                  CD4 T                             0
#>                    predicted_prediction.score.CD8.T
#> CAAGGCCTCTAGCGTG-1                       0.03249353
#> CAAGGCTGTGTTTGCT-1                       0.41327998
#> CAAGGGAGTGGTTCCC-1                       0.02760854
#> CAAGGTAAGCCACATG-1                       0.43976034
#> CAAGTAACAACAGGAT-1                       0.00000000
#> CAAGTAACACAACAGG-1                       0.41342177
#>                    predicted_prediction.score.Mono predicted_prediction.score.B
#> CAAGGCCTCTAGCGTG-1                      0.00000000                    0.8678555
#> CAAGGCTGTGTTTGCT-1                      0.05208446                    0.0000000
#> CAAGGGAGTGGTTCCC-1                      0.00000000                    0.8809738
#> CAAGGTAAGCCACATG-1                      0.00000000                    0.0000000
#> CAAGTAACAACAGGAT-1                      0.55003767                    0.0000000
#> CAAGTAACACAACAGG-1                      0.00000000                    0.0000000
#>                    predicted_prediction.score.CD4.T
#> CAAGGCCTCTAGCGTG-1                       0.00000000
#> CAAGGCTGTGTTTGCT-1                       0.53424081
#> CAAGGGAGTGGTTCCC-1                       0.01323434
#> CAAGGTAAGCCACATG-1                       0.56023966
#> CAAGTAACAACAGGAT-1                       0.00000000
#> CAAGTAACACAACAGG-1                       0.58657823
#>                    predicted_prediction.score.DC predicted_prediction.score.max
#> CAAGGCCTCTAGCGTG-1                  0.0996509484                      0.8678555
#> CAAGGCTGTGTTTGCT-1                  0.0003947406                      0.5342408
#> CAAGGGAGTGGTTCCC-1                  0.0781832706                      0.8809738
#> CAAGGTAAGCCACATG-1                  0.0000000000                      0.5602397
#> CAAGTAACAACAGGAT-1                  0.4499623257                      0.5500377
#> CAAGTAACACAACAGG-1                  0.0000000000                      0.5865782
```
