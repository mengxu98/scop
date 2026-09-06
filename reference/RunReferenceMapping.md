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

  A `Seurat` object.

- reference:

  RNA reference `Seurat` object used for mapping.

- assay:

  Assay to use. `NULL` uses the default assay.

- prefix:

  Prefix used to resolve ATAC reductions.

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
pbmcmultiome_sub <- pbmcmultiome_sub[, seq_len(200)]
pbmcmultiome_sub <- RunStandardWorkflow(
  pbmcmultiome_sub,
  assay = "RNA",
  linear_reduction_dims = 10,
  verbose = FALSE
)
#> ℹ [2026-09-06 22:31:35] Skip `log1p()` because `layer = data` is not "counts"
reference <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[1:120])
query <- subset(pbmcmultiome_sub, cells = colnames(pbmcmultiome_sub)[121:200])
query <- RunStandardWorkflow(
  query,
  assay = "peaks",
  normalization_method = "TFIDF",
  linear_reduction_dims = 10,
  verbose = FALSE
)
#> ! [2026-09-06 22:31:43] Only one cluster found
query <- RunReferenceMapping(
  srt = query,
  reference = reference,
  assay = "peaks",
  reference_assay = "RNA",
  reference_reduction = "Standardpca",
  ref_umap = "StandardUMAP2D",
  reference_label = "CellType",
  reference_dims = 1:10,
  dims = 2:5,
  verbose = FALSE
)
#> Requested to reuse weights matrix, but no weights found. Computing new weights.
#> 
#> Integrating dataset 2 with reference dataset
#> Finding integration vectors
#> Finding integration vector weights
#> Integrating data
#> ℹ [2026-09-06 22:31:54] No UMAP model detected. Set the `projection_method` to "knn"
#> ℹ [2026-09-06 22:31:54] Use the reduction to calculate distance metric
#> ℹ [2026-09-06 22:31:54] Use raw method to find neighbors
#> ℹ [2026-09-06 22:31:54] Predicting cell types based on ref_group
head(query[[]][, grep("^predicted_", colnames(query[[]]), value = TRUE), drop = FALSE])
#>                    predicted_predicted.id predicted_prediction.score.NK
#> CAAGGCCTCTAGCGTG-1                      B                    0.03359597
#> CAAGGCTGTGTTTGCT-1                  CD4 T                    0.09557635
#> CAAGGGAGTGGTTCCC-1                      B                    0.02860680
#> CAAGGTAAGCCACATG-1                  CD4 T                    0.02570682
#> CAAGTAACAACAGGAT-1                   Mono                    0.00000000
#> CAAGTAACACAACAGG-1                  CD4 T                    0.04576571
#>                    predicted_prediction.score.CD8.T
#> CAAGGCCTCTAGCGTG-1                        0.0000000
#> CAAGGCTGTGTTTGCT-1                        0.3184703
#> CAAGGGAGTGGTTCCC-1                        0.0000000
#> CAAGGTAAGCCACATG-1                        0.3904068
#> CAAGTAACAACAGGAT-1                        0.0000000
#> CAAGTAACACAACAGG-1                        0.4038763
#>                    predicted_prediction.score.Mono predicted_prediction.score.B
#> CAAGGCCTCTAGCGTG-1                      0.00000000                   0.80292686
#> CAAGGCTGTGTTTGCT-1                      0.01509085                   0.00000000
#> CAAGGGAGTGGTTCCC-1                      0.00000000                   0.81837625
#> CAAGGTAAGCCACATG-1                      0.00000000                   0.00000000
#> CAAGTAACAACAGGAT-1                      0.56249814                   0.02870852
#> CAAGTAACACAACAGG-1                      0.00000000                   0.00000000
#>                    predicted_prediction.score.CD4.T
#> CAAGGCCTCTAGCGTG-1                        0.0000000
#> CAAGGCTGTGTTTGCT-1                        0.5394940
#> CAAGGGAGTGGTTCCC-1                        0.0000000
#> CAAGGTAAGCCACATG-1                        0.5838864
#> CAAGTAACAACAGGAT-1                        0.0000000
#> CAAGTAACACAACAGG-1                        0.5503580
#>                    predicted_prediction.score.DC predicted_prediction.score.max
#> CAAGGCCTCTAGCGTG-1                    0.16347717                      0.8029269
#> CAAGGCTGTGTTTGCT-1                    0.03136851                      0.5394940
#> CAAGGGAGTGGTTCCC-1                    0.15301695                      0.8183762
#> CAAGGTAAGCCACATG-1                    0.00000000                      0.5838864
#> CAAGTAACAACAGGAT-1                    0.40879333                      0.5624981
#> CAAGTAACACAACAGG-1                    0.00000000                      0.5503580
#>                    predicted_ref_group
#> CAAGGCCTCTAGCGTG-1                Mono
#> CAAGGCTGTGTTTGCT-1                   B
#> CAAGGGAGTGGTTCCC-1                Mono
#> CAAGGTAAGCCACATG-1                   B
#> CAAGTAACAACAGGAT-1                  DC
#> CAAGTAACACAACAGG-1                   B
```
