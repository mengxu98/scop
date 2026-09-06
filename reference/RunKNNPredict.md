# Run KNN prediction

Performs KNN prediction to annotate cell types based on reference
scRNA-seq or bulk RNA-seq data.

## Usage

``` r
RunKNNPredict(
  srt_query,
  srt_ref = NULL,
  bulk_ref = NULL,
  query_group = NULL,
  ref_group = NULL,
  query_assay = NULL,
  ref_assay = NULL,
  query_reduction = NULL,
  ref_reduction = NULL,
  query_dims = 1:30,
  ref_dims = 1:30,
  query_collapsing = !is.null(query_group),
  ref_collapsing = TRUE,
  return_full_distance_matrix = FALSE,
  features = NULL,
  features_type = c("HVF", "DE"),
  feature_source = "both",
  nfeatures = 2000,
  DEtest_param = list(max.cells.per.ident = 200, test.use = "wilcox"),
  DE_threshold = "p_val_adj < 0.05",
  nn_method = NULL,
  distance_metric = "cosine",
  k = 30,
  filter_lowfreq = 0,
  prefix = "KNNPredict",
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- bulk_ref:

  A cell atlas matrix, where cell types are represented by columns and
  genes are represented by rows. Either `srt_ref` or `bulk_ref` must be
  provided.

- query_group:

  Column name in the `srt_query` metadata that represents the cell
  grouping.

- ref_group:

  Column name in the `srt_ref` metadata that represents the cell
  grouping.

- query_assay:

  Assay to be used for the query data. Default is the default assay of
  the `srt_query` object.

- ref_assay:

  Assay to be used for the reference data. Default is the default assay
  of the `srt_ref` object.

- query_reduction:

  Dimensionality reduction method used for the query data. If NULL, the
  function will use the default reduction method specified in the
  `srt_query` object.

- ref_reduction:

  Dimensionality reduction method used for the reference data. If NULL,
  the function will use the default reduction method specified in the
  `srt_ref` object.

- query_dims:

  Dimensions to be used for the query data. Default is the first `30`
  dimensions.

- ref_dims:

  Dimensions to be used for the reference data. Default is the first
  `30` dimensions.

- query_collapsing:

  A boolean value indicating whether the query data should be collapsed
  to group-level average expression values. If TRUE, the function will
  calculate the average expression values for each group in the query
  data and the annotation will be performed separately for each group.
  Otherwise it will use the raw expression values for each cell.

- ref_collapsing:

  A boolean value indicating whether the reference data should be
  collapsed to group-level average expression values. If TRUE, the
  function will calculate the average expression values for each group
  in the reference data and the annotation will be performed separately
  for each group. Otherwise it will use the raw expression values for
  each cell.

- return_full_distance_matrix:

  A boolean value indicating whether the full distance matrix should be
  returned. If TRUE, the function will return the distance matrix used
  for the KNN prediction, otherwise it will only return the annotated
  cell types.

- features:

  Features to be used for the KNN prediction. If `NULL`, all the
  features in the query and reference data will be used.

- features_type:

  Type of features to be used for the KNN prediction. Must be one of
  "HVF" (highly variable features) or "DE" (differentially expressed
  features).

- feature_source:

  The source of the features to be used. Must be one of "both", "query",
  or "ref".

- nfeatures:

  Maximum number of features to be used for the KNN prediction.

- DEtest_param:

  A list of parameters to be passed to the differential expression test
  function if `features_type` is set to "DE". Default is
  `list(max.cells.per.ident = 200, test.use = "wilcox")`.

- DE_threshold:

  Threshold used to filter the DE features. If using "roc" test,
  `DE_threshold` should be needs to be reassigned. e.g. "power \> 0.5".

- nn_method:

  Nearest neighbor search method to use. Options are "raw", "annoy", and
  "rann". If "raw" is selected, the function will use the brute-force
  method to find the nearest neighbors. If "annoy" is selected, the
  function will use the Annoy library for approximate nearest neighbor
  search. If "rann" is selected, the function will use the RANN library
  for approximate nearest neighbor search. If not provided, the function
  will choose the search method based on the size of the query and
  reference datasets. For finite dense inputs using cosine or Euclidean
  distance, the raw path uses a bounded-memory compiled top-k kernel
  unless the full distance matrix is requested. Other metrics and sparse
  inputs retain the existing proxyC path.

- distance_metric:

  Distance metric to be used for calculating similarity between cells.
  Must be one of "cosine", "euclidean", "manhattan", or "hamming".

- k:

  A number of nearest neighbors to be considered for the KNN prediction.

- filter_lowfreq:

  Threshold for filtering low-frequency cell types from the predicted
  results. Cell types with a frequency lower than `filter_lowfreq` will
  be labelled as "unreliable". Default is `0`, which means no filtering
  will be performed.

- prefix:

  Prefix to be added to the resulting annotations.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md),
[RunSingleR](https://mengxu98.github.io/scop/reference/RunSingleR.md),
[CellCorHeatmap](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)

## Examples

``` r
data(panc8_sub)
keep <- panc8_sub$celltype %in% c("ductal", "alpha", "beta")
reference <- Seurat::NormalizeData(
  panc8_sub[, keep & panc8_sub$tech != "fluidigmc1"],
  verbose = FALSE
)
query <- Seurat::NormalizeData(
  panc8_sub[, keep & panc8_sub$tech == "fluidigmc1"],
  verbose = FALSE
)
genes <- head(intersect(rownames(reference), rownames(query)), 200)
query <- RunKNNPredict(
  query, srt_ref = reference, ref_group = "celltype",
  features = genes, nfeatures = length(genes),
  ref_collapsing = FALSE, k = 10, verbose = FALSE
)
query <- RunStandardWorkflow(query, verbose = FALSE, linear_reduction_dims = 10)
#> ℹ [2026-09-06 22:27:33] Skip `log1p()` because `layer = data` is not "counts"
CellDimPlot(
  query, group.by = "KNNPredict_classification",
  label = TRUE, legend.position = "bottom"
)

FeatureStatPlot(
  query,
  stat.by = "KNNPredict_prob",
  group.by = "KNNPredict_classification",
  plot_type = "box"
)
```
