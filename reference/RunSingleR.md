# Annotate single cells using SingleR

Annotate single cells using SingleR

## Usage

``` r
RunSingleR(
  srt_query,
  srt_ref,
  query_group = NULL,
  ref_group = NULL,
  query_assay = "RNA",
  ref_assay = "RNA",
  genes = "de",
  de.method = "wilcox",
  sd.thresh = 1,
  de.n = NULL,
  aggr.ref = FALSE,
  aggr.args = list(),
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  prune = TRUE,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- query_group:

  A character vector specifying the column name in the `srt_query`
  metadata that represents the cell grouping.

- ref_group:

  A character vector specifying the column name in the `srt_ref`
  metadata that represents the cell grouping.

- query_assay:

  A character vector specifying the assay to be used for the query data.
  Default is the default assay of the `srt_query` object.

- ref_assay:

  A character vector specifying the assay to be used for the reference
  data. Default is the default assay of the `srt_ref` object.

- genes:

  A string containing `"de"`, indicating that markers should be
  calculated from `ref`. For back compatibility, other string values are
  allowed but will be ignored with a deprecation warning.

  Alternatively, if `ref` is *not* a list, `genes` can be either:

  - A list of lists of character vectors containing DE genes between
    pairs of labels.

  - A list of character vectors containing marker genes for each label.

  If `ref` *is* a list, `genes` can be a list of length equal to `ref`.
  Each element of the list should be one of the two above choices
  described for non-list `ref`, containing markers for labels in the
  corresponding entry of `ref`.

- de.method:

  String specifying how DE genes should be detected between pairs of
  labels. Defaults to `"classic"`, which sorts genes by the log-fold
  changes and takes the top `de.n`. Other options are `"wilcox"` and
  `"t"`, see Details. Ignored if `genes` is a list of markers/DE genes.

- sd.thresh:

  Deprecated and ignored.

- de.n:

  An integer scalar specifying the number of DE genes to use when
  `genes="de"`. If `de.method="classic"`, defaults to
  `500 * (2/3) ^ log2(N)` where `N` is the number of unique labels.
  Otherwise, defaults to 10. Ignored if `genes` is a list of markers/DE
  genes.

- aggr.ref, aggr.args:

  Arguments controlling the aggregation of the references prior to
  annotation, see
  [`trainSingleR`](https://rdrr.io/pkg/SingleR/man/trainSingleR.html).

- quantile:

  "quantile" parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- fine.tune:

  `"fine.tune"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- tune.thresh:

  `"tune.thresh"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- prune:

  `"prune"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

An annotate `Seurat` object. The annotation results are stored in the
`singler_annotation` column of the meta data, and the corresponding
scores are stored in the `singler_score` column.

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md)

## Examples

``` r
data(panc8_sub)
# Simply convert genes from human to mouse and preprocess the data
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
#> ℹ [2026-04-22 09:04:37] Rename features for the assay: RNA
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-04-22 09:04:37] Start standard processing workflow...
#> ℹ [2026-04-22 09:04:37] Checking a list of <Seurat>...
#> ! [2026-04-22 09:04:37] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-22 09:04:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-22 09:04:40] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-22 09:04:40] Use the separate HVF from `srt_list`
#> ℹ [2026-04-22 09:04:41] Number of available HVF: 2000
#> ℹ [2026-04-22 09:04:41] Finished check
#> ℹ [2026-04-22 09:04:41] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-22 09:04:42] Perform pca linear dimension reduction
#> ℹ [2026-04-22 09:04:43] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-22 09:04:43] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-22 09:04:43] Reorder clusters...
#> ℹ [2026-04-22 09:04:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-22 09:04:44] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-22 09:04:44] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-22 09:04:49] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-22 09:04:54] Standard processing workflow completed

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-22 09:04:54] Start standard processing workflow...
#> ℹ [2026-04-22 09:04:55] Checking a list of <Seurat>...
#> ! [2026-04-22 09:04:55] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-22 09:04:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-22 09:04:57] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-22 09:04:58] Use the separate HVF from `srt_list`
#> ℹ [2026-04-22 09:04:58] Number of available HVF: 2000
#> ℹ [2026-04-22 09:04:58] Finished check
#> ℹ [2026-04-22 09:04:58] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-22 09:04:58] Perform pca linear dimension reduction
#> ℹ [2026-04-22 09:04:59] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-22 09:05:00] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-22 09:05:00] Reorder clusters...
#> ℹ [2026-04-22 09:05:00] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-22 09:05:00] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-22 09:05:00] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-22 09:05:05] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-22 09:05:10] Standard processing workflow completed
pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "Standardpca_SNN_res.0.6",
  ref_group = "celltype"
)
#> ℹ [2026-04-22 09:05:10] Start SingleR annotation
#> ℹ [2026-04-22 09:11:55] Data type is log-normalized
#> ℹ [2026-04-22 09:11:55] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-04-22 09:11:56] Data type is log-normalized
#> ℹ [2026-04-22 09:11:56] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-04-22 09:11:56] Perform "SingleRCluster"
#> ✔ [2026-04-22 09:11:57] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "SubCellType")
)


pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = NULL,
  ref_group = "celltype"
)
#> ℹ [2026-04-22 09:11:58] Start SingleR annotation
#> ℹ [2026-04-22 09:11:58] Data type is log-normalized
#> ℹ [2026-04-22 09:11:58] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-04-22 09:11:59] Data type is log-normalized
#> ℹ [2026-04-22 09:11:59] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-04-22 09:12:00] Perform "SingleRCell"
#> ✔ [2026-04-22 09:12:04] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "SubCellType"),
  label = TRUE
)


CellCorHeatmap(
  pancreas_sub,
  group.by = "singler_annotation",
  assay = "RNA",
  layer = "data",
  method = "spearman"
)
#> Error in CellCorHeatmap(pancreas_sub, group.by = "singler_annotation",     assay = "RNA", layer = "data", method = "spearman"): unused arguments (group.by = "singler_annotation", assay = "RNA", layer = "data", method = "spearman")

ht1 <- CellCorHeatmap(
  srt_query = pancreas_sub,
  srt_ref = pancreas_sub,
  query_group = "SubCellType",
  cluster_rows = TRUE,
  ref_group = "singler_annotation",
  cluster_columns = TRUE,
  width = 2,
  height = 2
)
#> ℹ [2026-04-22 09:12:05] Drop [1] 19 cells with NA in the ref_group
#> ℹ [2026-04-22 09:12:05] Use the HVF to calculate distance metric
#> ℹ [2026-04-22 09:12:05] Use [1] 2000 features to calculate distance.
#> ℹ [2026-04-22 09:12:05] Detected query data type: "log_normalized_counts"
#> ℹ [2026-04-22 09:12:05] Detected reference data type: "log_normalized_counts"
#> ℹ [2026-04-22 09:12:05] Calculate similarity...
#> ℹ [2026-04-22 09:12:05] Use raw method to find neighbors
#> ℹ [2026-04-22 09:12:05] Predict cell type...

ht1$plot
```
