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

  `"genes"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- de.method:

  `"de.method"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

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
#> ℹ [2026-04-03 04:35:34] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-04-03 04:35:34] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-03 04:35:34] Checking a list of <Seurat>...
#> ! [2026-04-03 04:35:35] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:35:35] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-04-03 04:35:36] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-04-03 04:35:37] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:35:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-04-03 04:35:39] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-04-03 04:35:39] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:35:39] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-04-03 04:35:41] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-04-03 04:35:41] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:35:41] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-04-03 04:35:43] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-04-03 04:35:43] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:35:43] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-04-03 04:35:45] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-03 04:35:46] Use the separate HVF from `srt_list`
#> ℹ [2026-04-03 04:35:46] Number of available HVF: 2000
#> ℹ [2026-04-03 04:35:47] Finished check

# Annotation
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-03 04:35:50] Start standard processing workflow...
#> ℹ [2026-04-03 04:35:51] Checking a list of <Seurat>...
#> ! [2026-04-03 04:35:51] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 04:35:51] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 04:35:53] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 04:35:54] Use the separate HVF from `srt_list`
#> ℹ [2026-04-03 04:35:54] Number of available HVF: 2000
#> ℹ [2026-04-03 04:35:54] Finished check
#> ℹ [2026-04-03 04:35:54] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-03 04:35:55] Perform pca linear dimension reduction
#> ℹ [2026-04-03 04:35:55] Use stored estimated dimensions 1:12 for Standardpca
#> ℹ [2026-04-03 04:35:56] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-03 04:35:56] Reorder clusters...
#> ℹ [2026-04-03 04:35:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-03 04:35:56] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-03 04:35:56] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ℹ [2026-04-03 04:36:01] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ✔ [2026-04-03 04:36:06] Standard processing workflow completed
pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "Standardpca_SNN_res.0.6",
  ref_group = "celltype"
)
#> ℹ [2026-04-03 04:36:06] Start SingleR annotation
#> ℹ [2026-04-03 04:42:43] Data type is log-normalized
#> ℹ [2026-04-03 04:42:43] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-04-03 04:42:46] Data type is log-normalized
#> ℹ [2026-04-03 04:42:46] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-04-03 04:42:51] Perform "SingleRCluster"
#> ✔ [2026-04-03 04:42:52] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "CellType")
)


pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = NULL,
  ref_group = "celltype"
)
#> ℹ [2026-04-03 04:42:52] Start SingleR annotation
#> ℹ [2026-04-03 04:42:53] Data type is log-normalized
#> ℹ [2026-04-03 04:42:53] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-04-03 04:42:55] Data type is log-normalized
#> ℹ [2026-04-03 04:42:55] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-04-03 04:43:00] Perform "SingleRCell"
#> ✔ [2026-04-03 04:43:04] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "CellType")
)
```
