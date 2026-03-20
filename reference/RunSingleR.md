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
#> ℹ [2026-03-20 09:41:33] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-03-20 09:41:34] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-03-20 09:41:34] Checking a list of <Seurat>...
#> ! [2026-03-20 09:41:34] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:41:34] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-03-20 09:41:36] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-03-20 09:41:36] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:41:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-03-20 09:41:38] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-03-20 09:41:39] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:41:39] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-03-20 09:41:40] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-03-20 09:41:41] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:41:41] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-03-20 09:41:42] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-03-20 09:41:43] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:41:43] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-03-20 09:41:45] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-03-20 09:41:45] Use the separate HVF from `srt_list`
#> ℹ [2026-03-20 09:41:45] Number of available HVF: 2000
#> ℹ [2026-03-20 09:41:46] Finished check

# Annotation
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-03-20 09:41:48] Start standard scop workflow...
#> ℹ [2026-03-20 09:41:49] Checking a list of <Seurat>...
#> ! [2026-03-20 09:41:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-20 09:41:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:41:51] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-20 09:41:52] Use the separate HVF from `srt_list`
#> ℹ [2026-03-20 09:41:52] Number of available HVF: 2000
#> ℹ [2026-03-20 09:41:52] Finished check
#> ℹ [2026-03-20 09:41:52] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-20 09:41:53] Perform pca linear dimension reduction
#> ℹ [2026-03-20 09:41:54] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-20 09:41:54] Reorder clusters...
#> ℹ [2026-03-20 09:41:54] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-20 09:41:54] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-20 09:41:58] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-20 09:42:03] Run scop standard workflow completed
pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "Standardpca_SNN_res.0.6",
  ref_group = "celltype"
)
#> ℹ [2026-03-20 09:42:03] Start SingleR annotation
#> ℹ [2026-03-20 09:48:38] Data type is log-normalized
#> ℹ [2026-03-20 09:48:38] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-20 09:48:39] Data type is log-normalized
#> ℹ [2026-03-20 09:48:39] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-20 09:48:43] Perform "SingleRCluster"
#> ✔ [2026-03-20 09:48:44] SingleR annotation completed
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
#> ℹ [2026-03-20 09:48:45] Start SingleR annotation
#> ℹ [2026-03-20 09:48:45] Data type is log-normalized
#> ℹ [2026-03-20 09:48:45] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-20 09:48:48] Data type is log-normalized
#> ℹ [2026-03-20 09:48:48] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-20 09:48:53] Perform "SingleRCell"
#> ✔ [2026-03-20 09:48:57] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "CellType")
)
```
