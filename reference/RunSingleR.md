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
#> ℹ [2026-03-08 08:06:53] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-03-08 08:06:53] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-03-08 08:06:53] Checking a list of <Seurat>...
#> ! [2026-03-08 08:06:54] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-08 08:06:54] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-03-08 08:06:55] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-03-08 08:06:56] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-08 08:06:56] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-03-08 08:06:58] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-03-08 08:06:58] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-08 08:06:58] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-03-08 08:07:00] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-03-08 08:07:00] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-08 08:07:00] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-03-08 08:07:02] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-03-08 08:07:02] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-03-08 08:07:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-03-08 08:07:04] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-03-08 08:07:04] Use the separate HVF from `srt_list`
#> ℹ [2026-03-08 08:07:05] Number of available HVF: 2000
#> ℹ [2026-03-08 08:07:05] Finished check

# Annotation
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-03-08 08:07:08] Start standard scop workflow...
#> ℹ [2026-03-08 08:07:08] Checking a list of <Seurat>...
#> ! [2026-03-08 08:07:08] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-08 08:07:08] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-08 08:07:11] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-08 08:07:11] Use the separate HVF from `srt_list`
#> ℹ [2026-03-08 08:07:11] Number of available HVF: 2000
#> ℹ [2026-03-08 08:07:11] Finished check
#> ℹ [2026-03-08 08:07:12] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-08 08:07:12] Perform pca linear dimension reduction
#> ℹ [2026-03-08 08:07:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-08 08:07:13] Reorder clusters...
#> ℹ [2026-03-08 08:07:13] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-08 08:07:13] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-08 08:07:18] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-08 08:07:23] Run scop standard workflow completed
pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "Standardpca_SNN_res.0.6",
  ref_group = "celltype"
)
#> ℹ [2026-03-08 08:07:23] Start SingleR annotation
#> ℹ [2026-03-08 08:14:14] Data type is log-normalized
#> ℹ [2026-03-08 08:14:14] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-08 08:14:15] Data type is log-normalized
#> ℹ [2026-03-08 08:14:15] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-08 08:14:19] Perform "SingleRCluster"
#> ✔ [2026-03-08 08:14:20] SingleR annotation completed
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
#> ℹ [2026-03-08 08:14:21] Start SingleR annotation
#> ℹ [2026-03-08 08:14:21] Data type is log-normalized
#> ℹ [2026-03-08 08:14:21] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-08 08:14:24] Data type is log-normalized
#> ℹ [2026-03-08 08:14:24] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-08 08:14:30] Perform "SingleRCell"
#> ✔ [2026-03-08 08:14:34] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "CellType")
)
```
