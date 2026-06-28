# Annotate single cells using scmap.

Annotate single cells using scmap.

## Usage

``` r
RunScmap(
  srt_query,
  srt_ref,
  ref_group = NULL,
  query_assay = "RNA",
  ref_assay = "RNA",
  method = "scmapCluster",
  nfeatures = 500,
  threshold = 0.5,
  k = 10,
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- ref_group:

  A character vector specifying the column name in the `srt_ref`
  metadata that represents the cell grouping.

- query_assay:

  A character vector specifying the assay to be used for the query data.
  Default is the default assay of the `srt_query` object.

- ref_assay:

  A character vector specifying the assay to be used for the reference
  data. Default is the default assay of the `srt_ref` object.

- method:

  The method to be used for scmap analysis. Can be any of
  `"scmapCluster"` or `"scmapCell"`. Default is `"scmapCluster"`.

- nfeatures:

  The number of top features to be selected. Default is `500`.

- threshold:

  The threshold value on similarity to determine if a cell is assigned
  to a cluster. This should be a value between `0` and `1`. Default is
  `0.5`.

- k:

  Number of clusters per group for k-means clustering when `method` is
  `"scmapCell"`. Default is `10`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-06-28 21:16:16] Start standard processing workflow...
#> ℹ [2026-06-28 21:16:17] Checking a list of <Seurat>...
#> ! [2026-06-28 21:16:17] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-28 21:16:17] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 21:16:17] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 21:16:17] Use the separate HVF from `srt_list`
#> ℹ [2026-06-28 21:16:17] Number of available HVF: 2000
#> ℹ [2026-06-28 21:16:17] Finished check
#> ℹ [2026-06-28 21:16:17] Perform `ScaleData()`
#> ℹ [2026-06-28 21:16:17] Perform pca linear dimension reduction
#> ℹ [2026-06-28 21:16:19] Use stored estimated dimensions 1:27 for Standardpca
#> ℹ [2026-06-28 21:16:19] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-28 21:16:19] Reorder clusters...
#> ℹ [2026-06-28 21:16:19] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-28 21:16:19] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-28 21:16:27] Standard processing workflow completed

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
#> ℹ [2026-06-28 21:16:27] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-06-28 21:16:28] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-06-28 21:16:29] Checking a list of <Seurat>...
#> ℹ [2026-06-28 21:16:30] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-06-28 21:16:30] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-06-28 21:16:30] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-06-28 21:16:30] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-06-28 21:16:30] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-06-28 21:16:30] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-06-28 21:16:30] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-06-28 21:16:30] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-06-28 21:16:31] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-06-28 21:16:31] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-06-28 21:16:31] Use the separate HVF from `srt_list`
#> ℹ [2026-06-28 21:16:31] Number of available HVF: 2000
#> ℹ [2026-06-28 21:16:32] Finished check
#> Warning: Key ‘StandardpcaUMAP2D_’ taken, using ‘standardpcaumap2d_’ instead

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-28 21:16:32] Start standard processing workflow...
#> ℹ [2026-06-28 21:16:33] Checking a list of <Seurat>...
#> ! [2026-06-28 21:16:33] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-28 21:16:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 21:16:33] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 21:16:33] Use the separate HVF from `srt_list`
#> ℹ [2026-06-28 21:16:33] Number of available HVF: 2000
#> ℹ [2026-06-28 21:16:33] Finished check
#> ℹ [2026-06-28 21:16:33] Perform `ScaleData()`
#> ℹ [2026-06-28 21:16:33] Perform pca linear dimension reduction
#> ℹ [2026-06-28 21:16:34] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-28 21:16:34] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-28 21:16:35] Reorder clusters...
#> ℹ [2026-06-28 21:16:35] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-28 21:16:35] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-28 21:16:42] Standard processing workflow completed
pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCluster"
)
#> ℹ [2026-06-28 21:17:16] Data type is log-normalized
#> ℹ [2026-06-28 21:17:16] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-06-28 21:17:17] Data type is log-normalized
#> ℹ [2026-06-28 21:17:17] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-06-28 21:17:18] Perform selectFeatures
#> ℹ [2026-06-28 21:17:18] Perform indexCluster
#> ℹ [2026-06-28 21:17:19] Perform scmapCluster
#> Warning: Features Mt-atp6, Mt-co1, Mt-co2, Mt-co3, Mt-nd1, Mt-nd2, Mt-nd4, Mt-nd4l, Mt-nd5 are not present in the 'SCESet' object and therefore were not set.
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)


pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCell"
)
#> ℹ [2026-06-28 21:17:20] Data type is log-normalized
#> ℹ [2026-06-28 21:17:20] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-06-28 21:17:20] Data type is log-normalized
#> ℹ [2026-06-28 21:17:20] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-06-28 21:17:21] Perform selectFeatures
#> ℹ [2026-06-28 21:17:22] Perform indexCell
#> ℹ [2026-06-28 21:17:22] Perform scmapCell
#> ℹ [2026-06-28 21:17:23] Perform scmapCell2Cluster
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)
```
