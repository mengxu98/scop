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
  k = 10
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

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-03-11 17:48:44] Start standard scop workflow...
#> ℹ [2026-03-11 17:48:44] Checking a list of <Seurat>...
#> ! [2026-03-11 17:48:45] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-11 17:48:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-11 17:48:47] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-11 17:48:48] Use the separate HVF from `srt_list`
#> ℹ [2026-03-11 17:48:48] Number of available HVF: 2000
#> ℹ [2026-03-11 17:48:48] Finished check
#> ℹ [2026-03-11 17:48:48] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-11 17:48:49] Perform pca linear dimension reduction
#> ℹ [2026-03-11 17:48:50] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-11 17:48:50] Reorder clusters...
#> ℹ [2026-03-11 17:48:50] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-11 17:48:50] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-11 17:48:55] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-11 17:49:01] Run scop standard workflow completed

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
#> ℹ [2026-03-11 17:49:01] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-03-11 17:49:01] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-03-11 17:49:02] Checking a list of <Seurat>...
#> ℹ [2026-03-11 17:49:02] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-03-11 17:49:02] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ℹ [2026-03-11 17:49:03] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-03-11 17:49:03] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ℹ [2026-03-11 17:49:03] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-03-11 17:49:03] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ℹ [2026-03-11 17:49:04] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-03-11 17:49:04] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ℹ [2026-03-11 17:49:04] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-03-11 17:49:05] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-03-11 17:49:05] Use the separate HVF from `srt_list`
#> ℹ [2026-03-11 17:49:05] Number of available HVF: 2000
#> ℹ [2026-03-11 17:49:06] Finished check

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-03-11 17:49:08] Start standard scop workflow...
#> ℹ [2026-03-11 17:49:09] Checking a list of <Seurat>...
#> ! [2026-03-11 17:49:09] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-11 17:49:09] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-11 17:49:11] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-11 17:49:11] Use the separate HVF from `srt_list`
#> ℹ [2026-03-11 17:49:12] Number of available HVF: 2000
#> ℹ [2026-03-11 17:49:12] Finished check
#> ℹ [2026-03-11 17:49:12] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-11 17:49:12] Perform pca linear dimension reduction
#> ℹ [2026-03-11 17:49:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-11 17:49:13] Reorder clusters...
#> ℹ [2026-03-11 17:49:13] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-11 17:49:13] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-11 17:49:18] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-11 17:49:23] Run scop standard workflow completed
pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCluster"
)
#> ℹ [2026-03-11 17:49:52] Data type is log-normalized
#> ℹ [2026-03-11 17:49:52] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-11 17:49:53] Data type is log-normalized
#> ℹ [2026-03-11 17:49:53] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-11 17:49:56] Perform selectFeatures
#> ℹ [2026-03-11 17:49:56] Perform indexCluster
#> ℹ [2026-03-11 17:49:56] Perform scmapCluster
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
#> ℹ [2026-03-11 17:49:57] Data type is log-normalized
#> ℹ [2026-03-11 17:49:57] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-03-11 17:49:59] Data type is log-normalized
#> ℹ [2026-03-11 17:49:59] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-03-11 17:50:01] Perform selectFeatures
#> ℹ [2026-03-11 17:50:01] Perform indexCell
#> ℹ [2026-03-11 17:50:02] Perform scmapCell
#> ℹ [2026-03-11 17:50:03] Perform scmapCell2Cluster
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)
```
