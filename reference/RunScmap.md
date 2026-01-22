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
#> ℹ [2026-01-22 04:12:10] Start standard scop workflow...
#> ℹ [2026-01-22 04:12:11] Checking a list of <Seurat>...
#> ! [2026-01-22 04:12:11] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-22 04:12:11] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:12:14] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:12:14] Use the separate HVF from srt_list
#> ℹ [2026-01-22 04:12:14] Number of available HVF: 2000
#> ℹ [2026-01-22 04:12:15] Finished check
#> ℹ [2026-01-22 04:12:15] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-22 04:12:15] Perform pca linear dimension reduction
#> ℹ [2026-01-22 04:12:17] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-22 04:12:17] Reorder clusters...
#> ℹ [2026-01-22 04:12:17] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-22 04:12:17] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-22 04:12:22] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-22 04:12:27] Run scop standard workflow completed

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
#> ℹ [2026-01-22 04:12:27] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-01-22 04:12:27] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-22 04:12:28] Checking a list of <Seurat>...
#> ℹ [2026-01-22 04:12:28] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-22 04:12:28] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ℹ [2026-01-22 04:12:29] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-22 04:12:29] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ℹ [2026-01-22 04:12:29] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-22 04:12:29] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ℹ [2026-01-22 04:12:30] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-22 04:12:30] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ℹ [2026-01-22 04:12:31] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-22 04:12:31] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2026-01-22 04:12:31] Use the separate HVF from srt_list
#> ℹ [2026-01-22 04:12:31] Number of available HVF: 2000
#> ℹ [2026-01-22 04:12:32] Finished check

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-22 04:12:34] Start standard scop workflow...
#> ℹ [2026-01-22 04:12:35] Checking a list of <Seurat>...
#> ! [2026-01-22 04:12:35] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-22 04:12:35] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:12:37] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:12:37] Use the separate HVF from srt_list
#> ℹ [2026-01-22 04:12:38] Number of available HVF: 2000
#> ℹ [2026-01-22 04:12:38] Finished check
#> ℹ [2026-01-22 04:12:38] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-22 04:12:38] Perform pca linear dimension reduction
#> ℹ [2026-01-22 04:12:39] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-22 04:12:39] Reorder clusters...
#> ℹ [2026-01-22 04:12:40] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-22 04:12:40] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-22 04:12:44] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-22 04:12:49] Run scop standard workflow completed
pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCluster"
)
#> ℹ [2026-01-22 04:13:18] Data type is log-normalized
#> ℹ [2026-01-22 04:13:18] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-01-22 04:13:20] Data type is log-normalized
#> ℹ [2026-01-22 04:13:20] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-01-22 04:13:22] Perform selectFeatures
#> ℹ [2026-01-22 04:13:22] Perform indexCluster
#> ℹ [2026-01-22 04:13:23] Perform scmapCluster
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
#> ℹ [2026-01-22 04:13:24] Data type is log-normalized
#> ℹ [2026-01-22 04:13:24] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-01-22 04:13:25] Data type is log-normalized
#> ℹ [2026-01-22 04:13:25] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-01-22 04:13:29] Perform selectFeatures
#> ℹ [2026-01-22 04:13:29] Perform indexCell
#> ℹ [2026-01-22 04:13:30] Perform scmapCell
#> ℹ [2026-01-22 04:13:31] Perform scmapCell2Cluster
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)
```
