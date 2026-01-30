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
#> ℹ [2026-01-30 17:24:38] Start standard scop workflow...
#> ℹ [2026-01-30 17:24:38] Checking a list of <Seurat>...
#> ! [2026-01-30 17:24:38] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:24:38] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:24:41] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:24:41] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:24:42] Number of available HVF: 2000
#> ℹ [2026-01-30 17:24:42] Finished check
#> ℹ [2026-01-30 17:24:42] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:24:43] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:24:44] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:24:44] Reorder clusters...
#> ℹ [2026-01-30 17:24:44] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:24:44] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:24:50] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:24:55] Run scop standard workflow completed

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
#> ℹ [2026-01-30 17:24:55] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2026-01-30 17:24:55] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2026-01-30 17:24:56] Checking a list of <Seurat>...
#> ℹ [2026-01-30 17:24:57] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:24:57] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ℹ [2026-01-30 17:24:57] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:24:57] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ℹ [2026-01-30 17:24:58] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:24:58] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ℹ [2026-01-30 17:24:58] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:24:58] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ℹ [2026-01-30 17:24:59] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2026-01-30 17:24:59] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2026-01-30 17:24:59] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:24:59] Number of available HVF: 2000
#> ℹ [2026-01-30 17:25:00] Finished check

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-30 17:25:02] Start standard scop workflow...
#> ℹ [2026-01-30 17:25:03] Checking a list of <Seurat>...
#> ! [2026-01-30 17:25:03] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:25:03] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:25:05] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:25:06] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:25:06] Number of available HVF: 2000
#> ℹ [2026-01-30 17:25:06] Finished check
#> ℹ [2026-01-30 17:25:07] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:25:07] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:25:08] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:25:08] Reorder clusters...
#> ℹ [2026-01-30 17:25:08] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:25:08] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:25:13] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:25:18] Run scop standard workflow completed
pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCluster"
)
#> ℹ [2026-01-30 17:25:48] Data type is log-normalized
#> ℹ [2026-01-30 17:25:48] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-01-30 17:25:50] Data type is log-normalized
#> ℹ [2026-01-30 17:25:50] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-01-30 17:25:52] Perform selectFeatures
#> ℹ [2026-01-30 17:25:53] Perform indexCluster
#> ℹ [2026-01-30 17:25:53] Perform scmapCluster
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
#> ℹ [2026-01-30 17:25:54] Data type is log-normalized
#> ℹ [2026-01-30 17:25:54] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-01-30 17:25:55] Data type is log-normalized
#> ℹ [2026-01-30 17:25:55] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-01-30 17:25:59] Perform selectFeatures
#> ℹ [2026-01-30 17:25:59] Perform indexCell
#> ℹ [2026-01-30 17:26:00] Perform scmapCell
#> ℹ [2026-01-30 17:26:01] Perform scmapCell2Cluster
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)
```
