# Rename clusters for the Seurat object

Rename clusters for the Seurat object

## Usage

``` r
RenameClusters(
  object,
  group.by,
  nameslist = list(),
  name = "newclusters",
  keep_levels = FALSE,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- group.by:

  The old group used to rename cells.

- nameslist:

  A named list of new cluster value.

- name:

  New cluster stored in the Seurat object.

- keep_levels:

  If the old group is a factor, keep the order of the levels.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Examples

``` r
data(pancreas_sub)

# Rename all clusters
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-13 22:03:05] Start standard processing workflow...
#> ℹ [2026-09-13 22:03:05] Checking a list of <Seurat>...
#> ! [2026-09-13 22:03:05] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-13 22:03:05] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:03:05] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-13 22:03:05] Use the separate HVF from `srt_list`
#> ℹ [2026-09-13 22:03:05] Number of available HVF: 2000
#> ℹ [2026-09-13 22:03:05] Finished check
#> ℹ [2026-09-13 22:03:05] Perform `ScaleData()`
#> ℹ [2026-09-13 22:03:05] Perform pca linear dimension reduction
#> ℹ [2026-09-13 22:03:05] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-13 22:03:06] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-13 22:03:06] Reorder clusters...
#> ℹ [2026-09-13 22:03:06] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-13 22:03:06] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-13 22:03:12] Standard processing workflow completed
levels(pancreas_sub@meta.data[["SubCellType"]]) <- unique(
  pancreas_sub@meta.data[["SubCellType"]]
)
pancreas_sub <- RenameClusters(
  pancreas_sub,
  group.by = "SubCellType",
  nameslist = letters[1:8]
)
CellDimPlot(pancreas_sub, "newclusters")


# Rename specified clusters
pancreas_sub <- RenameClusters(pancreas_sub,
  group.by = "SubCellType",
  nameslist = list("a" = "Alpha", "b" = "Beta")
)
CellDimPlot(pancreas_sub, "newclusters")


# Merge and rename clusters
pancreas_sub <- RenameClusters(
  pancreas_sub,
  group.by = "SubCellType",
  nameslist = list(
    "EndocrineClusters" = c("Alpha", "Beta", "Epsilon", "Delta")
  ),
  name = "Merged",
  keep_levels = TRUE
)
CellDimPlot(pancreas_sub, "Merged")
```
