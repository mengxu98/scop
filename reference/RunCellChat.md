# Run CellChat analysis

Run CellChat analysis

## Usage

``` r
RunCellChat(
  srt,
  group.by,
  species = c("Homo_sapiens", "Mus_musculus", "zebrafish"),
  split.by = NULL,
  annotation_selected = NULL,
  group_column = NULL,
  group_cmp = NULL,
  thresh = 0.05,
  min.cells = 10,
  do.fast = FALSE,
  backend = c("cpp", "r"),
  assay = NULL,
  layer = "data",
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- species:

  The species of the data, either `"Homo_sapiens"`, `"Mus_musculus"`, or
  `"zebrafish"`.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- annotation_selected:

  A vector of cell annotations of interest for running the `CellChat`
  analysis. If not provided, all cell types will be considered.

- group_column:

  Name of the metadata column in the `Seurat` object that defines
  conditions or groups.

- group_cmp:

  A list of pairwise condition comparisons for differential `CellChat`
  analysis.

- thresh:

  The threshold for computing centrality scores. Default is `0.05`.

- min.cells:

  the minmum number of expressed cells required for the genes that are
  considered for cell-cell communication analysis. Default is `10`.

- do.fast:

  Whether to use CellChat's fast Wilcoxon implementation backed by
  `presto`. Set to `TRUE` only when `presto` is installed.

- backend:

  Backend used for scop post-processing and unified CCC table
  aggregation. Upstream CellChat inference is unchanged.

- assay:

  Which assay to use. If `NULL`, the default assay of the `Seurat`
  object will be used.

- layer:

  The layer to use for the expression data. Default is `"data"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with `CellChat` results stored in
`srt@tools[["CellChat"]]`.

## References

[CellChat](https://github.com/jinworks/CellChat)

## See also

[CCCHeatmap](https://mengxu98.github.io/scop/reference/CCCHeatmap.md),
[CCCStatPlot](https://mengxu98.github.io/scop/reference/CCCStatPlot.md),
[CCCNetworkPlot](https://mengxu98.github.io/scop/reference/CCCNetworkPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-28 12:26:26] Start standard processing workflow...
#> ℹ [2026-06-28 12:26:27] Checking a list of <Seurat>...
#> ! [2026-06-28 12:26:27] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-28 12:26:27] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 12:26:27] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 12:26:28] Use the separate HVF from `srt_list`
#> ℹ [2026-06-28 12:26:28] Number of available HVF: 2000
#> ℹ [2026-06-28 12:26:28] Finished check
#> ℹ [2026-06-28 12:26:28] Perform `ScaleData()`
#> ℹ [2026-06-28 12:26:28] Perform pca linear dimension reduction
#> ℹ [2026-06-28 12:26:29] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-28 12:26:29] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-28 12:26:29] Reorder clusters...
#> ℹ [2026-06-28 12:26:29] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-28 12:26:29] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-28 12:26:36] Standard processing workflow completed
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "Mus_musculus"
)
#> ℹ [2026-06-28 12:26:36] Start CellChat analysis
#> Error in loadNamespace(name): there is no package called ‘CellChat’

CCCNetworkPlot(
  pancreas_sub,
  method = "CellChat",
  plot_type = "bipartite"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCHeatmap(
  pancreas_sub,
  method = "CellChat",
  plot_type = "heatmap"
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`

CCCStatPlot(
  pancreas_sub,
  method = "CellChat",
  plot_type = "violin",
  top_n = 50
)
#> Error in get_dataset_object(srt, condition = condition, dataset = dataset): Unable to determine which CellChat object to plot. Please specify
#> `condition`
```
