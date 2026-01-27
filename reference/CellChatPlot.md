# Plot CellChat analysis results

CellChatPlot creates various visualizations for CellChat analysis
results stored in a Seurat object.

## Usage

``` r
CellChatPlot(
  srt,
  plot_type = "aggregate",
  condition = NULL,
  pathway = NULL,
  dirpath = NULL,
  output_format = "pdf",
  top_n = 10,
  base_height = 1,
  base_width = 1,
  res = 300,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object that has been processed with
  [RunCellChat](https://mengxu98.github.io/scop/reference/RunCellChat.md).

- plot_type:

  Type of plot to create. Options: `"aggregate"`, `"pathway"`,
  `"comparison"`, `"heatmap"`, `"circle"`, `"bubble"`, `"gene"`.

- condition:

  Condition to plot (if multiple conditions exist).

- pathway:

  Specific pathway to visualize (for pathway, bubble, and gene plots).
  If `NULL`, uses top pathways.

- dirpath:

  Directory to save plots.

- output_format:

  Format of output figure: `"png"` or `"pdf"`. Default is `"png"`.

- top_n:

  Number of top pathways to use for plotting. Default is `10`.

- base_height:

  Base height multiplier for all plots. Default is `1`.

- base_width:

  Base width multiplier for all plots. Default is `1`.

- res:

  Resolution for PNG output. Default is `300`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[RunCellChat](https://mengxu98.github.io/scop/reference/RunCellChat.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-27 07:11:22] Start standard scop workflow...
#> ℹ [2026-01-27 07:11:23] Checking a list of <Seurat>...
#> ! [2026-01-27 07:11:23] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-27 07:11:23] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 07:11:25] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 07:11:26] Use the separate HVF from srt_list
#> ℹ [2026-01-27 07:11:26] Number of available HVF: 2000
#> ℹ [2026-01-27 07:11:27] Finished check
#> ℹ [2026-01-27 07:11:27] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-27 07:11:28] Perform pca linear dimension reduction
#> ℹ [2026-01-27 07:11:29] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-27 07:11:29] Reorder clusters...
#> First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#> This message is displayed once every 8 hours.
#> ℹ [2026-01-27 07:11:29] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-27 07:11:29] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-27 07:11:32] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-27 07:11:34] Run scop standard workflow completed
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "Mus_musculus"
)
#> ℹ [2026-01-27 07:11:34] Start CellChat analysis
#> Error in loadNamespace(x): there is no package called ‘CellChat’

CellChatPlot(pancreas_sub, plot_type = "aggregate")
#> Error in CellChatPlot(pancreas_sub, plot_type = "aggregate"): No CellChat results found in <Seurat>. Please run `RunCellChat()` first

CellChatPlot(pancreas_sub, plot_type = "pathway")
#> Error in CellChatPlot(pancreas_sub, plot_type = "pathway"): No CellChat results found in <Seurat>. Please run `RunCellChat()` first

CellChatPlot(pancreas_sub, plot_type = "bubble")
#> Error in CellChatPlot(pancreas_sub, plot_type = "bubble"): No CellChat results found in <Seurat>. Please run `RunCellChat()` first

CellChatPlot(pancreas_sub, plot_type = "gene")
#> Error in CellChatPlot(pancreas_sub, plot_type = "gene"): No CellChat results found in <Seurat>. Please run `RunCellChat()` first

CellChatPlot(pancreas_sub, plot_type = "heatmap")
#> Error in CellChatPlot(pancreas_sub, plot_type = "heatmap"): No CellChat results found in <Seurat>. Please run `RunCellChat()` first
```
