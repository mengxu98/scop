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
options(log_message.verbose = FALSE)
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#> This message is displayed once every 8 hours.
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "mouse"
)
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
