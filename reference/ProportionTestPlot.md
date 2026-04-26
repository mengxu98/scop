# Proportion Test Plot

Generate differential abundance visualizations from outputs of
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).
Main display styles are effect-size view and UMAP projection view, while
legacy plot types remain callable for compatibility.

## Usage

``` r
ProportionTestPlot(
  srt,
  comparison = NULL,
  proportion_method = NULL,
  result_level = c("group", "neighborhood"),
  plot_type = c(
    "effect",
    "umap",
    "volcano",
    "manhattan",
    "ring",
    "milo_beeswarm",
    "milo_graph",
    "sccoda_forest"
  ),
  umap_mode = c("discrete", "continuous"),
  reduction = "UMAP",
  projection_args = list(),
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  order_by = c("value", "name"),
  palette = "RdBu",
  palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.sig = "red",
  cols.ns = "grey",
  cols.increase = "#d7301f",
  cols.decrease = "#2b8cbe",
  effect_color_mode = c("directional", "classic"),
  nlabel = 5,
  features_label = NULL,
  label = FALSE,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  aspect.ratio = NULL,
  xlab = "Cell Type",
  ylab = "log2 (FD)",
  theme_use = "theme_scop",
  theme_args = list(),
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.title = "Significance",
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object containing proportion test results.

- comparison:

  A character vector specifying which comparisons to plot. If `NULL`,
  all stored comparisons are plotted.

- proportion_method:

  Optional method to select from method-layer results. If `NULL`,
  active/most recent method is used.

- result_level:

  Result level for plotting: `"group"` or `"neighborhood"`.

- plot_type:

  Plot type. Main recommended values are `"effect"` and `"umap"`. Legacy
  values are retained for compatibility: `"volcano"`, `"manhattan"`,
  `"ring"`, `"milo_beeswarm"`, `"milo_graph"`, and `"sccoda_forest"`.

- umap_mode:

  UMAP mapping mode when `plot_type = "umap"`. `"discrete"` maps cells
  to `Increased/Decreased/NS`; `"continuous"` maps group-level
  `obs_log2FD`.

- reduction:

  Reduction name used for UMAP projection.

- projection_args:

  Additional arguments passed through to `CellDimPlot` (discrete mode)
  or `FeatureDimPlot` (continuous mode).

- FDR_threshold:

  FDR value cutoff for significance.

- log2FD_threshold:

  Absolute value of log2FD cutoff for significance.

- order_by:

  Method to order clusters. Options: `"name"`, `"value"`.

- palette:

  Palette name for continuous effect coloring.

- palcolor:

  Custom colors for `palette`.

- group_palette:

  Palette name for cluster/group coloring.

- group_palcolor:

  Custom colors for `group_palette`.

- pt.size:

  Point size.

- pt.alpha:

  Point transparency.

- cols.sig:

  Legacy color for significant/credible points and intervals.

- cols.ns:

  Color for non-significant points and intervals.

- cols.increase:

  Color for increased DA groups in directional effect mode.

- cols.decrease:

  Color for decreased DA groups in directional effect mode.

- effect_color_mode:

  Effect color style: `"directional"` (default) or `"classic"`.

- nlabel:

  Number of labels when `label = TRUE` and `features_label = NULL`.

- features_label:

  Character vector specifying labels to draw.

- label:

  Whether to draw labels.

- label.fg:

  Label foreground color.

- label.bg:

  Label background color.

- label.bg.r:

  Label background rounding radius.

- label.size:

  Label size.

- aspect.ratio:

  Aspect ratio of each panel.

- xlab:

  X-axis label.

- ylab:

  Y-axis label.

- theme_use:

  Theme used. Can be a character string or a theme function.

- theme_args:

  Additional arguments passed to `theme_use`.

- legend.position:

  Legend position.

- legend.direction:

  Legend direction.

- legend.title:

  Legend title.

- combine:

  Whether to combine into one patchwork object.

- nrow:

  Number of rows for combined plot.

- ncol:

  Number of columns for combined plot.

- byrow:

  Arrange panels by row.

- seed:

  Random seed used for jitter and label ranking ties.

## See also

[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)
