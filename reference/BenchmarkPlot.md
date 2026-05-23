# Plot benchmark metrics

Visualize benchmark results stored in a `Seurat` object or a summary
`data.frame`. Per-cell metrics such as LISI are shown as feature plots
and/or boxplots. Summary metrics such as integration or mapping
benchmark scores are shown as barplots or a `funkyheatmap`.

## Usage

``` r
BenchmarkPlot(
  srt = NULL,
  data = NULL,
  features = NULL,
  metrics = NULL,
  tool_name = NULL,
  reduction = NULL,
  plot_type = c("auto", "feature", "boxplot", "bar", "funkyheatmap"),
  plot_boxplot = TRUE,
  boxplot_jitter = FALSE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- data:

  Optional summary benchmark `data.frame` containing at least `metric`
  and `value`, and optionally `method` and `workflow`.

- features:

  Metadata columns containing per-cell benchmark scores. Default is
  `NULL`.

- metrics:

  One or more summary metric names to visualize. Default is `NULL`,
  which uses all available summary metrics.

- tool_name:

  Tool entries created by benchmark-related workflows. This can be a
  character vector. For per-cell metrics, benchmark columns are resolved
  from tool entries that contain `colnames`; for summary metrics,
  entries containing `summary` or `metrics$summary` are used.

- reduction:

  Dimensional reduction used for per-cell feature plots. Default is
  `NULL`, which uses the reduction stored in `tool_name` when available,
  otherwise
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- plot_type:

  Plot type. One of `"auto"`, `"feature"`, `"boxplot"`, `"bar"`, or
  `"funkyheatmap"`.

- plot_boxplot:

  Whether to add the summary boxplot when per-cell metrics are shown.
  Default is `TRUE`.

- boxplot_jitter:

  Whether to overlay jittered points on the boxplot. Default is `FALSE`.

- combine:

  Combine plots into a single `patchwork` object. If `FALSE`, return a
  list of ggplot objects.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

- pt.size:

  The size of the points in the plot.

- pt.alpha:

  The transparency of the data points. Default is `1`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  The message to print.

## Value

A ggplot, patchwork plot, or `funkyheatmap` object depending on the
selected mode. If `combine = FALSE` in per-cell mode, a named list of
plots is returned.

## Examples

``` r
metrics_df <- data.frame(
  method = c("Raw", "Raw", "Harmony", "Harmony"),
  metric = c("batch_ASW_mixing", "celltype_ASW", "batch_ASW_mixing", "celltype_ASW"),
  value = c(0.42, 0.71, 0.68, 0.66)
)
BenchmarkPlot(
  data = metrics_df,
  plot_type = "bar"
)


data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub[["MethodA_batch_LISI"]] <-
  seq_len(ncol(pbmcmultiome_sub)) / ncol(pbmcmultiome_sub)
pbmcmultiome_sub[["MethodB_batch_LISI"]] <-
  rev(pbmcmultiome_sub[["MethodA_batch_LISI", drop = TRUE]])
BenchmarkPlot(
  pbmcmultiome_sub,
  features = c("MethodA_batch_LISI", "MethodB_batch_LISI"),
  plot_type = "boxplot"
)

```
