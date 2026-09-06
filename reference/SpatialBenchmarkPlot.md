# Plot spatial domain clustering benchmarks

Visualize a `spatial_benchmark_result` from
[`RunSpatialBenchmark()`](https://mengxu98.github.io/scop/reference/RunSpatialBenchmark.md),
or a summary `data.frame`. Spatial results default to a
publication-oriented overview that pairs clustering quality with runtime
and peak-memory efficiency. For scRNA integration metrics such as LISI,
use
[`IntegrationBenchmarkPlot()`](https://mengxu98.github.io/scop/reference/IntegrationBenchmarkPlot.md).

## Usage

``` r
SpatialBenchmarkPlot(
  srt = NULL,
  data = NULL,
  features = NULL,
  metrics = NULL,
  tool_name = NULL,
  reduction = NULL,
  plot_type = c("auto", "overview", "quality", "efficiency", "heatmap", "feature",
    "boxplot", "bar", "funkyheatmap"),
  sort_by = c("quality", "method", "runtime", "memory"),
  show_values = TRUE,
  show_status = TRUE,
  resource_scale = c("auto", "linear", "log10"),
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

  Optional `spatial_benchmark_result` / `benchmark_result` object or
  summary `data.frame` containing at least `metric` and `value`, and
  optionally `method`, `workflow`, and `direction`.

- features:

  Metadata columns containing per-cell benchmark scores.

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

  Plot type. `"overview"`, `"quality"`, `"efficiency"`, and `"heatmap"`
  consume spatial benchmark results. Existing `"feature"`, `"boxplot"`,
  `"bar"`, and `"funkyheatmap"` modes remain supported.

- sort_by:

  Method ordering for spatial benchmark plots. `"quality"` sorts by the
  mean selected quality metric; other choices sort by method, runtime,
  or peak memory.

- show_values:

  Whether to print raw metric values on quality and heatmap panels.

- show_status:

  Whether the overview should add a status strip for failed,
  unavailable, or timed-out methods.

- resource_scale:

  Resource-axis transformation. `"auto"` independently uses log10 for
  runtime or memory when positive values span at least tenfold.

- plot_boxplot:

  Whether to add the summary boxplot when per-cell metrics are shown.

- boxplot_jitter:

  Whether to overlay jittered points on the boxplot.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  The message to print.

## Value

A ggplot, patchwork plot, or `funkyheatmap` object depending on the
selected mode. If `combine = FALSE` in per-cell mode, a named list of
plots is returned.

## See also

[RunSpatialBenchmark](https://mengxu98.github.io/scop/reference/RunSpatialBenchmark.md),
[IntegrationBenchmarkPlot](https://mengxu98.github.io/scop/reference/IntegrationBenchmarkPlot.md)

## Examples

``` r
metrics_df <- data.frame(
  method = c("Raw", "Raw", "Harmony", "Harmony"),
  metric = c("batch_ASW_mixing", "celltype_ASW", "batch_ASW_mixing", "celltype_ASW"),
  value = c(0.42, 0.71, 0.68, 0.66)
)
SpatialBenchmarkPlot(
  data = metrics_df,
  plot_type = "bar"
)


data("pbmcmultiome_sub", package = "scop")
pbmcmultiome_sub[["MethodA_batch_LISI"]] <-
  seq_len(ncol(pbmcmultiome_sub)) / ncol(pbmcmultiome_sub)
pbmcmultiome_sub[["MethodB_batch_LISI"]] <-
  rev(pbmcmultiome_sub[["MethodA_batch_LISI", drop = TRUE]])
SpatialBenchmarkPlot(
  pbmcmultiome_sub,
  features = c("MethodA_batch_LISI", "MethodB_batch_LISI"),
  plot_type = "boxplot"
)

```
