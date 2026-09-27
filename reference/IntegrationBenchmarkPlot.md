# Plot integration benchmark results

Visualize a `Seurat` object returned by
[`RunIntegrationBenchmark()`](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md).
`"box"` compares per-cell LISI scores across methods, `"heatmap"` shows
scaled summary metrics, `"scatter"` plots the bio-vs-batch trade-off,
and `"umap"` draws batch and cell-type embeddings. `"auto"` returns all
of these plots.

## Usage

``` r
IntegrationBenchmarkPlot(
  object,
  plot_type = c("auto", "box", "heatmap", "scatter", "umap"),
  tool_name = "IntegrationBenchmark",
  metrics = NULL,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  pt.size = NULL,
  pt.alpha = 1,
  boxplot_jitter = FALSE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  verbose = TRUE,
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object from
  [`RunIntegrationBenchmark()`](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md).
  `"box"` can also plot any metadata columns ending in `"_LISI"`.

- plot_type:

  Plot type, or `"auto"` for a named list of all views.

- tool_name:

  `srt@tools` entry created by
  [`RunIntegrationBenchmark()`](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md).

- metrics:

  Optional metric names used in the heatmap.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- boxplot_jitter:

  Whether to overlay jittered points on the LISI boxes.

- combine:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- nrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- ncol:

  Number of columns of the combined plot.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  The message to print.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A ggplot, patchwork, or a named list of plots when `plot_type = "auto"`.

## See also

[RunIntegrationBenchmark](https://mengxu98.github.io/scop/reference/RunIntegrationBenchmark.md),
[RunLISI](https://mengxu98.github.io/scop/reference/RunLISI.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(panc8_sub)
panc8_sub <- RunIntegrationBenchmark(
  panc8_sub,
  batch = "tech",
  celltype = "celltype",
  methods = c("Uncorrected", "Harmony"),
  nHVF = 500,
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:10,
  perplexity = 10
)
IntegrationBenchmarkPlot(panc8_sub)
IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")
IntegrationBenchmarkPlot(panc8_sub, plot_type = "heatmap")
IntegrationBenchmarkPlot(panc8_sub, plot_type = "scatter")
} # }
```
