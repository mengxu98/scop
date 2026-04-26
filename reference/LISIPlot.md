# Plot LISI scores

Backward-compatible wrapper around
[`BenchmarkPlot()`](https://mengxu98.github.io/scop/reference/BenchmarkPlot.md)
for LISI scores. Visualize LISI scores on a dimensional reduction and
compare methods with a summary boxplot.

## Usage

``` r
LISIPlot(
  srt,
  features = NULL,
  tool_name = NULL,
  reduction = NULL,
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

- features:

  Metadata columns containing LISI scores. Default is `NULL`, which will
  use columns stored in `tool_name`, or all metadata columns ending with
  `"_LISI"` when `tool_name` is `NULL`.

- tool_name:

  Tool entry created by
  [`RunLISI()`](https://mengxu98.github.io/scop/reference/RunLISI.md).
  Default is `NULL`.

- reduction:

  Dimensional reduction used for feature plots. If `NULL`, the reduction
  recorded in `tool_name` is used when available; otherwise
  [`DefaultReduction()`](https://mengxu98.github.io/scop/reference/DefaultReduction.md)
  is used.

- plot_boxplot:

  Whether to add boxplots. Default is `TRUE`.

- boxplot_jitter:

  Whether to overlay jittered points on boxplots. Default is `FALSE`.

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

If `combine = TRUE`, returns a combined `patchwork` plot. If
`combine = FALSE`, returns a named list of ggplot objects.

## See also

[RunLISI](https://mengxu98.github.io/scop/reference/RunLISI.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
