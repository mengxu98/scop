# Proportion Test Plot

Generate proportion test plots based on the results from
[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md).

## Usage

``` r
ProportionTestPlot(
  srt,
  comparison = NULL,
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  order_by = c("value", "name"),
  pt.size = 1,
  pt.alpha = 1,
  cols.sig = "red",
  cols.ns = "grey",
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
  byrow = TRUE
)
```

## Arguments

- srt:

  A Seurat object containing proportion test results.

- comparison:

  A character string specifying which comparison to plot. If NULL, plots
  all comparisons.

- FDR_threshold:

  FDR value cutoff for significance.

- log2FD_threshold:

  Absolute value of log2FD cutoff for significance.

- order_by:

  Method to order clusters. Options: "name" (alphabetical), "value" (by
  log2FD value).

- pt.size:

  The size of the points.

- pt.alpha:

  The transparency of the points.

- cols.sig:

  Color for significant points and confidence intervals.

- cols.ns:

  Color for non-significant points and confidence intervals.

- aspect.ratio:

  The aspect ratio of the plot.

- xlab:

  A character string specifying the x-axis label.

- ylab:

  A character string specifying the y-axis label.

- theme_use:

  A character string specifying the theme to use for the plot.

- theme_args:

  A list of theme arguments to pass to the `theme_use` function.

- legend.position:

  Position of the legend.

- legend.direction:

  Direction of the legend.

- legend.title:

  Title of the legend.

- combine:

  Whether to combine the plots for each comparison into a single plot.

- nrow:

  An integer value specifying the number of rows in the combined plot.

- ncol:

  An integer value specifying the number of columns in the combined
  plot.

- byrow:

  Whether to arrange the plots by row in the combined plot.

## See also

[RunProportionTest](https://mengxu98.github.io/scop/reference/RunProportionTest.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunProportionTest(
  pancreas_sub,
  group.by = "CellType",
  split.by = "Phase"
)

ProportionTestPlot(pancreas_sub)


# Plot specific comparisons
ProportionTestPlot(
  pancreas_sub,
  comparison = c("G2M_vs_G1", "G2M_vs_S")
)


# Plot paired comparisons using list format
ProportionTestPlot(
  pancreas_sub,
  comparison = list(c("G2M", "G1"))
)


ProportionTestPlot(
  pancreas_sub,
  cols.sig = "blue",
  comparison = list(c("G2M", "G1"))
)
```
