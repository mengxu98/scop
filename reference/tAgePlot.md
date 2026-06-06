# Plot tAge transcriptomic aging-clock predictions

Plot tAge transcriptomic aging-clock predictions

## Usage

``` r
tAgePlot(
  object,
  score_col = NULL,
  group.by = NULL,
  split.by = NULL,
  tool_name = "tAge",
  plot_type = c("box"),
  palette = "Chinese",
  palcolor = NULL,
  alpha = 1,
  point_size = 2.2,
  point_alpha = 0.85,
  box_alpha = 0.15,
  flip = TRUE,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  legend.position = "none",
  legend.direction = "vertical",
  grid_major = TRUE,
  grid_major_colour = "grey80",
  grid_major_linetype = 2,
  grid_major_linewidth = 0.3,
  aspect.ratio = NULL,
  seed = 11,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `Seurat` object returned by
  [`RuntAge()`](https://mengxu98.github.io/scop/reference/RuntAge.md) or
  a data frame with tAge prediction columns.

- score_col:

  Prediction column to plot. If `NULL`, the first column ending in
  `"_tAge"` is used.

- group.by:

  Column used to group predictions on the y axis.

- split.by:

  An optional column to split the plot into facets.

- tool_name:

  Name of the Seurat tool entry containing tAge results.

- plot_type:

  Plot type. Currently `"box"` draws a boxplot with jittered pseudobulk
  samples.

- palette, palcolor:

  Palette forwarded to
  [`thisplot::palette_colors()`](https://mengxu98.github.io/thisplot/reference/palette_colors.html).

- alpha:

  Overall point and box alpha. `point_alpha` and `box_alpha` take
  precedence when set.

- point_size, point_alpha:

  Jittered point size and alpha.

- box_alpha:

  Boxplot fill alpha.

- flip:

  Whether to flip coordinates so group labels run along the y axis.
  Default is `TRUE`.

- title, subtitle, xlab, ylab:

  Plot labels.

- theme_use:

  Theme function name. Default is `"theme_scop"`, which maps to
  [`thisplot::theme_this()`](https://mengxu98.github.io/thisplot/reference/theme_this.html),
  matching the style used by
  [`thisplot::StatPlot()`](https://mengxu98.github.io/thisplot/reference/StatPlot.html).

- theme_args:

  Additional arguments passed to `theme_use`.

- legend.position, legend.direction:

  Legend position and direction. Default is `"none"`.

- grid_major:

  Whether to show major panel grid lines.

- grid_major_colour, grid_major_linetype, grid_major_linewidth:

  Appearance of major panel grid lines.

- aspect.ratio:

  Aspect ratio of the plot (`y / x`).

- seed:

  RNG seed for jitter reproducibility.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Reserved for future use.

## Value

A `ggplot` object.

## See also

[RuntAge](https://mengxu98.github.io/scop/reference/RuntAge.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pancreas_sub <- RuntAge(pancreas_sub, group.by = "CellType")
tAgePlot(pancreas_sub, group.by = "CellType")
} # }
```
