# Plot deconvolution results

Plot deconvolution results

## Usage

``` r
DeconvolutionPlot(
  object = NULL,
  res = NULL,
  plot_type = c("bar", "heatmap", "box"),
  sample_order = NULL,
  cell_type_order = NULL,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "YlGnBu",
  heatmap_palcolor = NULL,
  sample_annotation = NULL,
  sample_annotation_palette = "Chinese",
  sample_annotation_palcolor = NULL,
  sample_split = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  theme_use = "theme_scop",
  theme_args = list(),
  legend.position = "right",
  legend.direction = "vertical",
  grid_major = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `SummarizedExperiment` object containing deconvolution results in
  `metadata(object)[["Deconvolution"]]`.

- res:

  A deconvolution result data frame. When provided, `object` is ignored.

- plot_type:

  Plot type. One of `"bar"`, `"heatmap"`, or `"box"`.

- sample_order:

  Optional sample order.

- cell_type_order:

  Optional cell-type order.

- palette:

  Palette used for discrete cell-type colors in `"bar"` and `"box"`
  modes.

- palcolor:

  Optional custom colors for `palette`.

- heatmap_palette:

  Palette used for continuous proportions in `"heatmap"` mode.

- heatmap_palcolor:

  Optional custom colors for `heatmap_palette`.

- sample_annotation:

  Character vector of `colData(object)` columns used as top annotations
  in `"heatmap"` mode.

- sample_annotation_palette:

  Palette(s) used for `sample_annotation`.

- sample_annotation_palcolor:

  Optional custom colors for `sample_annotation_palette`. Use a list
  when multiple annotations are provided.

- sample_split:

  Optional `colData(object)` column used to split heatmap columns for
  grouped comparison.

- cluster_rows, cluster_columns:

  Whether to cluster rows or columns in `"heatmap"` mode.

- show_row_names, show_column_names:

  Whether to show row or column names in `"heatmap"` mode.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- grid_major:

  Whether to show major grid lines for `"bar"` and `"box"` plots.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Arguments passed to other functions.

## Value

A `ggplot` object. For `plot_type = "heatmap"`, returns a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object.

## Details

Bar plots use
[ImmuneAbundancePlot](https://mengxu98.github.io/scop/reference/ImmuneAbundancePlot.md)
with stacked stored proportions. Box plots use the shared
expression-statistics renderer. Neither view estimates cell-type
proportions again. The constructed result table below demonstrates
plotting previously computed proportions.

## See also

[RunDeconvolution](https://mengxu98.github.io/scop/reference/RunDeconvolution.md),
[ImmuneAbundancePlot](https://mengxu98.github.io/scop/reference/ImmuneAbundancePlot.md),
[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)

## Examples

``` r
proportions <- data.frame(
  sample = rep(paste0("Sample", 1:4), each = 2),
  cell_type = rep(c("Alpha", "Beta"), 4),
  proportion = c(0.7, 0.3, 0.6, 0.4, 0.4, 0.6, 0.3, 0.7),
  method = "CIBERSORT"
)
DeconvolutionPlot(res = proportions, plot_type = "bar")

DeconvolutionPlot(res = proportions, plot_type = "box")


data(islet_bulk)
islet_bulk <- RunDeconvolution(
  islet_bulk,
  method = "CIBERSORT",
  backend = "cpp",
  perm = 0
)
#> ℹ [2026-11-11 04:50:02] Use 400 shared genes for CIBERSORT
DeconvolutionPlot(islet_bulk, plot_type = "bar")


DeconvolutionPlot(islet_bulk, plot_type = "box")


ht <- DeconvolutionPlot(
  islet_bulk,
  plot_type = "heatmap",
  sample_annotation = "condition",
  sample_split = "condition"
)
ComplexHeatmap::draw(ht)
```
