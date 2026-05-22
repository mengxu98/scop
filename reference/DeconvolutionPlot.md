# Plot deconvolution results

Visualize bulk or pseudobulk deconvolution results from a
`SummarizedExperiment` object or a supplied result table.

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

- theme_use:

  Theme function name. Default is `"theme_scop"`.

- theme_args:

  Additional theme arguments passed to `theme_use`.

- legend.position:

  Legend position. Default is `"right"`.

- legend.direction:

  Legend direction. Default is `"vertical"`.

- grid_major:

  Whether to show major grid lines for `"bar"` and `"box"` plots.
  Default is `FALSE`.

- ...:

  Reserved for future use.

## Value

A `ggplot` object. For `plot_type = "heatmap"`, returns a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object.

## See also

[RunDeconvolution](https://mengxu98.github.io/scop/reference/RunDeconvolution.md)

## Examples

``` r
data(islet_bulk)
data(panc8_sub)
islet_bulk <- RunDeconvolution(
  islet_bulk,
  reference = panc8_sub,
  method = "MuSiC",
  group.by = "celltype"
)
#> Creating Relative Abudance Matrix...
#> Creating Variance Matrix...
#> Creating Library Size Matrix...
#> Used 12276 common genes...
#> Used 12 cell types in deconvolution...
#> HI-070719-DMSO_S28 has common genes 12192 ...
#> HI-070719-BFA_S29 has common genes 12028 ...
#> HI-043019-BFA_S25 has common genes 11854 ...
#> HI-043019-DMSO_S24 has common genes 12062 ...
#> HI-061119-BFA_S27 has common genes 11957 ...
#> HI-061119-DMSO_S26 has common genes 12093 ...
#> HI-042519-BFA_S23 has common genes 12029 ...
#> HI-042519-DMSO_S22 has common genes 12103 ...
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
