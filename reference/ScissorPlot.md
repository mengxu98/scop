# Plot Scissor results

Plot Scissor results

## Usage

``` r
ScissorPlot(
  srt,
  plot_type = c("umap", "heatmap", "bar", "upset", "rose", "ring", "pie", "dot"),
  reduction = NULL,
  prefix = "Scissor",
  group.by = NULL,
  split.by = NULL,
  features = NULL,
  nfeatures = 50,
  feature_method = c("variance", "status_diff", "coef_cor", "input_order"),
  tool_name = "Scissor",
  status = c("Scissor+", "Scissor-"),
  include.background = TRUE,
  upset_top_n = NULL,
  cells = NULL,
  layer = "data",
  assay = NULL,
  max_cells = 100,
  cell_order = NULL,
  exp_method = "zscore",
  stat_type = c("percent", "count"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "RdBu",
  group_palette = "Chinese",
  group_palcolor = NULL,
  cell_annotation = NULL,
  cell_annotation_palette = "Chinese",
  cell_annotation_palcolor = NULL,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object after
  [RunScissor](https://mengxu98.github.io/scop/reference/RunScissor.md).

- plot_type:

  Plot type. `"umap"` shows embedding panels, `"heatmap"` shows a
  `FeatureHeatmap`, and statistical views such as `"bar"` and `"upset"`
  are drawn with
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html).

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- prefix:

  Prefix used by
  [RunScissor](https://mengxu98.github.io/scop/reference/RunScissor.md).

- group.by:

  Optional metadata column shown together with Scissor status. For
  `"heatmap"`, the default is the Scissor status column. For statistical
  plots, it is passed to
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html),
  except that `"upset"` uses it to split Scissor status distributions by
  group.

- split.by:

  Metadata column to facet by.

- features:

  Features to plot.

- nfeatures:

  Number of features to show when `plot_type = "heatmap"` and
  `features = NULL`.

- feature_method:

  Method used to rank heatmap features when `features = NULL`.
  `"variance"` ranks genes by variance in selected cells,
  `"status_diff"` ranks by the largest mean-expression difference
  between Scissor status groups, `"coef_cor"` ranks by absolute
  correlation with Scissor coefficients, and `"input_order"` keeps the
  [RunScissor](https://mengxu98.github.io/scop/reference/RunScissor.md)
  input order.

- tool_name:

  Name of the `srt@tools` entry created by
  [RunScissor](https://mengxu98.github.io/scop/reference/RunScissor.md).

- status:

  Scissor status levels included in the heatmap.

- include.background:

  Whether to include background cells in heatmap and background status
  in upset plots.

- upset_top_n:

  Maximum number of `group.by` levels to show in `plot_type = "upset"`.
  The most frequent levels are kept. `NULL` keeps all.

- cells:

  Cell names to include.

- layer:

  Assay layer to use.

- assay:

  Assay to use. `NULL` uses the default assay.

- max_cells:

  An integer, maximum number of cells to sample per group.

- cell_order:

  A vector of cell names defining the order of cells.

- exp_method:

  Expression transform: `"zscore"`, `"raw"`, `"fc"`, `"log2fc"`, or
  `"log1p"`.

- stat_type:

  The type of statistic to compute for the plot. Can be one of
  `"percent"`, `"count"`, or `"value"`. Continuous value mode supports
  `plot_type = "bar"` and `plot_type = "dot"`.

- combine:

  Whether to combine UMAP panels or StatPlot panels.

- nrow, ncol, byrow:

  Layout parameters passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
  or
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html).

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- heatmap_palette:

  Palette used for CNV heatmap values.

- group_palette:

  Palette for cell types (groups) in Manhattan plot.

- group_palcolor:

  Custom colors for cell types (groups) in Manhattan plot.

- cell_annotation:

  Cell annotations. Palette length should match `cell_annotation`.

- cell_annotation_palette:

  Color palette for cell-type annotations.

- cell_annotation_palcolor:

  Custom colors for cell-type annotations.

- show_row_names:

  Whether to draw row/column names for the heatmap body.

- show_column_names:

  Whether to draw row/column names for the heatmap body.

- cluster_rows:

  Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.

- cluster_columns:

  Whether to cluster heatmap rows/columns. Defaults are both `FALSE`.

- theme_use, theme_args:

  Theme used by
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html).

- verbose:

  Whether to print messages.

- ...:

  Additional arguments passed to
  [CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
  [FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
  [FeatureHeatmap](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md),
  or
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html),
  depending on `plot_type`.

## Value

A ggplot/patchwork object for embedding and statistical plots, or a list
returned by
[FeatureHeatmap](https://mengxu98.github.io/scop/reference/FeatureHeatmap.md)
for `plot_type = "heatmap"`.

## See also

ScissorPlot

## Examples

``` r
data(panc8_sub)
data(islet_bulk)
panc8_sub <- RunScissor(
  panc8_sub,
  bulk_dataset = islet_bulk,
  condition.by = "condition",
  positive = "bfa",
  family = "binomial",
  features = head(intersect(
    rownames(panc8_sub),
    rownames(SummarizedExperiment::assay(islet_bulk, "counts"))
  ), 1000),
  alpha = 0.2,
  cutoff = 0.5
)
#> ℹ [2026-09-06 22:42:04] Build a temporary Scissor-style SNN graph
#> ℹ [2026-09-06 22:42:08] Scissor alpha 0.2 selected 1 positive and 501 negative cells (31.375%)
#> ✔ [2026-09-06 22:42:08] Scissor stored 1 Scissor+ and 501 Scissor- cells
panc8_sub <- RunStandardWorkflow(panc8_sub, verbose = FALSE)
#> ℹ [2026-09-06 22:42:10] Skip `log1p()` because `layer = data` is not "counts"

ScissorPlot(
  panc8_sub,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)


ScissorPlot(
  panc8_sub,
  group.by = "celltype",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)


ht <- ScissorPlot(
  panc8_sub,
  plot_type = "heatmap",
  group.by = "celltype"
)
ht$plot


ScissorPlot(
  panc8_sub,
  plot_type = "bar",
  group.by = "celltype"
)


ScissorPlot(
  panc8_sub,
  plot_type = "upset"
)
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?


ScissorPlot(
  panc8_sub,
  plot_type = "rose",
  label = TRUE
)


ScissorPlot(
  panc8_sub,
  plot_type = "ring",
  label = TRUE
)


ScissorPlot(
  panc8_sub,
  plot_type = "pie",
  label = TRUE
)


ScissorPlot(
  panc8_sub,
  plot_type = "dot",
  label = TRUE
)
```
