# Immune abundance plots

Visualize immune abundance or score matrices from deconvolution results,
spatial metadata, GSVA-like scores, or user-provided matrices.

## Usage

``` r
ImmuneAbundancePlot(
  object = NULL,
  immune.data = NULL,
  immune.cols = NULL,
  group.by = NULL,
  group.data = NULL,
  plot_type = c("heatmap", "bar", "box", "cor"),
  sample_order = NULL,
  cell_type_order = NULL,
  scale = c("none", "row", "column"),
  cor_method = c("spearman", "pearson", "kendall"),
  bar_position = c("stack", "fill"),
  show_sample_names = FALSE,
  show_cor_label = FALSE,
  cor_label_size = 3,
  add_stat = FALSE,
  comparisons = NULL,
  pairwise_method = "wilcox.test",
  sig_label = c("p.signif", "p.format"),
  sig_labelsize = 3.5,
  box_width = 0.64,
  box_alpha = 0.92,
  pt.size = 1.35,
  pt.alpha = 0.72,
  jitter.width = 0.12,
  grid_major = FALSE,
  palette = "Chinese",
  palcolor = NULL,
  heatmap_palette = "YlGnBu",
  heatmap_palcolor = NULL,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list()
)
```

## Arguments

- object:

  Optional `SummarizedExperiment`, `Seurat`, deconvolution bundle, or
  matrix-like object.

- immune.data:

  Optional abundance matrix with samples/spots/cells in rows and immune
  cell types or signatures in columns.

- immune.cols:

  Metadata columns to extract from a `Seurat` object. If `NULL`, columns
  matching `RCTD_prop_*`, `*_prop_*`, or `*_frac_*` are used.

- group.by:

  Optional grouping column for `bar` and `box` plots.

- group.data:

  Optional named vector or data frame containing sample groups when
  groups cannot be inferred from `object`.

- plot_type:

  Plot type: `"heatmap"`, `"bar"`, `"box"`, or `"cor"`.

- sample_order:

  Optional sample order.

- cell_type_order:

  Optional immune cell type or signature order.

- scale:

  Scaling for heatmap values: `"none"`, `"row"`, or `"column"`.

- cor_method:

  Correlation method for `plot_type = "cor"`.

- bar_position:

  Bar position for `plot_type = "bar"`. `"stack"` shows raw abundance;
  `"fill"` rescales each sample to 1.

- show_sample_names:

  Whether to show sample names in heatmap and bar plots.

- show_cor_label:

  Whether to print correlation values in `cor` plots.

- cor_label_size:

  Text size for correlation labels.

- add_stat:

  Whether to add group comparison labels to `box` plots. Requires
  `ggpubr`.

- comparisons:

  Optional group comparisons passed to
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- pairwise_method:

  Statistical test for box plot comparisons.

- sig_label:

  Significance label type for
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- sig_labelsize:

  Significance label size.

- box_width, box_alpha:

  Box plot width and alpha.

- pt.size, pt.alpha, jitter.width:

  Point size, alpha, and jitter width.

- grid_major:

  Whether to show major y grid lines for bar and box plots.

- palette:

  Discrete palette for grouped plots.

- palcolor:

  Optional custom discrete colors.

- heatmap_palette:

  Continuous palette name for heatmaps.

- heatmap_palcolor:

  Optional custom continuous colors.

- title, subtitle:

  Plot title and subtitle.

- legend.position:

  Legend position.

- legend.direction:

  Legend direction.

- theme_use:

  Theme function name.

- theme_args:

  Additional theme arguments.

## Value

A `ggplot` object.
