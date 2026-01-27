# Transcript-averaged cell scoring (TACS)

TACS is a method for plotting a FACS-like plot for two features based on
sc-RNA-seq data. For each of two query features, 100 features with
similar expression patterns are selected and ranked by their Pearson
correlation with the query. In a process akin to compensation, the
intersection of the feature lists is removed from each list. The log
normalized expression of the resulting features are then averaged within
each cell, and the resulting quantities are plotted. This function is
based on a simple scheme: choose features similar to the ones specified
and average them to reduce the noise.

## Usage

``` r
TACSPlot(
  srt,
  ref_srt = NULL,
  assay = "RNA",
  layer = "data",
  group.by = NULL,
  feature1,
  feature2,
  cutoffs = NULL,
  density = FALSE,
  palette = "Paired",
  num_features_add = 100,
  features_predetermined = FALSE,
  aggregator = "sum",
  remove_outliers = FALSE,
  aspect.ratio = 1,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  suffix = " expression level",
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  include_all = FALSE,
  all_color = "grey20",
  quadrants_line_color = "grey30",
  quadrants_line_type = "solid",
  quadrants_line_width = 0.3,
  quadrants_label_size = 3,
  density_alpha = NULL,
  bins = 20,
  h = NULL,
  nrow = NULL,
  ncol = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- ref_srt:

  A Seurat object. If your dataset is perturbed in a way that would
  substantially alter feature-feature correlations, for example if
  different time points are present or certain cell types are mostly
  depleted, you can feed in a reference srt, and TACS will choose axes
  based on the reference data. Default is `NULL`.

- assay:

  Which assay to use. Default is `"RNA"`.

- layer:

  Which layer to use. Default is `data`.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- feature1:

  Horizontal axis on plot mimics this feature. Character, usually length
  1 but possibly longer.

- feature2:

  Vertical axis on plot mimics this feature. Character, usually length 1
  but possibly longer.

- cutoffs:

  If given, divide plot into four quadrants and annotate with
  percentages. Can be a numeric vector of length 1 or 2, or a list of
  two numeric vectors for x and y axes respectively.

- density:

  If `TRUE`, plot contours instead of points.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Paired"`.

- num_features_add:

  Each axis shows a simple sum of similar features. This is how many
  (before removing overlap).

- features_predetermined:

  If `FALSE`, plot the sum of many features similar to feature1 instead
  of feature1 alone (same for feature2). See
  [GetSimilarFeatures](https://mengxu98.github.io/scop/reference/GetSimilarFeatures.md).
  If `TRUE`, plot the sum of only the features given.

- aggregator:

  How to combine correlations when finding similar features. Options:
  `"sum"` (default), `"min"` (for "and"-like filter), `"max"`, or
  `"mean"`.

- remove_outliers:

  If `TRUE`, remove outliers from the plot. Default is `FALSE`.

- aspect.ratio:

  Aspect ratio of the panel. Default is `1`.

- title:

  The text for the title. Default is `NULL`.

- subtitle:

  The text for the subtitle for the plot which will be displayed below
  the title. Default is `NULL`.

- xlab:

  The x-axis label of the plot. Default is `NULL`.

- ylab:

  The y-axis label of the plot. Default is `NULL`.

- suffix:

  The suffix of the axis labels. Default is `" expression level"`.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

- legend.direction:

  The direction of the legend in the plot. Can be one of `"vertical"` or
  `"horizontal"`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- include_all:

  If `TRUE`, include a panel with all cells. Default is `FALSE`.

- all_color:

  The color of the all cells panel. Default is `"grey20"`.

- quadrants_line_color:

  The color of the quadrants lines. Default is `"grey30"`.

- quadrants_line_type:

  The type of the quadrants lines. Default is `"solid"`.

- quadrants_line_width:

  The width of the quadrants lines. Default is `0.3`.

- quadrants_label_size:

  The size of the quadrants labels. Default is `3`.

- density_alpha:

  The alpha of the density plot. Default is `NULL`.

- bins:

  Number of bins for density plot. Default is `20`.

- h:

  Bandwidth for density plot. Default is `NULL`.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to
  [ggplot2::stat_density2d](https://ggplot2.tidyverse.org/reference/geom_density_2d.html).

## References

[Kernfeld et al. paper](https://doi.org/10.1016/j.immuni.2018.04.015),
[Github](https://github.com/maehrlab/thymusatlastools2/blob/f8b51ad684d56b2eeda780787eb9ad4ff3003eef/R/data_handling_seurat.R#L271)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-27 08:35:40] Start standard scop workflow...
#> ℹ [2026-01-27 08:35:41] Checking a list of <Seurat>...
#> ! [2026-01-27 08:35:41] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-27 08:35:41] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 08:35:43] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 08:35:43] Use the separate HVF from srt_list
#> ℹ [2026-01-27 08:35:43] Number of available HVF: 2000
#> ℹ [2026-01-27 08:35:44] Finished check
#> ℹ [2026-01-27 08:35:44] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-27 08:35:44] Perform pca linear dimension reduction
#> ℹ [2026-01-27 08:35:45] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-27 08:35:45] Reorder clusters...
#> ℹ [2026-01-27 08:35:45] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-27 08:35:45] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-27 08:35:50] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-27 08:35:54] Run scop standard workflow completed
TACSPlot(
  pancreas_sub,
  feature1 = "H3f3b",
  feature2 = "Eif1",
  group.by = "CellType"
)


TACSPlot(
  pancreas_sub,
  feature1 = "H3f3b",
  feature2 = "Eif1",
  group.by = "CellType",
  density = TRUE,
  include_all = TRUE,
  cutoffs = c(3, 2.5)
)


TACSPlot(
  pancreas_sub,
  feature1 = "H3f3b",
  feature2 = "Eif1",
  group.by = "CellType",
  density = TRUE,
  cutoffs = list(x = c(2, 3), y = c(2.5))
)


TACSPlot(
  pancreas_sub,
  feature1 = "H3f3b",
  feature2 = "Eif1",
  group.by = "SubCellType",
  density = TRUE
)
```
