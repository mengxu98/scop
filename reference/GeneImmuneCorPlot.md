# Gene-immune correlation butterfly plot

Draw a linkET-style correlation butterfly plot. The heatmap shows
immune-cell or immune-signature correlations, and curved links show
correlations between target genes and immune abundances.

## Usage

``` r
GeneImmuneCorPlot(
  object = NULL,
  gene.data = NULL,
  immune.data = NULL,
  features,
  immune.cols = NULL,
  assay = NULL,
  layer = "data",
  cor_method = c("spearman", "pearson", "kendall"),
  p_cutoff = 0.05,
  abs_cor_breaks = c(0.2, 0.4),
  heatmap_colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"),
  link_colors = c(Positive = "#E71D36", Negative = "#0073C2", Not = "grey82"),
  link_sizes = c(">= 0.4" = 1.5, "0.2 - 0.4" = 0.8, "< 0.2" = 0.3),
  title = NULL,
  subtitle = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE
)
```

## Arguments

- object:

  Optional `Seurat`, `SummarizedExperiment`, deconvolution bundle, or
  matrix-like object.

- gene.data:

  Optional gene expression matrix with genes in rows and samples in
  columns.

- immune.data:

  Optional immune abundance matrix with samples in rows and immune cell
  types or signatures in columns.

- features:

  Target genes to correlate with immune abundance.

- immune.cols:

  Metadata columns to extract from a `Seurat` object.

- assay, layer:

  Expression assay and layer settings.

- cor_method:

  Correlation method.

- p_cutoff:

  P-value cutoff used for positive/negative link categories.

- abs_cor_breaks:

  Two numeric breakpoints for link width categories.

- heatmap_colors:

  Continuous colors for immune-cell correlations.

- link_colors:

  Colors for `Positive`, `Negative`, and `Not` links.

- link_sizes:

  Link sizes for strong, medium, and weak absolute correlations.

- title, subtitle:

  Plot title and subtitle.

- theme_use, theme_args:

  Theme settings.

- verbose:

  Whether to print progress messages.

## Value

A `ggplot` object.
