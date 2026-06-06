# Plot scFEA flux Cohen's d volcano plots

Plot scFEA flux Cohen's d volcano plots

## Usage

``` r
scFEAVolcanoPlot(
  srt,
  group.by,
  ident.1 = NULL,
  ident.2 = NULL,
  assay = "scFEAflux",
  layer = "data",
  label_by = c("reaction", "module", "module_reaction"),
  pathways = NULL,
  p_adj_cutoff = 0.001,
  cohen_cutoff = 0.2,
  combine = TRUE,
  width = 12,
  height = 10.4,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object returned by \[RunscFEA()\].

- group.by:

  Metadata column defining groups. If \`ident.1\` and \`ident.2\` are
  both \`NULL\`, each group is compared against all remaining cells.

- ident.1, ident.2:

  Group names to compare. Cohen's d is \`mean(ident.1) - mean(ident.2)\`
  divided by pooled standard deviation. If both are \`NULL\`,
  one-vs-rest contrasts are generated automatically.

- assay:

  Flux assay name.

- layer:

  Flux assay layer.

- label_by:

  Label mode for significant pathway modules.

- pathways:

  Optional \`SM_anno\` classes to plot.

- p_adj_cutoff:

  Adjusted p-value cutoff for highlighting, labels, and the horizontal
  threshold line.

- cohen_cutoff:

  Absolute Cohen's d cutoff for highlighting, labels, and the vertical
  threshold lines.

- combine:

  Whether to combine pathway plots into 2 x 2 paged panels. If
  \`FALSE\`, returns one plot per pathway.

- width, height:

  Suggested export size. Stored as attributes on the returned plot
  object.

- verbose:

  Whether to print messages.

## Value

For a single contrast, a paged list of \`ggplot\` objects when \`combine
= TRUE\`, otherwise a named list of pathway \`ggplot\` objects.
Statistics are stored in the \`"data"\` attribute. For automatic
one-vs-rest mode, a named list of single-contrast results is returned
and combined statistics are stored in the outer \`"data"\` attribute.
