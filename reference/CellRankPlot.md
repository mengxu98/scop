# Plot stored CellRank results

Plot CellRank outputs without rerunning a Python backend.

## Usage

``` r
CellRankPlot(
  object,
  plot_type = c("fate", "states", "circular", "drivers", "trends", "clusters",
    "enrichment", "projection", "random_walks"),
  lineage = NULL,
  database = NULL,
  reduction = NULL,
  group.by = NULL,
  top_n = 20L,
  n_sims = 100L,
  max_iter = 500L,
  seed = 0L,
  palette = "Chinese",
  palcolor = NULL,
  feature_palette = "Spectral",
  feature_palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A Seurat object returned by \[RunCellRank\].

- plot_type:

  One of \`"fate"\`, \`"states"\`, \`"circular"\`, \`"drivers"\`,
  \`"trends"\`, \`"clusters"\`, \`"enrichment"\`, \`"projection"\`, or
  \`"random_walks"\`.

- lineage:

  Lineage used for driver/trend plots.

- database:

  Enrichment database used when \`plot_type = "enrichment"\`. If
  \`NULL\`, the first stored non-empty database is used.

- reduction:

  Reduction used for cell-space plots.

- group.by:

  Group column used for state/random-walk plots.

- top_n:

  Number of driver genes to display.

- n_sims:

  Number of random walks.

- max_iter:

  Maximum length of each random walk.

- seed:

  Random seed.

- palette, palcolor:

  Discrete SCOP palette for states, modules, and cell groups.

- feature_palette, feature_palcolor:

  Continuous SCOP palette for fate, trends, pseudotime, and enrichment
  strength.

- theme_use, theme_args:

  SCOP plot theme and arguments passed to it.

- ...:

  Arguments passed to the underlying SCOP plotting function.

- srt:

  Deprecated alias for \`object\`; supply exactly one of the two. It
  will be removed in scop 1.0.0.

## Value

A ggplot object or a SCOP plot object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunCellRank(
  pancreas_sub,
  group.by = "SubCellType",
  kernel_type = "pseudotime",
  time_key = "palantir_pseudotime",
  schur_method = "brandts",
  backend = "python"
)
CellRankPlot(pancreas_sub, plot_type = "fate", reduction = "UMAP")
CellRankPlot(pancreas_sub, plot_type = "circular")
} # }
```
