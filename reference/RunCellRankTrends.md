# Fit and cluster CellRank lineage trends

Refit CellRank GAM trends for a selected lineage and store the
normalized trend matrix and gene modules in
\`srt@tools\$CellRank\$trends\`.

## Usage

``` r
RunCellRankTrends(
  srt,
  lineage,
  features = NULL,
  top_n = 500L,
  heatmap_n = 50L,
  time_key = "palantir_pseudotime",
  assay = "RNA",
  layer = "data",
  n_points = 200L,
  norm = TRUE,
  resolution = 0.7,
  max_iter = 1000L,
  spline_order = 3L,
  n_knots = 8L,
  cores = 1L,
  random_state = 0L,
  output_dir = NULL,
  envname = NULL,
  conda = "auto",
  min_expressed_cells = 20L,
  distribution = "gamma",
  link = "log",
  fallback_distribution = "normal",
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object returned by \[RunCellRank\].

- lineage:

  A lineage name in the stored fate-probability matrix.

- features:

  Optional genes. If \`NULL\`, positive lineage drivers are used.

- top_n:

  Maximum number of driver genes used for clustering.

- heatmap_n:

  Number of genes retained for the compact heatmap view.

- time_key:

  Metadata column containing pseudotime.

- assay:

  Assay containing the expression layer.

- layer:

  Expression layer used for GAM fitting.

- n_points:

  Number of points used to predict each trend.

- norm:

  Whether to z-normalize each gene trend before clustering.

- resolution:

  Leiden resolution for trend modules.

- max_iter:

  Maximum GAM iterations.

- spline_order:

  Spline order passed to CellRank GAM.

- n_knots:

  Number of GAM knots.

- cores:

  Number of Python workers.

- random_state:

  Reproducibility seed.

- output_dir:

  Optional directory for CSV exports.

- envname:

  Optional Python environment name.

- conda:

  Conda-compatible executable used by \[PrepareEnv\].

- min_expressed_cells:

  Minimum cells with positive finite expression required for a GAM
  feature.

- distribution:

  Primary CellRank GAM distribution.

- link:

  Link function for the primary GAM distribution.

- fallback_distribution:

  Optional distribution used when the primary model has a fatal fit
  failure; the selected distribution is recorded.

- verbose:

  Whether to print progress messages.

## Value

The Seurat object with a \`CellRank\$trends\[\[lineage\]\]\` result
bundle.
