# Benchmark single-cell integration methods

Run one or more
[`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md)
methods on the same input, score each latent embedding with scIB-style
batch-correction and bio-conservation metrics (scaled iLISI/cLISI, ASW,
graph connectivity, ARI/NMI), and store the tables in
`srt@tools$IntegrationBenchmark`. The object itself is returned so later
plotting can reuse the embeddings and per-cell LISI scores.

## Usage

``` r
RunIntegrationBenchmark(
  srt,
  batch,
  celltype = NULL,
  methods = c("Uncorrected", "Harmony"),
  linear_reduction = "pca",
  nHVF = 2000,
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  nonlinear_reduction = "umap",
  perplexity = 30,
  k_graph = 15,
  bio_weight = 0.6,
  batch_weight = 0.4,
  skip_failed = TRUE,
  tool_name = "IntegrationBenchmark",
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object containing a batch column.

- batch:

  Metadata column used as the technical batch.

- celltype:

  Optional metadata column used as biological labels.

- methods:

  Integration methods passed to
  [`RunIntegration()`](https://mengxu98.github.io/scop/reference/RunIntegration.md).

- perplexity:

  Neighborhood size for
  [`RunLISI()`](https://mengxu98.github.io/scop/reference/RunLISI.md).

- k_graph:

  Neighbor size for graph connectivity.

- bio_weight, batch_weight:

  Weights used for the overall score. Defaults follow scIB (`0.6` bio
  conservation, `0.4` batch correction).

- skip_failed:

  Whether a failed method is recorded and skipped.

- tool_name:

  Name of the `srt@tools` entry used to store the tables.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  The message to print.

## Value

The input `Seurat` object with integration reductions, per-cell LISI
columns, and `srt@tools[[tool_name]]` containing `summary`, `metrics`,
and `runs`. Print the tables with
[`thisplot::print_colored_table()`](https://mengxu98.github.io/thisplot/reference/print_colored_table.html).

## See also

[IntegrationBenchmarkPlot](https://mengxu98.github.io/scop/reference/IntegrationBenchmarkPlot.md),
[RunIntegration](https://mengxu98.github.io/scop/reference/RunIntegration.md),
[RunLISI](https://mengxu98.github.io/scop/reference/RunLISI.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(panc8_sub)
panc8_sub <- RunIntegrationBenchmark(
  panc8_sub,
  batch = "tech",
  celltype = "celltype",
  methods = c("Uncorrected", "Harmony"),
  nHVF = 500,
  linear_reduction_dims = 20,
  linear_reduction_dims_use = 1:10,
  perplexity = 10
)
thisplot::print_colored_table(
  panc8_sub@tools$IntegrationBenchmark$summary,
  by = "row",
  palette = "Chinese"
)
thisplot::print_colored_table(
  panc8_sub@tools$IntegrationBenchmark$metrics,
  by = "col",
  palette = "Paired"
)
IntegrationBenchmarkPlot(panc8_sub)
IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")
} # }
```
