# Run STdeconvolve reference-free spatial deconvolution

Estimate spot-level topic proportions from a spatial `Seurat` object
using the optional `STdeconvolve` package.

## Usage

``` r
RunSTdeconvolve(
  srt,
  assay = NULL,
  layer = "counts",
  features = NULL,
  k = NULL,
  k_candidates = 2:9,
  opt = "min",
  clean_counts = TRUE,
  clean_counts_params = list(),
  restrict_corpus = TRUE,
  restrict_corpus_params = list(),
  fit_lda_params = list(),
  get_beta_theta_params = list(),
  prefix = "STdeconvolve",
  tool_name = "STdeconvolve",
  store_results = TRUE,
  round_counts = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  Spatial `Seurat` object used as the RCTD query.

- assay:

  Assay used in `srt`. If `NULL`, the default assay is used.

- layer:

  Assay layer used as STdeconvolve input.

- features:

  Features used for RCTD. If `NULL`, shared features are used.

- k:

  Number of topics. If `NULL`, models are fit over `k_candidates` and
  [`STdeconvolve::optimalModel()`](https://rdrr.io/pkg/STdeconvolve/man/optimalModel.html)
  is used to choose a model.

- k_candidates:

  Candidate topic numbers used when `k = NULL`.

- opt:

  Optimal model selector passed to
  [`STdeconvolve::optimalModel()`](https://rdrr.io/pkg/STdeconvolve/man/optimalModel.html).

- clean_counts:

  Whether to call
  [`STdeconvolve::cleanCounts()`](https://rdrr.io/pkg/STdeconvolve/man/cleanCounts.html).

- clean_counts_params:

  Additional parameters passed to
  [`STdeconvolve::cleanCounts()`](https://rdrr.io/pkg/STdeconvolve/man/cleanCounts.html).

- restrict_corpus:

  Whether to call
  [`STdeconvolve::restrictCorpus()`](https://rdrr.io/pkg/STdeconvolve/man/restrictCorpus.html).

- restrict_corpus_params:

  Additional parameters passed to
  [`STdeconvolve::restrictCorpus()`](https://rdrr.io/pkg/STdeconvolve/man/restrictCorpus.html).

- fit_lda_params:

  Additional parameters passed to
  [`STdeconvolve::fitLDA()`](https://rdrr.io/pkg/STdeconvolve/man/fitLDA.html).

- get_beta_theta_params:

  Additional parameters passed to
  [`STdeconvolve::getBetaTheta()`](https://rdrr.io/pkg/STdeconvolve/man/getBetaTheta.html).

- prefix:

  Prefix for metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store detailed RCTD results in `srt@tools`.

- round_counts:

  Whether to round non-integer counts before model fitting.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional parameters passed to the RCTD run step.

## Value

A `Seurat` object with topic proportions in metadata and detailed
results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- RunSTdeconvolve(
  visium_human_pancreas_sub,
  assay = "Spatial",
  k = 5
)
SpatialSpotPlot(spatial, group.by = "STdeconvolve_dominant_type")
STdeconvolvePlot(spatial, plot_type = "pie")
} # }
```
