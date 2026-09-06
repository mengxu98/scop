# Run STdeconvolve reference-free spatial deconvolution

Estimate spot-level topic proportions from a spatial `Seurat` object
using the optional `STdeconvolve` package. The producer example is a
non-executing template;
[`STdeconvolvePlot()`](https://mengxu98.github.io/scop/reference/STdeconvolvePlot.md)
demonstrates a validated stored result without rerunning the backend.

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
  `STdeconvolve::optimalModel()` is used to choose a model.

- k_candidates:

  Candidate topic numbers used when `k = NULL`.

- opt:

  Optimal model selector passed to `STdeconvolve::optimalModel()`.

- clean_counts:

  Whether to call `STdeconvolve::cleanCounts()`.

- clean_counts_params:

  Additional parameters passed to `STdeconvolve::cleanCounts()`.

- restrict_corpus:

  Whether to call `STdeconvolve::restrictCorpus()`.

- restrict_corpus_params:

  Additional parameters passed to `STdeconvolve::restrictCorpus()`.

- fit_lda_params:

  Additional parameters passed to `STdeconvolve::fitLDA()`.

- get_beta_theta_params:

  Additional parameters passed to `STdeconvolve::getBetaTheta()`.

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
thisutils::check_r("JEFworks-Lab/STdeconvolve", verbose = FALSE)
data(visium_human_pancreas_sub)
keep_spots <- unique(round(seq(
  1,
  ncol(visium_human_pancreas_sub),
  length.out = 120
)))
spatial <- visium_human_pancreas_sub[, keep_spots]
spatial <- RunSTdeconvolve(
  spatial,
  assay = "Spatial",
  layer = "counts",
  features = head(rownames(spatial), 300),
  k = 3,
  verbose = FALSE
)
STdeconvolvePlot(
  spatial,
  topics = 1L,
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
} # }
```
