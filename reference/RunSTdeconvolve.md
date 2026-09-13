# Run STdeconvolve reference-free spatial deconvolution

Infer expression topics and their spot-level proportions using
STdeconvolve.

## Usage

``` r
RunSTdeconvolve(
  object,
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
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object containing spatial expression data.

- assay:

  Assay used as STdeconvolve input. If `NULL`, the default assay is
  used.

- layer:

  Assay layer used as STdeconvolve input.

- features:

  Features passed to the topic model. If `NULL`, all features in the
  selected assay are used.

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

  Prefix for topic proportion metadata columns.

- tool_name:

  Name used to store detailed results in `srt@tools`.

- store_results:

  Whether to store detailed topic matrices in `srt@tools`.

- round_counts:

  Whether to round non-integer counts before model fitting. This is a
  preprocessing choice and does not recover original counts.

- verbose:

  Whether to print messages.

- ...:

  Additional arguments passed to the STdeconvolve fit steps.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object with topic proportions in metadata and detailed
results stored in `srt@tools[[tool_name]]` when `store_results = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- visium_human_pancreas_sub
counts <- GetAssayData5(spatial, assay = "Spatial", layer = "counts")
features <- names(sort(Matrix::rowSums(counts), decreasing = TRUE))[
  seq_len(min(300L, nrow(counts)))
]
spatial <- RunSTdeconvolve(
  spatial,
  assay = "Spatial",
  layer = "counts",
  features = features,
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
