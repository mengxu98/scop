# Find features with expression patterns similar to provided features

Find features with expression patterns similar to provided features

## Usage

``` r
GetSimilarFeatures(
  object,
  features,
  n,
  features_use = rownames(srt),
  anticorr = FALSE,
  aggregator = "sum",
  assay = "RNA",
  layer = "data",
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- features:

  Feature names.

- n:

  An integer; number of results to return.

- features_use:

  Features eligible to be returned.

- anticorr:

  Whether to allow negatively correlated features.

- aggregator:

  How to combine correlations when finding similar features. Options:
  `"sum"` (default), `"min"` (for "and"-like filter), `"max"`, or
  `"mean"`.

- assay:

  Which assay to use. Default is `"RNA"`.

- layer:

  Assay layer to use.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

character vector.
