# Find features with expression patterns similar to provided features

Find features with expression patterns similar to provided features

## Usage

``` r
GetSimilarFeatures(
  srt,
  features,
  n,
  features_use = rownames(srt),
  anticorr = FALSE,
  aggregator = "sum",
  assay = "RNA",
  layer = "data",
  verbose = TRUE
)
```

## Arguments

- srt:

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

## Value

character vector.
