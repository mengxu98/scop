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
  layer = "data"
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A character vector of feature names.

- n:

  An integer; number of results to return.

- features_use:

  A character vector of features eligible to be returned.

- anticorr:

  Whether to allow negatively correlated features. Default is `FALSE`.

- aggregator:

  How to combine correlations when finding similar features. Options:
  `"sum"` (default), `"min"` (for "and"-like filter), `"max"`, or
  `"mean"`.

- assay:

  Which assay to use. Default is `"RNA"`.

- layer:

  Which layer to use. Default is `data`.

## Value

character vector.
