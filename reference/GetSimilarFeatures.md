# Find genes with expression patterns similar to the genes you've specified.

Given a Seurat object and a list of feature names, this function returns
features that are strongly correlated with those markers.

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

  A character vector; giving feature names.

- n:

  An integer; number of results to return.

- features_use:

  A character vector; list of features eligible to be returned.

- anticorr:

  Whether to allow negatively correlated genes. Default is `FALSE`.

- aggregator:

  How to combine correlations when finding similar features. Options:
  `"sum"` (default), `"min"` (for "and"-like filter), `"max"`, or
  `"mean"`.

- assay:

  The assay to use. Default is `"RNA"`.

- layer:

  The layer to use. Default is `"data"`.

## Value

character vector.
