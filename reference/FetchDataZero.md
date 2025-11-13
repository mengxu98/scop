# FetchData but with zeroes for unavailable genes

FetchData but with zeroes for unavailable genes

## Usage

``` r
FetchDataZero(
  srt,
  assay = "RNA",
  layer = "data",
  features,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The assay to use. Default is `"RNA"`.

- layer:

  The layer to use. Default is `"data"`.

- features:

  A character vector of feature names.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Other arguments to pass to FetchData
