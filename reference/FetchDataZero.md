# FetchData but with zeroes for unavailable genes

FetchData but with zeroes for unavailable genes

## Usage

``` r
FetchDataZero(
  object,
  features,
  assay = "RNA",
  layer = "data",
  verbose = TRUE,
  srt = NULL,
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- features:

  Feature names.

- assay:

  Which assay to use. Default is `"RNA"`.

- layer:

  Assay layer to use.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

- ...:

  Arguments passed to other functions.
