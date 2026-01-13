# FetchData but with zeroes for unavailable genes

FetchData but with zeroes for unavailable genes

## Usage

``` r
FetchDataZero(
  srt,
  features,
  assay = "RNA",
  layer = "data",
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A character vector of feature names.

- assay:

  Which assay to use. Default is `"RNA"`.

- layer:

  Which layer to use. Default is `data`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Other arguments to pass to
  [Seurat::FetchData](https://satijalab.org/seurat/reference/reexports.html).
