# Check and report the type of data in `Seurat` object

This function checks and returns a string indicating the type of data.
It checks for the presence of infinite values, negative values, and
whether the values are floats or integers.

## Usage

``` r
CheckDataType(object, ...)

# S3 method for class 'Seurat'
CheckDataType(object, layer = "data", assay = NULL, verbose = TRUE, ...)

# Default S3 method
CheckDataType(object, verbose = TRUE, ...)
```

## Arguments

- object:

  A `Seurat` object or a matrix.

- ...:

  The message to print.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A string indicating the type of data. Possible values are:
`"raw_counts"`, `"log_normalized_counts"`, `"raw_normalized_counts"`, or
`"unknown"`.

## Examples

``` r
data(pancreas_sub)
CheckDataType(pancreas_sub)
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-01-29 12:38:27] Infinite values detected
#> [1] "unknown"
```
