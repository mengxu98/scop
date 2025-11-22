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

  The layer in the `srt` object from which to extract the data. Default
  is `"data"`.

- assay:

  The assay to extract the data from. If not provided, the default assay
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
#> [1] "unknown"
```
