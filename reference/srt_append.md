# Append a Seurat object to another

Append a Seurat object to another

## Usage

``` r
srt_append(
  srt_raw,
  srt_append,
  slots = methods::slotNames(srt_append),
  pattern = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt_raw:

  A Seurat object to be appended.

- srt_append:

  New Seurat object to append.

- slots:

  slots names.

- pattern:

  A character string containing a regular expression. All data with
  matching names will be considered for appending.

- overwrite:

  Whether to overwrite.

- verbose:

  Whether to print the message. Default is `TRUE`.
