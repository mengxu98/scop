# Run doublet-calling with scDblFinder

Run doublet-calling with scDblFinder

## Usage

``` r
RunscDblFinder(
  srt,
  assay = "RNA",
  db_rate = ncol(srt)/1000 * 0.01,
  data_type = NULL,
  ...,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- assay:

  Assay used for doublet calling.

- db_rate:

  Expected doublet rate.

- data_type:

  Optional
  [CheckDataType](https://mengxu98.github.io/scop/reference/CheckDataType.md)
  result, used internally to avoid rescanning the count matrix.

- ...:

  Additional arguments to be passed to
  [`scDblFinder::scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.html).

- verbose:

  Whether to print messages.
