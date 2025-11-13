# Run doublet-calling with scDblFinder

This function performs doublet-calling using the scDblFinder package on
a Seurat object.

## Usage

``` r
db_scDblFinder(srt, assay = "RNA", db_rate = ncol(srt)/1000 * 0.01, ...)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for doublet-calling. Default is
  `"RNA"`.

- db_rate:

  The expected doublet rate. Default is calculated as
  `ncol(srt) / 1000 * 0.01`.

- ...:

  Additional arguments to be passed to
  [`scDblFinder::scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.html).
