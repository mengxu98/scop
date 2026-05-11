# Read an `.h5ad` file and convert to a `Seurat`

Read an `.h5ad` file and convert to a `Seurat`

## Usage

``` r
h5ad_to_srt(path, verbose = TRUE, prepare_for_reticulate = TRUE)
```

## Arguments

- path:

  Path to an `.h5ad` file (passed to `anndata.read_h5ad()`).

- verbose:

  Whether to print the message. Default is `TRUE`.

- prepare_for_reticulate:

  If `TRUE` (default), coerces `X` and each layer matrix to CSR
  `float64` in Python (avoids invalid `dgRMatrix` conversion via
  reticulate). Layers that still fail in
  [`adata_to_srt()`](https://mengxu98.github.io/scop/reference/adata_to_srt.md)
  are skipped and reported. Set to `FALSE` for a plain `read_h5ad` then
  convert.

## Value

A `Seurat` object.

## See also

[adata_to_srt](https://mengxu98.github.io/scop/reference/adata_to_srt.md),
[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
srt <- h5ad_to_srt("path/to/data.h5ad")
srt
} # }
```
