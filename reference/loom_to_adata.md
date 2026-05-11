# Read a `.loom` file as an AnnData object

Read a `.loom` file as an AnnData object

## Usage

``` r
loom_to_adata(path, verbose = TRUE, ...)
```

## Arguments

- path:

  Path to a `.loom` file (passed to `scanpy.read_loom()`).

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to `scanpy.read_loom()`.

## Value

A Python `anndata.AnnData` object.

## Details

This is a Python-backed wrapper and requires `reticulate` plus a Python
environment with `scanpy` and `loompy` available. It is independent from
[`loom_to_srt()`](https://mengxu98.github.io/scop/reference/loom_to_srt.md),
which reads loom files directly in R without initializing Python.

## See also

[loom_to_srt](https://mengxu98.github.io/scop/reference/loom_to_srt.md),
[adata_to_srt](https://mengxu98.github.io/scop/reference/adata_to_srt.md),
[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
adata <- loom_to_adata("path/to/data.loom")
adata
} # }
```
