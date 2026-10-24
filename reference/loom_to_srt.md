# Read a `.loom` file and convert to a `Seurat`

Read a `.loom` file and convert to a `Seurat`

## Usage

``` r
loom_to_srt(
  path,
  layers = c("spliced", "unspliced"),
  verbose = TRUE,
  chunk_rows = 1000
)
```

## Arguments

- path:

  Path to a `.loom` file.

- layers:

  Character vector of loom layers to import as additional Seurat assays.
  Missing layers are skipped with a warning. Default is
  `c("spliced", "unspliced")`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- chunk_rows:

  Number of feature rows to read from each matrix dataset per chunk.
  Larger values can be faster but use more memory.

## Value

A `Seurat` object.

## Details

This function reads loom/HDF5 files directly in R using `rhdf5`; it does
not call `reticulate`, `scanpy`, or
[`loom_to_adata()`](https://mengxu98.github.io/scop/reference/loom_to_adata.md).
The loom `/matrix` dataset is imported as the `RNA` assay. Requested
`/layers/*` datasets are imported as additional assays, which is useful
for velocity-style loom files with `spliced` and `unspliced` layers.

## See also

[loom_to_adata](https://mengxu98.github.io/scop/reference/loom_to_adata.md),
[adata_to_srt](https://mengxu98.github.io/scop/reference/adata_to_srt.md),
[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
srt <- loom_to_srt("path/to/data.loom")
srt
} # }
```
