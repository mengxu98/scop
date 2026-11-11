# Convert Giotto to Seurat

Convert a Giotto object into a new Seurat object. Requires GiottoClass.

## Usage

``` r
giotto_to_srt(giotto, ...)
```

## Arguments

- giotto:

  A Giotto object. For a round trip, normalize the Seurat object with
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
  before converting it to Giotto.

- ...:

  Arguments passed to other functions.

## Value

A new `Seurat` object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
giotto <- srt_to_giotto(spatial, image = "slice1")
giotto_to_srt(giotto)
} # }
```
