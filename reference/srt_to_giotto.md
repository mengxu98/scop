# Convert Seurat to a native Giotto object

Convert one Seurat spatial image into a native Giotto object without
running a Giotto workflow or changing the input object.

## Usage

``` r
srt_to_giotto(srt, image = NULL, ...)
```

## Arguments

- srt:

  A \`Seurat\` object.

- image:

  Seurat image name. Multi-image objects require an explicit name.

- ...:

  Additional arguments passed to the Giotto converter.

## Value

A native Giotto object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
if (all(unlist(
  thisutils::check_r("drieslab/Giotto", verbose = FALSE),
  use.names = FALSE
))) {
  giotto <- srt_to_giotto(spatial, image = "slice1")
  print(giotto)
}
} # }
```
