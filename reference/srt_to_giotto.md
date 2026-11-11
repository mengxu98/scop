# Convert Seurat to a Giotto object

Convert one Seurat spatial image into a Giotto object without running a
Giotto workflow or changing the input object. Requires GiottoClass.

## Usage

``` r
srt_to_giotto(object, image = NULL, srt = NULL, ...)
```

## Arguments

- object:

  A \`Seurat\` object.

- image:

  Seurat image name. Multi-image objects require an explicit name.

- srt:

  Deprecated alias for \`object\`; supply exactly one of the two. It
  will be removed in scop 1.0.0.

- ...:

  Arguments passed to other functions.

## Value

A Giotto object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
spatial <- Seurat::NormalizeData(
  visium_human_pancreas_sub,
  assay = "Spatial",
  verbose = FALSE
)
srt_to_giotto(spatial, image = "slice1")
} # }
```
