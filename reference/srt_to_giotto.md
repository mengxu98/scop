# Convert Seurat to a native Giotto object

Convert one Seurat spatial image into a native Giotto object without
running a Giotto workflow or changing the input object.

## Usage

``` r
srt_to_giotto(object, image = NULL, ..., srt = NULL)
```

## Arguments

- object:

  A \`Seurat\` object.

- image:

  Seurat image name. Multi-image objects require an explicit name.

- ...:

  Additional arguments passed to the Giotto converter.

- srt:

  Deprecated alias for \`object\`; supply exactly one of the two. It
  will be removed in scop 1.0.0.

## Value

A native Giotto object.

## Details

Converters live in GiottoClass. If that package is already installed
from any source, the bridge reuses it and does not reinstall
\`drieslab/Giotto\`. A GitHub reinstall is requested only when
GiottoClass is missing. Conversion does not initialize Giotto's optional
Python/conda environment.

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
print(giotto)
} # }
```
