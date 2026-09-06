# Convert Giotto to Seurat

Convert a native Giotto object into a new Seurat object.

## Usage

``` r
giotto_to_srt(giotto, ...)
```

## Arguments

- giotto:

  A native Giotto object.

- ...:

  Additional arguments passed to the Giotto converter.

## Value

A new \`Seurat\` object.

## Details

The bridge selects the GiottoClass v4 or v5 converter according to the
installed Seurat major version. For round trips produced by
\[srt_to_giotto()\], normalize the Seurat input first (for example with
\[Seurat::NormalizeData()\]); the official GiottoClass converter records
an empty normalized layer otherwise, which its reverse converter cannot
read back.

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
  roundtrip <- giotto_to_srt(giotto)
  print(roundtrip)
}
} # }
```
