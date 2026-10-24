# Get an upstream CellChat-family object

Get an upstream CellChat-family object

## Usage

``` r
GetCCCObject(
  object,
  method = c("CellChat", "SpatialCellChat"),
  result.name = NULL,
  sample = NULL
)
```

## Arguments

- object:

  A \`Seurat\` object.

- method:

  Either \`"CellChat"\` or \`"SpatialCellChat"\`.

- result.name:

  A CellChat condition or SpatialCellChat named result.

- sample:

  Spatial sample for SpatialCellChat results.

## Value

An upstream backend object.
