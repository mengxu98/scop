# Access stored SpatialDM results

Access stored SpatialDM results

## Usage

``` r
GetSpatialDMResult(
  object,
  result.name = NULL,
  type = c("global", "local", "weights"),
  pair = NULL
)
```

## Arguments

- object:

  A Seurat object with SpatialDM results.

- result.name:

  Stored result name.

- type:

  Result type to return (\`"global"\`, \`"local"\`, or \`"weights"\`).

- pair:

  Optional LR interaction name for local results.

## Value

The requested stored SpatialDM result payload.
