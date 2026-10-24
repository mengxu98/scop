# Cell-cycle scoring with native sparse module scoring

This follows Seurat's cell-cycle classification contract while routing
the two module scores through \[AddModuleScore()\] when its native
sparse contract is met. Unsupported scoring branches transparently use
Seurat.

## Usage

``` r
CellCycleScoring(
  object,
  s.features,
  g2m.features,
  ctrl = NULL,
  set.ident = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- s.features, g2m.features:

  Feature vectors for S and G2M phases.

- ctrl:

  Number of control features; defaults to the smaller feature set.

- set.ident:

  Whether to set cell identities to the inferred phase.

- ...:

  Arguments passed to \[AddModuleScore()\].

## Value

A Seurat object containing \`S.Score\`, \`G2M.Score\`, and \`Phase\`.
