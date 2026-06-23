# Run Giotto HMRF spatial domains

Run Giotto HMRF spatial domains

## Usage

``` r
GiottoHMRF(
  x,
  spatial_genes = NULL,
  network_name = "Delaunay_full",
  k = 20,
  betas = c(0, 10, 20),
  hmrf_name = "scop_HMRF",
  params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- spatial_genes:

  Spatial genes used by Giotto HMRF.

- network_name:

  Name for the Giotto spatial network.

- k:

  Number of HMRF domains.

- betas:

  HMRF beta values.

- hmrf_name:

  Name for the HMRF result.

- params:

  Additional parameters passed to \`Giotto::doHMRF()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.
