# Plot stored COMMOT results

Plot stored COMMOT group communication, cluster matrices, or
raw-coordinate direction vectors without rerunning Python.

## Usage

``` r
COMMOTPlot(
  object,
  result.name = NULL,
  plot_type = c("network", "matrix", "direction"),
  key = NULL,
  ...
)
```

## Arguments

- object:

  A \`Seurat\` object with COMMOT results.

- result.name:

  Stored COMMOT result name.

- plot_type:

  Network, cluster matrix, or direction-vector view.

- key:

  Stored cluster or direction selection key.

- ...:

  Arguments passed to the SCOP network plot for \`"network"\`.

## Value

A \`ggplot\` or compatible plot object.
