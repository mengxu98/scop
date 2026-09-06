# Plot stored Statial Kontextual results

Plot contextual relationship scores across radii from a result produced
by \[RunStatialKontextual()\] without rerunning Statial.

## Usage

``` r
StatialKontextualPlot(object = NULL, res = NULL, tests = NULL, images = NULL)
```

## Arguments

- object:

  Optional \`Seurat\` object containing the result.

- res:

  Optional result list, usually \`object@tools\$StatialKontextual\`.

- tests:

  Optional relationship names to retain.

- images:

  Optional image identifiers to retain.

## Value

A \`ggplot\` object.

## Examples

``` r
data(visium_human_pancreas_results_sub)
StatialKontextualPlot(
  res = visium_human_pancreas_results_sub@tools$StatialKontextual
)
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
```
