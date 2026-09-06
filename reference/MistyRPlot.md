# Plot stored MISTy results

Plot model improvements or view contributions from a result produced by
\[RunMistyR()\] without rerunning the backend.

## Usage

``` r
MistyRPlot(
  object = NULL,
  res = NULL,
  type = c("improvements", "contributions"),
  top_n = 20,
  target = NULL
)
```

## Arguments

- object:

  Optional \`Seurat\` object containing \`MistyR\` results.

- res:

  Optional result list, usually \`object@tools\$MistyR\`.

- type:

  Result table to plot.

- top_n:

  Maximum number of records shown after ranking by absolute value.

- target:

  Optional target feature filter.

## Value

A \`ggplot\` object.

## Examples

``` r
data(visium_human_pancreas_results_sub)
MistyRPlot(
  res = visium_human_pancreas_results_sub@tools$MistyR,
  type = "improvements",
  top_n = 10
)
```
