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
  target = NULL,
  measure = NULL
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

- measure:

  Numeric result column to display. For \`type = "improvements"\` use
  \`"gain.R2"\` or \`"gain.RMSE"\`; for \`type = "contributions"\`, use
  \`"contribution"\` or \`"importance"\` as provided by the backend.
  \`NULL\` defaults to \`"gain.R2"\` for improvements and
  \`"contribution"\` for contributions, and errors if that column is not
  present.

## Value

A \`ggplot\` object.

## See also

\[RunMistyR()\]
