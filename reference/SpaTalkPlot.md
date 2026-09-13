# Plot stored SpaTalk results

Plot a stored SpaTalk communication result without rerunning the
optional backend.

## Usage

``` r
SpaTalkPlot(
  object,
  result.name = NULL,
  plot_type = c("network", "bubble", "pathway", "tf"),
  theme_use = "theme_scop",
  theme_args = list(),
  ...
)
```

## Arguments

- object:

  A \`Seurat\` object with SpaTalk results.

- result.name:

  Stored SpaTalk result name.

- plot_type:

  Network, ligand-receptor bubble, pathway bubble, or receptor-TF view.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- ...:

  Arguments passed to the corresponding SCOP CCC plot.

## Value

A \`ggplot\` or compatible plot object.
