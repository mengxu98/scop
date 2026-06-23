# Preprocess an internal Giotto workflow object

Preprocess an internal Giotto workflow object

## Usage

``` r
GiottoPreprocess(
  x,
  filter_params = list(),
  norm_params = list(),
  stat_params = list(),
  hvf_params = list(),
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- x:

  A \`giotto2\` workflow object.

- filter_params:

  Additional parameters reserved for future filtering.

- norm_params:

  Additional parameters passed to \`Giotto::normalizeGiotto()\`.

- stat_params:

  Additional parameters passed to \`Giotto::addStatistics()\`.

- hvf_params:

  Additional parameters passed to \`Giotto::calculateHVF()\`.

- verbose:

  Whether to print progress messages.

- seed:

  Random seed for reproducible Giotto calls.

## Value

A \`giotto2\` workflow object.
