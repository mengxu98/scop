# List conda environments

List conda environments

## Usage

``` r
ListEnv(conda = "auto")
```

## Arguments

- conda:

  The path to a `conda` executable. Use `"auto"` to allow `reticulate`
  to automatically find an appropriate `conda` binary. See **Finding
  Conda** and
  [`conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  for more details.

## Value

A data frame of conda environments.
