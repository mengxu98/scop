# Print environment information

Print environment information

## Usage

``` r
env_info(conda, envname)
```

## Arguments

- conda:

  The path to a `conda` executable. Use `"auto"` to allow `reticulate`
  to automatically find an appropriate `conda` binary. See **Finding
  Conda** and
  [`conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  for more details.

- envname:

  The name of the conda environment to remove. If `NULL`, uses the
  default scop environment name.
