# Find the path to Python associated with a conda environment

Find the path to Python associated with a conda environment

## Usage

``` r
conda_python(envname = NULL, conda = "auto", all = FALSE)
```

## Arguments

- envname:

  The name of, or path to, a conda environment.

- conda:

  The path to a `conda` executable. Use `"auto"` to allow `reticulate`
  to automatically find an appropriate `conda` binary. See **Finding
  Conda** and
  [`conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  for more details.

- all:

  Boolean; report all instances of Python found?
