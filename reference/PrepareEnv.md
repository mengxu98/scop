# Prepare the python environment

Prepare the python environment by installing the required dependencies
and setting up the environment.

## Usage

``` r
PrepareEnv(
  envname = NULL,
  conda = "auto",
  miniconda_repo = "https://repo.anaconda.com/miniconda",
  version = "3.10-1",
  force = FALSE,
  ...
)
```

## Arguments

- envname:

  The name of the conda environment. If `NULL`, the environment name
  will be set to `"scop_env"`. Default is `NULL`.

- conda:

  The path to a conda executable. Use `"auto"` to allow automatically
  finding an appropriate conda binary.

- miniconda_repo:

  Repository URL for miniconda. Default is
  <https://repo.anaconda.com/miniconda>.

- version:

  The Python version. Default is `"3.10-1"`.

- force:

  Whether to force recreation of the environment. If `TRUE`, the
  existing environment will be removed and recreated. Default is
  `FALSE`.

- ...:

  Additional arguments passed to package installation functions.
