# Check and install python packages

Check and install python packages

## Usage

``` r
check_python(
  packages,
  envname = NULL,
  conda = "auto",
  force = FALSE,
  pip = TRUE,
  pip_options = character(),
  verbose = TRUE,
  ...
)
```

## Arguments

- packages:

  A character vector of package names to check and install. Use
  `"<package>==<version>"` to request a specific version.

- envname:

  The name of the conda environment. If `NULL`, the environment name
  will be set to `"scop_env"`. Default is `NULL`.

- conda:

  The path to a conda executable. Use `"auto"` to allow automatically
  finding an appropriate conda binary.

- force:

  Whether to force package reinstallation. Default is `FALSE`.

- pip:

  Whether to use `pip`/`uv` (`TRUE`) or `conda` (`FALSE`) for
  installation. Default is `TRUE`. When `TRUE`, uv is used as the
  primary installer with pip as fallback.

- pip_options:

  Additional command line arguments to be passed to `uv`/`pip` when
  `pip = TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Other arguments to be passed to other functions.

## Examples

``` r
if (FALSE) { # \dontrun{
PrepareEnv()

# Then check/install packages
check_python(
  packages = c("numpy", "pandas")
)

check_python(
  packages = "numpy==1.26.4",
  envname = "scop_env",
  pip_options = "-i https://pypi.tuna.tsinghua.edu.cn/simple"
)
} # }
```
