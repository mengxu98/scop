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

  The name of the conda-compatible Python environment. If `NULL`, the
  environment name will be set to `"scop_env"`. Default is `NULL`.

- conda:

  The path or command name of a conda-compatible executable (`conda`,
  `mamba`, or `micromamba`). Use `"auto"` to allow automatically finding
  an appropriate environment manager. If `"micromamba"` is requested and
  micromamba is not available on `PATH`, a package-managed micromamba is
  downloaded automatically.

- force:

  Whether to force package reinstallation. Default is `FALSE`.

- pip:

  Whether to use `pip`/`uv` (`TRUE`) or the configured conda-compatible
  environment manager (`FALSE`) for installation. Default is `TRUE`.
  When `TRUE`, uv is used as the primary installer with pip as fallback.

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
