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

  A character vector, indicating package names which should be installed
  or removed. Use `"<package>==<version>"` to request the installation
  of a specific version of a package.

- envname:

  The name of a conda environment.

- conda:

  The path to a conda executable. Use `"auto"` to allow scop to
  automatically find an appropriate conda binary.

- force:

  Whether to force package installation. Default is `FALSE`.

- pip:

  Whether to use pip for package installation. Default is `TRUE`,
  packages are installed from the active conda channels.

- pip_options:

  An optional character vector of additional command line arguments to
  be passed to `pip`. Only relevant when `pip = TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Other arguments passed to
  [reticulate::conda_install](https://rstudio.github.io/reticulate/reference/conda-tools.html)

## Examples

``` r
if (FALSE) { # \dontrun{
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
