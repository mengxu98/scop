# Enhanced conda installation

Enhanced conda installation

## Usage

``` r
conda_install(
  envname = NULL,
  packages,
  forge = TRUE,
  channel = character(),
  pip = FALSE,
  pip_options = character(),
  pip_ignore_installed = FALSE,
  conda = "auto",
  python_version = NULL,
  ...
)
```

## Arguments

- envname:

  The name of, or path to, a conda environment.

- packages:

  A character vector, indicating package names which should be installed
  or removed. Use `<package>==<version>` to request the installation of
  a specific version of a package. A `NULL` value for
  [`conda_remove()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  will be interpretted to `"--all"`, removing the entire environment.

- forge:

  Boolean; include the [conda-forge](https://conda-forge.org/)
  repository?

- channel:

  An optional character vector of conda channels to include. When
  specified, the `forge` argument is ignored. If you need to specify
  multiple channels, including the conda forge, you can use
  `c("conda-forge", <other channels>)`.

- pip:

  Boolean; use `pip` for package installation? By default, packages are
  installed from the active conda channels.

- pip_options:

  An optional character vector of additional command line arguments to
  be passed to `pip`. Only relevant when `pip = TRUE`.

- pip_ignore_installed:

  Ignore already-installed versions when using pip? (defaults to
  `FALSE`). Set this to `TRUE` so that specific package versions can be
  installed even if they are downgrades. The `FALSE` option is useful
  for situations where you don't want a pip install to attempt an
  overwrite of a conda binary package (e.g. SciPy on Windows which is
  very difficult to install via pip due to compilation requirements).

- conda:

  The path to a `conda` executable. Use `"auto"` to allow `reticulate`
  to automatically find an appropriate `conda` binary. See **Finding
  Conda** and
  [`conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
  for more details.

- python_version:

  The version of Python to be installed. Set this if you'd like to
  change the version of Python associated with a particular conda
  environment.

- ...:

  Optional arguments, reserved for future expansion.
