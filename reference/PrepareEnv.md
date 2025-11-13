# Prepare the virtual environment

Prepare the virtual environment by installing the required dependencies
and setting up the environment. This function prepares the virtual
environment by checking if conda is installed, creating a new conda
environment if needed, installing the required packages, and setting up
the Python environment for use with scop. In order to create the
environment, this function requires the path to the conda binary. If
`conda` is set to `"auto"`, it will attempt to automatically find the
conda binary. If a conda environment with the specified name already
exists and `force` is set to `FALSE`, the function will use the existing
environment. If `force` set to `TRUE`, the existing environment will be
recreated. Note that recreating the environment will remove any existing
data in the environment. The function also checks if the package
versions in the environment meet the requirements specified by the
`version` parameter.

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

  The name of a conda environment.

- conda:

  The path to a conda executable. Use `"auto"` to allow scop to
  automatically find an appropriate conda binary.

- miniconda_repo:

  Repositories for miniconda. Default is
  <https://repo.anaconda.com/miniconda>.

- version:

  A character vector specifying the version of the environment. Default
  is `"3.10-1"`.

- force:

  Whether to force a new environment to be created. If `TRUE`, the
  existing environment will be recreated. Default is `FALSE`.

- ...:

  Other arguments passed to
  [reticulate::conda_install](https://rstudio.github.io/reticulate/reference/conda-tools.html)
