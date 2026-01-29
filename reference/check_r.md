# Check and install R packages

Check and install R packages

## Usage

``` r
check_r(
  packages,
  lib = .libPaths()[1],
  dependencies = TRUE,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- packages:

  Package to be installed. Package source can be *CRAN*, *Bioconductor*
  or *Github*. By default, the package name is extracted according to
  the `packages` parameter.

- lib:

  The location of the library directories where to install the packages.

- dependencies:

  Whether to install dependencies of the packages. Default is `TRUE`.

- force:

  Whether to force the installation of packages. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

Package installation status.

## Examples

``` r
check_r(c("ggplot2", "dplyr"))
#> ✔ [2026-01-29 13:44:35] ggplot2 and dplyr installed successfully
```
