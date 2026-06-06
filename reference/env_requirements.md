# Python environment requirements

The function returns a list of requirements including the required
Python version, package versions, and package name aliases for
platform-specific packages. All packages will be installed using uv as
the primary tool.

## Usage

``` r
env_requirements(
  version = "3.10-1",
  include_optional = FALSE,
  modules = NULL,
  verbose = TRUE
)
```

## Arguments

- version:

  The Python version of the environment. Default is `"3.10-1"`.

- include_optional:

  Whether to include optional Python dependencies.

- modules:

  Optional requirement modules to include. Supported values are
  `"scanpy"`, `"scvi"`, `"scanorama"`, `"bbknn"`, `"celltypist"`,
  `"cellphonedb"`, `"magic"`, `"scrublet"`, `"doubletdetection"`,
  `"sccoda"`, `"doublet"`, `"palantir"`, `"scvelo"`, `"cellrank"`,
  `"wot"`, `"phate"`, `"pacmap"`, `"trimap"`, `"multimap"`, `"scomm"`,
  `"scenic"`, `"seacells"`, and `"tage"`. If `NULL`, the default
  environment is returned. The default excludes `"sccoda"`, `"scomm"`,
  and `"scenic"` because these workflows require dependency stacks that
  should be prepared explicitly. The `"scenic"` module is standalone and
  always uses Python `"3.10-1"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list containing:

- python:

  Python version string

- packages:

  Named vector of package version specifications

- package_aliases:

  Named list mapping logical package names to actual installed names

## Examples

``` r
env_requirements("3.10-1")
```
