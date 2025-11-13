# Check if the python package exists in the environment

Check if the python package exists in the environment

## Usage

``` r
exist_python_pkgs(packages, envname = NULL, conda = "auto")
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
