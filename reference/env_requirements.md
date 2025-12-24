# Python environment requirements

The function returns a list of requirements including the required
Python version, package versions, and package name aliases for
platform-specific packages. All packages will be installed using uv as
the primary tool.

## Usage

``` r
env_requirements(version = "3.10-1")
```

## Arguments

- version:

  The Python version of the environment. Default is `"3.10-1"`.

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
#> $python
#> [1] "3.10-1"
#> 
#> $packages
#>               leidenalg                     tbb           python-igraph 
#>     "leidenalg==0.10.2"         "tbb==2022.2.0" "python-igraph==0.11.9" 
#>              matplotlib                   numba                llvmlite 
#>    "matplotlib==3.10.8"         "numba==0.59.1"      "llvmlite==0.42.0" 
#>                   numpy               packaging                palantir 
#>         "numpy==1.26.4"       "packaging>=24.0"       "palantir==1.4.1" 
#>                  pandas                  scanpy            scikit-learn 
#>         "pandas==2.0.3"        "scanpy==1.11.3"   "scikit-learn==1.7.0" 
#>                   scipy                  scvelo                     wot 
#>         "scipy==1.15.3"         "scvelo==0.3.3"      "wot==1.0.8.post2" 
#>                  trimap                  pacmap                   phate 
#>         "trimap==1.1.4"         "pacmap==0.8.0"         "phate==1.0.11" 
#>                   bbknn               scanorama              scvi-tools 
#>          "bbknn==1.6.0"      "scanorama==1.7.4"     "scvi-tools==1.2.1" 
#>                cellrank              celltypist 
#>       "cellrank==2.0.7"            "celltypist" 
#> 
#> $install_methods
#>     leidenalg           tbb python-igraph    scvi-tools    matplotlib 
#>       "conda"       "conda"       "conda"       "conda"         "pip" 
#>         numba      llvmlite         numpy     packaging      palantir 
#>         "pip"         "pip"         "pip"         "pip"         "pip" 
#>        pandas        scanpy  scikit-learn         scipy        scvelo 
#>         "pip"         "pip"         "pip"         "pip"         "pip" 
#>           wot        trimap        pacmap         phate         bbknn 
#>         "pip"         "pip"         "pip"         "pip"         "pip" 
#>     scanorama      cellrank    celltypist 
#>         "pip"         "pip"         "pip" 
#> 
#> $package_aliases
#> $package_aliases$`python-igraph`
#> [1] "igraph"
#> 
#> 
```
