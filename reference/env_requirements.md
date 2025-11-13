# Python environment requirements

The function returns a list of requirements including the required
Python version and a list of packages with their corresponding versions.

## Usage

``` r
env_requirements(version = "3.10-1")
```

## Arguments

- version:

  A character vector specifying the version of the environment. Default
  is "3.10-1".

## Value

A list of requirements for the specified version.

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
#>    "matplotlib==3.10.3"         "numba==0.59.1"      "llvmlite==0.42.0" 
#>                   numpy                palantir                  pandas 
#>         "numpy==1.26.4"       "palantir==1.4.1"         "pandas==2.0.3" 
#>                  scanpy            scikit-learn                   scipy 
#>        "scanpy==1.11.3"   "scikit-learn==1.7.0"         "scipy==1.15.3" 
#>                  scvelo                     wot                  trimap 
#>         "scvelo==0.3.3"      "wot==1.0.8.post2"         "trimap==1.1.4" 
#>                  pacmap                   phate                   bbknn 
#>         "pacmap==0.8.0"         "phate==1.0.11"          "bbknn==1.6.0" 
#>               scanorama              scvi-tools                cellrank 
#>      "scanorama==1.7.4"     "scvi-tools==1.2.1"       "cellrank==2.0.7" 
#> 
#> $install_methods
#>     leidenalg           tbb python-igraph    matplotlib         numba 
#>       "conda"       "conda"       "conda"         "pip"         "pip" 
#>      llvmlite         numpy      palantir        pandas        scanpy 
#>         "pip"         "pip"         "pip"         "pip"         "pip" 
#>  scikit-learn         scipy        scvelo           wot        trimap 
#>         "pip"         "pip"         "pip"         "pip"         "pip" 
#>        pacmap         phate         bbknn     scanorama    scvi-tools 
#>         "pip"         "pip"         "pip"         "pip"         "pip" 
#>      cellrank 
#>         "pip" 
#> 
```
