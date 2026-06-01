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
  and `"scenic"`. If `NULL`, the default environment is returned. The
  default excludes `"sccoda"`, `"scomm"`, and `"scenic"` because these
  workflows require dependency stacks that should be prepared
  explicitly. The `"scenic"` module is standalone and always uses Python
  `"3.10-1"`.

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
#> $python
#> [1] "3.10-1"
#> 
#> $packages
#>                                     setuptools 
#>                                "setuptools<81" 
#>                                      leidenalg 
#>                            "leidenalg==0.10.2" 
#>                                            tbb 
#>                                "tbb==2022.2.0" 
#>                                  python-igraph 
#>                        "python-igraph==0.11.9" 
#>                                     matplotlib 
#>                           "matplotlib==3.10.8" 
#>                                          numba 
#>                                "numba==0.59.1" 
#>                                       llvmlite 
#>                             "llvmlite==0.42.0" 
#>                                          numpy 
#>                                "numpy==1.26.4" 
#>                                      packaging 
#>                              "packaging>=24.0" 
#>                                         pandas 
#>                                "pandas==2.0.3" 
#>                                   scikit-learn 
#>                          "scikit-learn==1.7.0" 
#>                                          scipy 
#>                                "scipy==1.15.3" 
#>                                         scanpy 
#>                               "scanpy==1.11.3" 
#>                                         loompy 
#>                                       "loompy" 
#>                                     scvi-tools 
#>                            "scvi-tools==1.2.1" 
#>                                            jax 
#>                             "jax[cpu]==0.4.38" 
#>                                         scglue 
#>                                "scglue==0.4.0" 
#>                                       bedtools 
#>                                     "bedtools" 
#>                                      scanorama 
#>                             "scanorama==1.7.4" 
#>                                          bbknn 
#>                                 "bbknn==1.6.0" 
#>                                     celltypist 
#>                            "celltypist==1.7.1" 
#>                                    cellphonedb 
#>                           "cellphonedb==5.0.1" 
#>                                   magic-impute 
#>                          "magic-impute==3.0.0" 
#>                                       scrublet 
#>                              "scrublet==0.2.3" 
#>                               doubletdetection 
#>                "doubletdetection==4.3.0.post1" 
#>                                        louvain 
#>                               "louvain==0.8.2" 
#>                                       palantir 
#>                              "palantir==1.4.1" 
#>                                         scvelo 
#>                                "scvelo==0.3.3" 
#>                                       cellrank 
#>                              "cellrank==2.0.7" 
#>                                            wot 
#>                             "wot==1.0.8.post2" 
#>                                          phate 
#>                                "phate==1.0.11" 
#>                                         pacmap 
#>                                "pacmap==0.8.0" 
#>                                         trimap 
#>                                "trimap==1.1.4" 
#>                                       multimap 
#> "git+https://github.com/Teichlab/MultiMAP.git" 
#> 
#> $install_methods
#>       setuptools        leidenalg              tbb    python-igraph 
#>            "pip"          "conda"          "conda"          "conda" 
#>       matplotlib            numba         llvmlite            numpy 
#>            "pip"            "pip"            "pip"            "pip" 
#>        packaging           pandas     scikit-learn            scipy 
#>            "pip"            "pip"            "pip"            "pip" 
#>           scanpy           loompy       scvi-tools              jax 
#>            "pip"            "pip"          "conda"            "pip" 
#>           scglue         bedtools        scanorama            bbknn 
#>            "pip"          "conda"            "pip"            "pip" 
#>       celltypist      cellphonedb     magic-impute         scrublet 
#>            "pip"            "pip"            "pip"            "pip" 
#> doubletdetection          louvain         palantir           scvelo 
#>            "pip"            "pip"            "pip"            "pip" 
#>         cellrank              wot            phate           pacmap 
#>            "pip"            "pip"            "pip"            "pip" 
#>           trimap         multimap 
#>            "pip"            "pip" 
#> 
#> $package_aliases
#> $package_aliases$`python-igraph`
#> [1] "igraph"
#> 
#> $package_aliases$multimap
#> [1] "MultiMAP"
#> 
#> $package_aliases$MultiMAP
#> [1] "multimap"
#> 
#> 
```
