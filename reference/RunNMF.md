# Run NMF (non-negative matrix factorization)

Run NMF (non-negative matrix factorization)

## Usage

``` r
RunNMF(object, ...)

# S3 method for class 'Seurat'
RunNMF(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nbes = 50,
  nmf.method = "RcppML",
  tol = 1e-05,
  maxit = 100,
  rev.nmf = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "nmf",
  reduction.key = "BE_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay'
RunNMF(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nbes = 50,
  nmf.method = "RcppML",
  tol = 1e-05,
  maxit = 100,
  rev.nmf = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "BE_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay5'
RunNMF(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nbes = 50,
  nmf.method = "RcppML",
  tol = 1e-05,
  maxit = 100,
  rev.nmf = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "BE_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunNMF(
  object,
  assay = NULL,
  layer = "data",
  nbes = 50,
  nmf.method = "RcppML",
  tol = 1e-05,
  maxit = 100,
  rev.nmf = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "BE_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, an Assay object, or a
  matrix-like object.

- ...:

  Additional arguments passed to
  [RcppML::nmf](https://rdrr.io/pkg/RcppML/man/nmf.html) or
  [NMF::nmf](https://rdrr.io/pkg/NMF/man/nmf.html).

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `data`.

- features:

  A character vector of features to use. Default is `NULL`.

- nbes:

  The number of basis vectors (components) to be computed. Default is
  `50`.

- nmf.method:

  The NMF algorithm to be used. Currently supported values are
  `"RcppML"` and `"NMF"`. Default is `"RcppML"`.

- tol:

  The tolerance for convergence (only applicable when nmf.method is
  `"RcppML"`). Default is `1e-5`.

- maxit:

  The maximum number of iterations for convergence (only applicable when
  nmf.method is `"RcppML"`). Default is `100`.

- rev.nmf:

  Whether to perform reverse NMF (i.e., transpose the input matrix)
  before running the analysis. Default is `FALSE`.

- ndims.print:

  The dimensions (number of basis vectors) to print in the output.
  Default is `1:5`.

- nfeatures.print:

  The number of features to print in the output. Default is `30`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"nmf"`.

- reduction.key:

  The prefix for the column names of the basis vectors. Default is
  `"BE_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
library(Matrix)
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunNMF(pancreas_sub)
#> 
#> iter |      tol 
#> ---------------
#>    1 | 6.94e-01
#>    2 | 9.76e-02
#>    3 | 3.18e-02
#>    4 | 1.49e-02
#>    5 | 8.15e-03
#>    6 | 5.06e-03
#>    7 | 3.40e-03
#>    8 | 2.43e-03
#>    9 | 1.80e-03
#>   10 | 1.37e-03
#>   11 | 1.07e-03
#>   12 | 8.59e-04
#>   13 | 7.01e-04
#>   14 | 5.85e-04
#>   15 | 5.07e-04
#>   16 | 4.47e-04
#>   17 | 4.03e-04
#>   18 | 3.67e-04
#>   19 | 3.32e-04
#>   20 | 2.98e-04
#>   21 | 2.72e-04
#>   22 | 2.48e-04
#>   23 | 2.28e-04
#>   24 | 2.11e-04
#>   25 | 1.99e-04
#>   26 | 1.91e-04
#>   27 | 1.85e-04
#>   28 | 1.76e-04
#>   29 | 1.72e-04
#>   30 | 1.68e-04
#>   31 | 1.62e-04
#>   32 | 1.54e-04
#>   33 | 1.45e-04
#>   34 | 1.35e-04
#>   35 | 1.24e-04
#>   36 | 1.15e-04
#>   37 | 1.06e-04
#>   38 | 9.81e-05
#>   39 | 9.10e-05
#>   40 | 8.53e-05
#>   41 | 8.05e-05
#>   42 | 7.62e-05
#>   43 | 7.25e-05
#>   44 | 6.94e-05
#>   45 | 6.73e-05
#>   46 | 6.61e-05
#>   47 | 6.49e-05
#>   48 | 6.36e-05
#>   49 | 6.22e-05
#>   50 | 6.04e-05
#>   51 | 5.84e-05
#>   52 | 5.62e-05
#>   53 | 5.39e-05
#>   54 | 5.05e-05
#>   55 | 4.69e-05
#>   56 | 4.39e-05
#>   57 | 4.11e-05
#>   58 | 3.85e-05
#>   59 | 3.60e-05
#>   60 | 3.38e-05
#>   61 | 3.18e-05
#>   62 | 3.01e-05
#>   63 | 2.88e-05
#>   64 | 2.77e-05
#>   65 | 2.69e-05
#>   66 | 2.59e-05
#>   67 | 2.50e-05
#>   68 | 2.40e-05
#>   69 | 2.30e-05
#>   70 | 2.20e-05
#>   71 | 2.12e-05
#>   72 | 2.05e-05
#>   73 | 2.00e-05
#>   74 | 1.96e-05
#>   75 | 1.94e-05
#>   76 | 1.90e-05
#>   77 | 1.87e-05
#>   78 | 1.82e-05
#>   79 | 1.72e-05
#>   80 | 1.62e-05
#>   81 | 1.54e-05
#>   82 | 1.45e-05
#>   83 | 1.37e-05
#>   84 | 1.31e-05
#>   85 | 1.25e-05
#>   86 | 1.20e-05
#>   87 | 1.16e-05
#>   88 | 1.11e-05
#>   89 | 1.07e-05
#>   90 | 1.02e-05
#>   91 | 9.91e-06
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "nmf"
)
```
