# Standard workflow for scop

This function performs a standard single-cell analysis workflow.

## Usage

``` r
standard_scop(
  srt,
  prefix = "Standard",
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_method = "vst",
  nHVF = 2000,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- prefix:

  A prefix to add to the names of intermediate objects created by the
  function. Default is `"Standard"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- do_normalization:

  Whether to perform normalization. If `NULL`, normalization will be
  performed if the specified assay does not have scaled data.

- normalization_method:

  The method to use for normalization. Options are `"LogNormalize"`,
  `"SCT"`, or `"TFIDF"`. Default is `"LogNormalize"`.

- do_HVF_finding:

  Whether to perform high variable feature finding. If `TRUE`, the
  function will force to find the highly variable features (HVF) using
  the specified HVF method.

- HVF_method:

  The method to use for finding highly variable features. Options are
  `"vst"`, `"mvp"`, or `"disp"`. Default is `"vst"`.

- nHVF:

  The number of highly variable features to select. If NULL, all highly
  variable features will be used. Default is `2000`.

- HVF:

  A vector of feature names to use as highly variable features. If NULL,
  the function will use the highly variable features identified by the
  HVF method.

- do_scaling:

  Whether to perform scaling. If `TRUE`, the function will force to
  scale the data using the
  [Seurat::ScaleData](https://satijalab.org/seurat/reference/ScaleData.html)
  function.

- vars_to_regress:

  A vector of feature names to use as regressors in the scaling step. If
  NULL, no regressors will be used.

- regression_model:

  The regression model to use for scaling. Options are `"linear"`,
  `"poisson"`, or `"negativebinomial"`. Default is `"linear"`.

- linear_reduction:

  The linear dimensionality reduction method to use. Options are
  `"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`, or `"glmpca"`. Default is
  `"pca"`.

- linear_reduction_dims:

  The number of dimensions to keep after linear dimensionality
  reduction. Default is `50`.

- linear_reduction_dims_use:

  The dimensions to use for downstream analysis. If `NULL`, all
  dimensions will be used.

- linear_reduction_params:

  A list of parameters to pass to the linear dimensionality reduction
  method.

- force_linear_reduction:

  Whether to force linear dimensionality reduction even if the specified
  reduction is already present in the Seurat object.

- nonlinear_reduction:

  The nonlinear dimensionality reduction method to use. Options are
  `"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`, `"pacmap"`,
  `"trimap"`, `"largevis"`, or `"fr"`. Default is `"umap"`.

- nonlinear_reduction_dims:

  The number of dimensions to keep after nonlinear dimensionality
  reduction. If a vector is provided, different numbers of dimensions
  can be specified for each method. Default is `c(2, 3)`.

- nonlinear_reduction_params:

  A list of parameters to pass to the nonlinear dimensionality reduction
  method.

- force_nonlinear_reduction:

  Whether to force nonlinear dimensionality reduction even if the
  specified reduction is already present in the Seurat object. Default
  is `TRUE`.

- neighbor_metric:

  The distance metric to use for finding neighbors. Options are
  `"euclidean"`, `"cosine"`, `"manhattan"`, or `"hamming"`. Default is
  `"euclidean"`.

- neighbor_k:

  The number of nearest neighbors to use for finding neighbors. Default
  is `20`.

- cluster_algorithm:

  The clustering algorithm to use. Options are `"louvain"`, `"slm"`, or
  `"leiden"`. Default is `"louvain"`.

- cluster_resolution:

  The resolution parameter to use for clustering. Larger values result
  in fewer clusters. Default is `0.6`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A `Seurat` object.

## Examples

``` r
library(Matrix)
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType"
)


# Use a combination of different linear
# or non-linear dimension reduction methods
linear_reductions <- c(
  "pca", "nmf", "mds"
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  linear_reduction = linear_reductions,
  nonlinear_reduction = "umap"
)
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
plist1 <- lapply(
  linear_reductions, function(lr) {
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = paste0(
        "Standard", lr, "UMAP2D"
      ),
      xlab = "", ylab = "", title = lr,
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plist1)


nonlinear_reductions <- c(
  "umap", "tsne", "dm", "fr"
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  linear_reduction = "pca",
  nonlinear_reduction = nonlinear_reductions
)
#> Error in tryCatchOne(expr, names, parentenv, handlers[[1L]]): Error when performing dm nonlinear dimension reduction. Skip it
plist2 <- lapply(
  nonlinear_reductions, function(nr) {
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = paste0(
        "Standardpca", nr, "2D"
      ),
      xlab = "", ylab = "", title = nr,
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plist2)
```
