# Run scTenifoldNet network comparison

Run scTenifoldNet network comparison

## Usage

``` r
RunscTenifoldNet(
  object,
  y = NULL,
  group.by = NULL,
  condition1 = NULL,
  condition2 = NULL,
  assay = NULL,
  layer = "counts",
  features = NULL,
  qc = TRUE,
  qc_min_library_size = 1000,
  qc_remove_outlier_cells = TRUE,
  qc_min_pct = 0.05,
  qc_max_mt_ratio = 0.1,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_symmetric = FALSE,
  nc_scaleScores = TRUE,
  nc_q = 0.05,
  td_K = 3,
  td_nDecimal = 1,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  ma_nDim = 30,
  cores = 1,
  store_networks = TRUE,
  store_manifold = TRUE,
  tool_name = "scTenifoldNet",
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object or a raw count matrix with genes in rows and cells
  in columns.

- y:

  A second raw count matrix. Required when `object` is a matrix and
  ignored when `object` is a `Seurat` object.

- group.by:

  Metadata column used to split a `Seurat` object into the two
  conditions being compared.

- condition1, condition2:

  Condition labels from `group.by`. If omitted, the first two group
  levels are used.

- assay, layer:

  Assay and layer used as the count matrix when `object` is a `Seurat`
  object.

- features:

  Optional genes to retain before running the comparison.

- qc:

  Whether to apply scTenifoldNet-style quality control.

- qc_min_library_size, qc_remove_outlier_cells, qc_min_pct,
  qc_max_mt_ratio:

  Quality-control parameters forwarded to
  [`scTenifoldNet::scQC()`](https://rdrr.io/pkg/scTenifoldNet/man/scQC.html).

- nc_nNet, nc_nCells, nc_nComp, nc_symmetric, nc_scaleScores, nc_q:

  Network construction parameters forwarded to
  [`scTenifoldNet::makeNetworks()`](https://rdrr.io/pkg/scTenifoldNet/man/makeNetworks.html).

- td_K, td_nDecimal, td_maxIter, td_maxError:

  Tensor decomposition parameters forwarded to
  [`scTenifoldNet::tensorDecomposition()`](https://rdrr.io/pkg/scTenifoldNet/man/tensorDecomposition.html).

- ma_nDim:

  Manifold-alignment dimension forwarded to
  [`scTenifoldNet::manifoldAlignment()`](https://rdrr.io/pkg/scTenifoldNet/man/manifoldAlignment.html).

- cores:

  Number of cores forwarded to
  [`scTenifoldNet::scTenifoldNet()`](https://rdrr.io/pkg/scTenifoldNet/man/scTenifoldNet.html).

- store_networks:

  Whether to keep tensor networks in the stored result when `object` is
  a `Seurat` object.

- store_manifold:

  Whether to keep manifold-alignment coordinates in the stored result
  when `object` is a `Seurat` object.

- tool_name:

  Name of the `object@tools` entry when `object` is a `Seurat` object.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A scTenifoldNet result list for matrix input, or a `Seurat` object with
results stored in `object@tools[[tool_name]]`.

## Examples

``` r
data(pancreas_sub)
counts <- GetAssayData5(pancreas_sub, assay = "RNA", layer = "counts")
detected <- names(sort(Matrix::rowSums(counts > 0), decreasing = TRUE))
features_use <- head(detected, 300)

pancreas_sub <- RunscTenifoldNet(
  pancreas_sub,
  group.by = "CellType",
  condition1 = "Ductal",
  condition2 = "Endocrine",
  features = features_use,
  qc = FALSE,
  nc_nNet = 3,
  nc_nCells = 200,
  td_maxIter = 200,
  ma_nDim = 2,
  store_networks = FALSE,
  store_manifold = TRUE
)
#> ℹ [2026-06-24 04:41:47] Run scTenifoldNet comparison using upstream implementation
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |======================================================================| 100%
#> ✔ [2026-06-24 04:41:51] scTenifoldNet results stored in `object@tools[[scTenifoldNet]]`

dr <- pancreas_sub@tools$scTenifoldNet$diffRegulation
head(dr)
#>         gene    distance        Z        FC      p.value       p.adj
#> 25  Hsp90ab1 0.004553045 1.660546 20.824368 5.033841e-06 0.001510152
#> 26      Rps2 0.004154439 1.590437 17.337751 3.129086e-05 0.004693629
#> 125    Hspa8 0.003136144 1.383049  9.880067 1.670788e-03 0.167078817
#> 47      Ftl1 0.002977463 1.346001  8.905553 2.843048e-03 0.193654167
#> 177 Hsp90aa1 0.002938524 1.336669  8.674141 3.227569e-03 0.193654167
#> 155   Sec61b 0.002786133 1.299164  7.797794 5.231006e-03 0.240336435

scTenifoldNetPlot(pancreas_sub, plot_type = "effect")
```
