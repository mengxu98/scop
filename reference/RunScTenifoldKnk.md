# Run scTenifoldKnk in-silico knockout analysis

Run scTenifoldKnk in-silico knockout analysis

## Usage

``` r
RunscTenifoldKnk(
  srt,
  gKO,
  assay = NULL,
  layer = "counts",
  features = NULL,
  qc = TRUE,
  qc_mt_threshold = 0.1,
  qc_min_library_size = 1000,
  qc_min_cells = 25,
  nc_lambda = 0,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q = 0.9,
  td_K = 3,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  td_nDecimal = 3,
  ma_nDim = 2,
  cores = 1,
  backend = c("r", "cpp"),
  store_networks = TRUE,
  store_manifold = TRUE,
  tool_name = "scTenifoldKnk",
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- gKO:

  Gene symbol or symbols to knock out. All genes must be present after
  optional feature and QC filtering.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the count matrix.

- features:

  Optional genes to retain before running network construction. If
  supplied, `gKO` is always retained when present in the input assay.

- qc:

  Whether to apply scTenifoldKnk-style quality control.

- qc_mt_threshold:

  Maximum mitochondrial read fraction per cell.

- qc_min_library_size:

  Minimum library size per cell.

- qc_min_cells:

  Minimum number of expressing cells required per gene.

- nc_lambda, nc_nNet, nc_nCells, nc_nComp, nc_scaleScores, nc_symmetric,
  nc_q:

  Network construction parameters forwarded to
  [`scTenifoldNet::makeNetworks()`](https://rdrr.io/pkg/scTenifoldNet/man/makeNetworks.html).

- td_K, td_maxIter, td_maxError, td_nDecimal:

  Tensor decomposition parameters forwarded to
  [`scTenifoldNet::tensorDecomposition()`](https://rdrr.io/pkg/scTenifoldNet/man/tensorDecomposition.html).

- ma_nDim:

  Manifold-alignment dimension forwarded to
  [`scTenifoldNet::manifoldAlignment()`](https://rdrr.io/pkg/scTenifoldNet/man/manifoldAlignment.html).

- cores:

  Number of cores used by native network-construction workers and
  forwarded to downstream linear algebra where applicable.

- backend:

  `r` calls
  [`scTenifoldKnk::scTenifoldKnk()`](https://rdrr.io/pkg/scTenifoldKnk/man/scTenifoldKnk.html)
  directly and is the default high-consistency path. `cpp` follows the
  upstream `scTenifoldNet`/`scTenifoldKnk` network construction, tensor
  decomposition, manifold alignment, and differential-regulation steps
  while keeping input handling and result storage inside `scop`.

- store_networks:

  Whether to keep WT/KO tensor networks in `srt@tools`.

- store_manifold:

  Whether to keep manifold-alignment coordinates in `srt@tools`.

- tool_name:

  Name of the `srt@tools` entry.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with scTenifoldKnk results stored in
`srt@tools[[tool_name]]`.

## Examples

``` r
data(pancreas_sub)
gene_use <- "Pdx1"
counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)
detected <- names(
  sort(Matrix::rowSums(counts > 0),
    decreasing = TRUE
  )
)
features_use <- unique(c(gene_use, head(detected, 300)))

pancreas_sub <- RunscTenifoldKnk(
  pancreas_sub,
  gKO = gene_use,
  features = features_use,
  qc = FALSE,
  nc_nNet = 3,
  nc_nCells = 200,
  td_maxIter = 200,
  store_networks = FALSE,
  store_manifold = TRUE
)
#> ℹ [2026-06-29 04:35:18] Run scTenifoldKnk knockout for "Pdx1" using "r" backend
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |======================================================================| 100%
#> ✔ [2026-06-29 04:35:30] scTenifoldKnk results stored in `srt@tools[[scTenifoldKnk]]`

dr <- pancreas_sub@tools$scTenifoldKnk$diffRegulation
head(dr)
#>        gene     distance        Z         FC      p.value        p.adj
#> 1      Pdx1 6.765145e-04 3.880684 7505.34751 0.000000e+00 0.000000e+00
#> 33     Gnas 4.408814e-05 1.882629   31.87572 1.643584e-08 2.473594e-06
#> 148    Ssr2 4.051410e-05 1.832606   26.91714 2.123672e-07 2.130751e-05
#> 156  Sec61b 3.769553e-05 1.790425   23.30217 1.384438e-06 1.041789e-04
#> 236    Dad1 3.334593e-05 1.719717   18.23485 1.952730e-05 1.175543e-03
#> 279 Rpl36al 3.128641e-05 1.683448   16.05196 6.162780e-05 3.091661e-03

scTenifoldKnkPlot(pancreas_sub, plot_type = "effect")
```
