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
  backend = c("cpp", "scTenifoldKnk"),
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

  `cpp` is a native equivalent covariance-based network construction,
  direct sparse network assembly, controlled per-gene eigensolver
  parallelism, and helpers for tensor decomposition, manifold matrix
  construction, directionality, and distance calculation.
  `scTenifoldKnk` calls
  [`scTenifoldKnk::scTenifoldKnk()`](https://rdrr.io/pkg/scTenifoldKnk/man/scTenifoldKnk.html)
  directly for comparison.

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
counts <- GetAssayData5(pancreas_sub, assay = "RNA", layer = "counts")
detected <- names(sort(Matrix::rowSums(counts > 0), decreasing = TRUE))
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
#> ℹ [2026-05-11 16:20:06] Run scTenifoldKnk knockout for "Pdx1" using "cpp" backend
#> ℹ [2026-05-11 16:20:09] Construct scTenifoldNet network ensemble
#> ℹ [2026-05-11 16:20:09] Denoise network ensemble with tensor decomposition
#> ℹ [2026-05-11 16:20:09] Denoise network ensemble with tensor decomposition  ■■ …
#> ℹ [2026-05-11 16:20:09] Align WT and KO network manifolds
#> ✔ [2026-05-11 16:20:10] scTenifoldKnk results stored in `srt@tools[[scTenifoldKnk]]`

dr <- pancreas_sub@tools$scTenifoldKnk$diffRegulation
head(dr)
#>     gene     distance        Z          FC      p.value        p.adj
#> 1   Pdx1 9.304567e-04 4.569367 12059.91897 0.000000e+00 0.000000e+00
#> 2   Gnas 4.887597e-05 2.058208    33.27690 7.992633e-09 1.202891e-06
#> 3   Ssr2 4.386627e-05 1.985523    26.80486 2.250708e-07 2.258210e-05
#> 4   Ssr4 4.132486e-05 1.945911    23.78893 1.074990e-06 8.089300e-05
#> 5 Sec61b 4.074745e-05 1.936623    23.12879 1.515051e-06 9.120607e-05
#> 6   Dad1 3.408756e-05 1.820516    16.18616 5.741213e-05 2.880175e-03

scTenifoldKnkPlot(pancreas_sub, plot_type = "effect")
```
