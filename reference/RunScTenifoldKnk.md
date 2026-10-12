# Run scTenifoldKnk in-silico knockout analysis

Run scTenifoldKnk in-silico knockout analysis

## Usage

``` r
RunscTenifoldKnk(
  object,
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
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- gKO:

  Gene symbol or symbols to knock out. All genes must be present after
  optional feature and QC filtering.

- assay:

  Assay to use. `NULL` uses the default assay.

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

  `r` runs the upstream `scTenifoldNet` pipeline with the
  `scTenifoldKnk` differential-regulation helpers and is the default
  high-consistency path. `cpp` runs all steps with scop's native kernels
  and is much faster on large gene sets.

- store_networks:

  Whether to keep WT/KO tensor networks in `object@tools`. This only
  controls what is stored; the network ensemble is always built the same
  way.

- store_manifold:

  Whether to keep manifold-alignment coordinates in `object@tools`.

- tool_name:

  Name of the `srt@tools` entry.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

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
#> ℹ [2026-10-11 23:47:07] Run scTenifoldKnk knockout for "Pdx1" using "r" backend
#> ℹ Building 3 gene regulatory networks (200 cells each)
#> ✔ Network construction complete: 3 networks
#> ℹ [X] Tensor: 301 x 301 x 3 (K=3)
#> ✔ [X] CP decomposition complete (norm explained: 39.9%)
#> ℹ Manifold alignment: 301 shared genes, d=2
#> ✔ Manifold alignment complete: 2 dimensions
#> ℹ Computing distances for 301 genes
#> ✔ Differential regulation complete: 5/301 significant genes (FDR < 0.05)
#> ✔ [2026-10-11 23:47:11] scTenifoldKnk results stored in `srt@tools[[scTenifoldKnk]]`

dr <- pancreas_sub@tools$scTenifoldKnk$diffRegulation
head(dr)
#>       gene     distance        Z         FC      p.value        p.adj
#> 1     Pdx1 6.754982e-04 3.963089 7458.81278 0.000000e+00 0.000000e+00
#> 292   Cd81 7.505946e-05 2.316055   92.09419 8.264751e-22 1.243845e-19
#> 256   Myl6 3.581292e-05 1.864745   20.96529 4.676796e-06 4.692386e-04
#> 84   Actg1 3.437486e-05 1.841072   19.31538 1.108104e-05 8.338484e-04
#> 148   Ssr2 2.780289e-05 1.720629   12.63575 3.784388e-04 2.278201e-02
#> 156 Sec61b 2.552481e-05 1.673092   10.64992 1.100757e-03 5.522129e-02

scTenifoldKnkPlot(pancreas_sub, plot_type = "effect")
```
