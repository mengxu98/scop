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
  cores = 0,
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
  cores = 0,
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
  cores = 0,
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
  cores = 0,
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

- cores:

  The number of threads to be used in `RcppML` functions that are
  parallelized with `OpenMP`. If `0`, the number of threads will be
  automatically determined by
  [`RcppML::setRcppMLthreads()`](https://rdrr.io/pkg/RcppML/man/setRcppMLthreads.html).
  Default is `0`.

## Examples

``` r
library(Matrix)
#> 
#> Attaching package: ‘Matrix’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     expand
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-27 08:22:02] Start standard scop workflow...
#> ℹ [2026-01-27 08:22:02] Checking a list of <Seurat>...
#> ! [2026-01-27 08:22:03] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-27 08:22:03] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 08:22:05] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-27 08:22:05] Use the separate HVF from srt_list
#> ℹ [2026-01-27 08:22:05] Number of available HVF: 2000
#> ℹ [2026-01-27 08:22:06] Finished check
#> ℹ [2026-01-27 08:22:06] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-27 08:22:06] Perform pca linear dimension reduction
#> ℹ [2026-01-27 08:22:07] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-27 08:22:07] Reorder clusters...
#> ℹ [2026-01-27 08:22:07] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-27 08:22:07] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-27 08:22:12] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-27 08:22:17] Run scop standard workflow completed
pancreas_sub <- RunNMF(pancreas_sub)
#> ℹ [2026-01-27 08:22:17] Running NMF...
#> ℹ BE_ 1 
#> ℹ Positive:  Ccnd1, Spp1, Mdk, Rps2, Ldha, Pebp1, Cd24a, Dlk1, Krt8, Mgst1 
#> ℹ      Clu, Gapdh, Eno1, Prdx1, Cldn10, Mif, Cldn7, Npm1, Dbi, Vim 
#> ℹ      Sox9, Rpl12, Aldh1b1, Rplp1, Wfdc2, Krt18, Tkt, Aldoa, Hspe1, Ptma 
#> ℹ Negative:  Tmem108, Poc1a, Epn3, Wipi1, Tmcc3, Nhsl1, Fgf12, Plekho1, Tecpr2, Zbtb4 
#> ℹ      Gm10941, Trf, Man1c1, Hmgcs1, Nipal1, Jam3, Pgap1, Alpl, Kcnip3, Tnr 
#> ℹ      Gm15915, Rbp2, Cbfa2t2, Sh2d4a, Bbc3, Megf6, Naaladl2, Fam46d, Hist2h2ac, Tox2 
#> ℹ BE_ 2 
#> ℹ Positive:  Spp1, Gsta3, Sparc, Vim, Atp1b1, Mt1, Dbi, Anxa2, Rps2, Id2 
#> ℹ      Rpl22l1, Rplp1, Mgst1, Clu, Sox9, Cldn6, Mdk, Pdzk1ip1, Bicc1, 1700011H14Rik 
#> ℹ      Rps12, S100a10, Cldn3, Rpl36a, Ppp1r1b, Adamts1, Serpinh1, Mt2, Ifitm2, Rpl39 
#> ℹ Negative:  Rpa3, Aacs, Tmem108, Poc1a, Epn3, Wipi1, B830012L14Rik, Tmcc3, Wsb1, Plekho1 
#> ℹ      Ppp2r2b, Tecpr2, Zbtb4, Haus8, Trf, Gm5420, Man1c1, Hmgcs1, Nipal1, Jam3 
#> ℹ      Tcerg1, Pgap1, Snrpa1, Alpl, Larp1b, Kcnip3, Tnr, Lsm12, Ptbp3, Gm15915 
#> ℹ BE_ 3 
#> ℹ Positive:  Cck, Mdk, Gadd45a, Neurog3, Selm, Sox4, Btbd17, Tmsb4x, Btg2, Cldn6 
#> ℹ      Cotl1, Ptma, Jun, Ppp1r14a, Rps2, Ifitm2, Neurod2, Igfbpl1, Gnas, Krt7 
#> ℹ      Nkx6-1, Aplp1, Ppp3ca, Lrpap1, Rplp1, Hn1, Rps12, Mfng, BC023829, Smarcd2 
#> ℹ Negative:  Elovl6, Tmem108, Poc1a, Epn3, Nop56, Wipi1, B830012L14Rik, Rrp15, Rfc1, Fgf12 
#> ℹ      Slc20a1, Ppp2r2b, Lama1, Tecpr2, Zbtb4, Eif1ax, Fam162a, P4ha3, Gm10941, Tenm4 
#> ℹ      Pde4b, Gm5420, Man1c1, Hmgcs1, Pgap1, Mgst2, Larp1b, Kcnip3, Tnr, Lsm12 
#> ℹ BE_ 4 
#> ℹ Positive:  Spp1, Cyr61, Krt18, Tpm1, Krt8, Myl12a, Vim, Jun, Anxa5, Tnfrsf12a 
#> ℹ      Csrp1, Sparc, Cldn7, Nudt19, Anxa2, Clu, Myl9, Atp1b1, Cldn3, Tagln2 
#> ℹ      S100a10, 1700011H14Rik, Cd24a, Rps2, Dbi, Id2, Lurap1l, Rplp1, Myl12b, Klf6 
#> ℹ Negative:  Rpa3, Elovl6, Aacs, Tmem108, Poc1a, Tmcc3, Rfc1, Plekho1, Slc20a1, Ppp2r2b 
#> ℹ      Lama1, Tecpr2, Gm10941, Tenm4, Pde4b, Man1c1, Nipal1, Jam3, Pgap1, Alpl 
#> ℹ      Mgst2, Kcnip3, Tnr, Ptbp3, Gm15915, Cntln, Ocln, Fras1, Rbp2, Cbfa2t2 
#> ℹ BE_ 5 
#> ℹ Positive:  2810417H13Rik, Rrm2, Hmgb2, Dut, Pcna, Lig1, H2afz, Tipin, Tuba1b, Tk1 
#> ℹ      Mcm5, Dek, Tyms, Gmnn, Ran, Tubb5, Rfc2, Srsf2, Ranbp1, Orc6 
#> ℹ      Mcm3, Uhrf1, Gins2, Dnajc9, Mcm6, Siva1, Rfc3, Mcm7, Rpa2, Ptma 
#> ℹ Negative:  1110002L01Rik, Aacs, Wipi1, B830012L14Rik, Tmcc3, Trib1, Fgf12, Plekho1, Ppp2r2b, Lama1 
#> ℹ      Tenm4, Trf, Gm5420, Man1c1, Jam3, Mgst2, Kcnip3, Tnr, Gm15915, Cbfa2t2 
#> ℹ      Sh2d4a, Bbc3, Fkbp9, Ano6, Prkcb, Megf6, Fam46d, Slc52a3, Ankrd2, Tox2 
#> ✔ [2026-01-27 08:22:21] NMF compute completed
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "nmf"
)
```
