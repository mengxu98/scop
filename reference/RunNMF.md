# Run NMF

Run NMF

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

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- features:

  Features used instead of a reduction.

- nbes:

  The number of basis vectors (components) to be computed.

- nmf.method:

  The NMF algorithm to be used. Currently supported values are
  `"RcppML"` and `"NMF"`.

- tol:

  The tolerance for convergence (only applicable when nmf.method is
  `"RcppML"`).

- maxit:

  The maximum number of iterations for convergence (only applicable when
  nmf.method is `"RcppML"`).

- rev.nmf:

  Whether to perform reverse NMF (i.e., transpose the input matrix)
  before running the analysis.

- ndims.print:

  The dimensions (number of basis vectors) to print in the output.

- nfeatures.print:

  The number of features to print in the output.

- reduction.name:

  Reduction to be stored in the Seurat object.

- reduction.key:

  The prefix for the column names of the basis vectors.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  Random seed.

- cores:

  The number of threads to be used in older `RcppML` releases that
  expose a global OpenMP thread setter. Newer releases manage their
  thread count internally.

## Examples

``` r
library(Matrix)
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:28:54] Start standard processing workflow...
#> ℹ [2026-09-06 22:28:55] Checking a list of <Seurat>...
#> ! [2026-09-06 22:28:55] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:28:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:28:55] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:28:55] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:28:55] Number of available HVF: 2000
#> ℹ [2026-09-06 22:28:55] Finished check
#> ℹ [2026-09-06 22:28:55] Perform `ScaleData()`
#> ℹ [2026-09-06 22:28:55] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:28:56] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:28:56] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:28:56] Reorder clusters...
#> ℹ [2026-09-06 22:28:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:28:56] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:29:04] Standard processing workflow completed
pancreas_sub <- RunNMF(pancreas_sub)
#> ℹ [2026-09-06 22:29:04] Running NMF...
#> ℹ BE_ 1 
#> ℹ Positive:  Ccnd1, Spp1, Eno1, Rps2, Mif, Ldha, Pebp1, Gapdh, Aldoa, Prdx1 
#> ℹ      Hspe1, Ybx1, Rplp1, Npm1, Mdk, Rpl12, Ptma, Cldn10, Clu, Mgst1 
#> ℹ      Tkt, Cd24a, Slc25a5, Rpl22l1, Krt8, Tpi1, Dbi, Ranbp1, Sox9, Rps12 
#> ℹ Negative:  Tmem108, Poc1a, Epn3, Wipi1, Tmcc3, Fgf12, Plekho1, Tecpr2, Zbtb4, Gm10941 
#> ℹ      Trf, Man1c1, Hmgcs1, Nipal1, Alpl, Larp1b, Kcnip3, Tnr, Gm15915, Cntln 
#> ℹ      Cbfa2t2, Sh2d4a, Bbc3, Megf6, Naaladl2, Fam46d, Hist2h2ac, Slc52a3, Prim2, 1810041L15Rik 
#> ℹ BE_ 2 
#> ℹ Positive:  Spp1, Atp1b1, Sparc, Id2, Vim, Anxa2, Dbi, Mgst1, Gsta3, Mt1 
#> ℹ      Pdzk1ip1, Bicc1, 1700011H14Rik, Acot1, Nudt19, Adamts1, Rpl22l1, Sox9, Clu, Rps2 
#> ℹ      Jun, Ppp1r1b, Rplp1, S100a10, Cldn3, Krt8, Rps12, Lurap1l, Krt18, Cldn6 
#> ℹ Negative:  Aacs, Tmem108, Poc1a, Epn3, B830012L14Rik, Tmcc3, Rfc1, Wsb1, Plekho1, Ppp2r2b 
#> ℹ      Tecpr2, Zbtb4, Haus8, Gm5420, Man1c1, Hmgcs1, Nipal1, Jam3, Tcerg1, Pgap1 
#> ℹ      Alpl, Larp1b, Kcnip3, Tnr, Ptbp3, Gm15915, Cntln, Rbp2, Uchl5, Ahctf1 
#> ℹ BE_ 3 
#> ℹ Positive:  Clu, Mt1, Ckb, Spp1, Mt2, Sparc, Gsta3, Tpi1, Aldoa, Cldn3 
#> ℹ      Cdkn1c, Krt18, 1700011H14Rik, Mgst1, Serpinh1, Prdx1, Tle6, Anxa2, Mif, Rbbp7 
#> ℹ      Muc1, Paics, Ambp, Ldha, Rps2, Myl12a, Litaf, Gapdh, Ptma, Cat 
#> ℹ Negative:  1110002L01Rik, Rpa3, Elovl6, Aacs, Tmem108, Poc1a, Nop56, B830012L14Rik, Tmcc3, Trib1 
#> ℹ      Fgf12, Plekho1, Slc20a1, Ppp2r2b, Zbtb4, Haus8, Gm10941, Trf, Pde4b, Gm5420 
#> ℹ      Man1c1, Nipal1, Jam3, Tcerg1, Pgap1, Snrpa1, Alpl, Larp1b, Kcnip3, Tnr 
#> ℹ BE_ 4 
#> ℹ Positive:  Tpm1, Anxa5, Spp1, Myl12a, Vim, Krt19, Krt18, Ccnd2, Krt8, Sparc 
#> ℹ      Tnfrsf12a, Csrp1, S100a10, Anxa2, Cyr61, Clu, Cldn7, Cdkn1c, Tagln2, Nudt19 
#> ℹ      Myl12b, Jun, Tmsb10, Cldn6, Hsp90aa1, Myl9, Cd24a, Tsc22d1, 1700011H14Rik, Dbi 
#> ℹ Negative:  Elovl6, Aacs, Tmem108, Poc1a, B830012L14Rik, Tmcc3, Rfc1, Plekho1, Slc20a1, Tecpr2 
#> ℹ      Zbtb4, Gm10941, Tenm4, Man1c1, Nipal1, Jam3, Pgap1, Alpl, Kcnip3, Tnr 
#> ℹ      Ptbp3, Gm15915, Cntln, Ocln, Blvrb, Fras1, Rbp2, Cbfa2t2, Ptgfrn, Ugt8a 
#> ℹ BE_ 5 
#> ℹ Positive:  2810417H13Rik, Hmgb2, Rrm2, Dut, Pcna, Tipin, Lig1, Tuba1b, Mcm5, H2afz 
#> ℹ      Tk1, Gmnn, Dek, Tyms, Mcm3, Rfc2, Ran, Orc6, Tubb5, Ranbp1 
#> ℹ      Mcm6, Uhrf1, Srsf2, Rfc3, Mcm7, Gins2, Rpa2, Siva1, Dnajc9, Smc2 
#> ℹ Negative:  1110002L01Rik, Aacs, Wipi1, B830012L14Rik, Tmcc3, Trib1, Fgf12, Plekho1, Ppp2r2b, Lama1 
#> ℹ      Tecpr2, Zbtb4, P4ha3, Tenm4, Trf, Gm5420, Man1c1, Nipal1, Jam3, Alpl 
#> ℹ      Mgst2, Kcnip3, Tnr, Gm15915, Cbfa2t2, Sh2d4a, Bbc3, Surf4, Ano6, Megf6 
#> ✔ [2026-09-06 22:29:14] NMF compute completed
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "nmf"
)


FeatureDimPlot(
  pancreas_sub,
  features = c("BE_1", "BE_2", "BE_3"),
  reduction = "UMAP",
  palette = "RdBu",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  theme_use = "theme_blank"
)


ht_cells <- NMFHeatmap(
  pancreas_sub,
  plot_type = "cells",
  cell_annotation = "CellType"
)
#> ℹ [2026-09-06 22:29:14] `NMFHeatmap()` input: 1000 cells x 50 NMF dimensions. Computing a 1000 x 1000 similarity matrix (~0.01 GiB dense numeric matrix).
#> ℹ [2026-09-06 22:29:14] Ordering `NMFHeatmap()` rows and columns ...
#> ℹ [2026-09-06 22:29:14] Building ComplexHeatmap object for `NMFHeatmap()` ...
#> ℹ [2026-09-06 22:29:14] Calculating `NMFHeatmap()` render size ...
#> ℹ [2026-09-06 22:29:14] Drawing `NMFHeatmap()`; this can take time for large similarity matrices ...
#> ℹ [2026-09-06 22:29:16] Assembling `NMFHeatmap()` plot object ...
ht_cells$plot


ht_features <- NMFHeatmap(
  pancreas_sub,
  plot_type = "features"
)
#> ℹ [2026-09-06 22:29:16] `NMFHeatmap()` input: 2000 features x 50 NMF dimensions. Computing a 2000 x 2000 similarity matrix (~0.03 GiB dense numeric matrix).
#> ℹ [2026-09-06 22:29:16] Ordering `NMFHeatmap()` rows and columns ...
#> ℹ [2026-09-06 22:29:17] Building ComplexHeatmap object for `NMFHeatmap()` ...
#> ℹ [2026-09-06 22:29:17] Calculating `NMFHeatmap()` render size ...
#> ℹ [2026-09-06 22:29:17] Drawing `NMFHeatmap()`; this can take time for large similarity matrices ...
#> ℹ [2026-09-06 22:29:22] Assembling `NMFHeatmap()` plot object ...
ht_features$plot
```
