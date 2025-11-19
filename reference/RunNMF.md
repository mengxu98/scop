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

  The assay to be used for the analysis. Default is `NULL`.

- layer:

  The layer to be used for the analysis. Default is `"data"`.

- features:

  The features to be used for the analysis. Default is `NULL`, which
  uses all variable features.

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

  The random seed to be used. Default is `11`.

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
#> ℹ [2025-11-19 14:58:58] Start standard scop workflow...
#> ℹ [2025-11-19 14:58:59] Checking a list of <Seurat> object...
#> ! [2025-11-19 14:58:59] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:58:59] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:59:01] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:59:02] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:59:02] Number of available HVF: 2000
#> ℹ [2025-11-19 14:59:02] Finished check
#> ℹ [2025-11-19 14:59:02] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 14:59:02] Perform pca linear dimension reduction
#> StandardPC_ 1 
#> Positive:  Aplp1, Cpe, Gnas, Fam183b, Map1b, Hmgn3, Pcsk1n, Chga, Tuba1a, Bex2 
#>     Syt13, Isl1, 1700086L19Rik, Pax6, Chgb, Scgn, Rbp4, Scg3, Gch1, Camk2n1 
#>     Cryba2, Pcsk2, Pyy, Tspan7, Mafb, Hist3h2ba, Dbpht2, Abcc8, Rap1b, Slc38a5 
#> Negative:  Spp1, Anxa2, Sparc, Dbi, 1700011H14Rik, Wfdc2, Gsta3, Adamts1, Clu, Mgst1 
#>     Bicc1, Ldha, Vim, Cldn3, Cyr61, Rps2, Mt1, Ptn, Phgdh, Nudt19 
#>     Smtnl2, Smco4, Habp2, Mt2, Col18a1, Rpl12, Galk1, Cldn10, Acot1, Ccnd1 
#> StandardPC_ 2 
#> Positive:  Rbp4, Tagln2, Tuba1b, Fkbp2, Pyy, Pcsk2, Iapp, Tmem27, Meis2, Tubb4b 
#>     Pcsk1n, Dbpht2, Rap1b, Dynll1, Tubb2a, Sdf2l1, Scgn, 1700086L19Rik, Scg2, Abcc8 
#>     Atp1b1, Hspa5, Fam183b, Papss2, Slc38a5, Scg3, Mageh1, Tspan7, Ppp1r1a, Ociad2 
#> Negative:  Neurog3, Btbd17, Gadd45a, Ppp1r14a, Neurod2, Sox4, Smarcd2, Mdk, Pax4, Btg2 
#>     Sult2b1, Hes6, Grasp, Igfbpl1, Gpx2, Cbfa2t3, Foxa3, Shf, Mfng, Tmsb4x 
#>     Amotl2, Gdpd1, Cdc14b, Epb42, Rcor2, Cotl1, Upk3bl, Rbfox3, Cldn6, Cer1 
#> StandardPC_ 3 
#> Positive:  Nusap1, Top2a, Birc5, Aurkb, Cdca8, Pbk, Mki67, Tpx2, Plk1, Ccnb1 
#>     2810417H13Rik, Incenp, Cenpf, Ccna2, Prc1, Racgap1, Cdk1, Aurka, Cdca3, Hmmr 
#>     Spc24, Kif23, Sgol1, Cenpe, Cdc20, Hist1h1b, Cdca2, Mxd3, Kif22, Ska1 
#> Negative:  Anxa5, Pdzk1ip1, Acot1, Tpm1, Anxa2, Dcdc2a, Capg, Sparc, Ttr, Pamr1 
#>     Clu, Cxcl12, Ndrg2, Hnf1aos1, Gas6, Gsta3, Krt18, Ces1d, Atp1b1, Muc1 
#>     Hhex, Acadm, Spp1, Enpp2, Bcl2l14, Sat1, Smtnl2, 1700011H14Rik, Tgm2, Fam159a 
#> StandardPC_ 4 
#> Positive:  Glud1, Tm4sf4, Akr1c19, Cldn4, Runx1t1, Fev, Pou3f4, Gm43861, Pgrmc1, Arx 
#>     Cd200, Lrpprc, Hmgn3, Ppp1r14c, Pam, Etv1, Tsc22d1, Slc25a5, Akap17b, Pgf 
#>     Fam43a, Emb, Jun, Krt8, Dnajc12, Mid1ip1, Ids, Rgs17, Uchl1, Alcam 
#> Negative:  Ins2, Ins1, Ppp1r1a, Nnat, Calr, Sytl4, Sdf2l1, Iapp, Pdia6, Mapt 
#>     G6pc2, C2cd4b, Npy, Gng12, P2ry1, Ero1lb, Adra2a, Papss2, Arhgap36, Fam151a 
#>     Dlk1, Creld2, Gip, Tmem215, Gm27033, Cntfr, Prss53, C2cd4a, Lyve1, Ociad2 
#> StandardPC_ 5 
#> Positive:  Pdx1, Nkx6-1, Npepl1, Cldn4, Cryba2, Fev, Jun, Chgb, Gng12, Adra2a 
#>     Mnx1, Sytl4, Pdk3, Gm27033, Nnat, Chga, Ins2, 1110012L19Rik, Enho, Krt7 
#>     Mlxipl, Tmsb10, Flrt1, Pax4, Tubb3, Prrg2, Gars, Frzb, BC023829, Gm2694 
#> Negative:  Irx2, Irx1, Gcg, Ctxn2, Tmem27, Ctsz, Tmsb15l, Nap1l5, Pou6f2, Gria2 
#>     Ghrl, Peg10, Smarca1, Arx, Lrpap1, Rgs4, Ttr, Gast, Tmsb15b2, Serpina1b 
#>     Slc16a10, Wnk3, Ly6e, Auts2, Sct, Arg1, Dusp10, Sphkap, Dock11, Edn3 
#> ℹ [2025-11-19 14:59:03] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 14:59:04] Reorder clusters...
#> ℹ [2025-11-19 14:59:04] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 14:59:04] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:59:04] UMAP will return its model
#> ℹ [2025-11-19 14:59:08] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:59:08] UMAP will return its model
#> ✔ [2025-11-19 14:59:13] Run scop standard workflow done
pancreas_sub <- RunNMF(pancreas_sub)
#> ✔ [2025-11-19 14:59:13] zdebruine/RcppML installed successfully
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
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "nmf"
)
```
