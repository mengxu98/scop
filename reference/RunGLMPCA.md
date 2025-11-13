# Run GLMPCA (generalized version of principal components analysis)

Run GLMPCA (generalized version of principal components analysis)

## Usage

``` r
RunGLMPCA(object, ...)

# S3 method for class 'Seurat'
RunGLMPCA(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "glmpca",
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay'
RunGLMPCA(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay5'
RunGLMPCA(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunGLMPCA(
  object,
  assay = NULL,
  layer = "counts",
  features = NULL,
  L = 5,
  fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
  rev.gmlpca = FALSE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "GLMPC_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, an assay object, or a
  matrix-like object.

- ...:

  Additional arguments to be passed to the
  [glmpca::glmpca](https://rdrr.io/pkg/glmpca/man/glmpca.html) function.

- assay:

  The assay to be used for the analysis. Default is `NULL`.

- layer:

  The layer to be used for the analysis. Default is `"counts"`.

- features:

  The features to be used for the analysis. Default is `NULL`, which
  uses all variable features.

- L:

  The number of components to be computed. Default is `5`.

- fam:

  The family of the generalized linear model to be used. Currently
  supported values are `"poi"`, `"nb"`, `"nb2"`, `"binom"`, `"mult"`,
  and `"bern"`. Default is `"poi"`.

- rev.gmlpca:

  Whether to perform reverse GLMPCA (i.e., transpose the input matrix)
  before running the analysis. Default is `FALSE`.

- ndims.print:

  The dimensions (number of components) to print in the output. Default
  is `1:5`.

- nfeatures.print:

  The number of features to print in the output. Default is `30`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"glmpca"`.

- reduction.key:

  The prefix for the column names of the basis vectors. Default is
  `"GLMPC_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  The random seed to be used. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 12:23:04] Start standard scop workflow...
#> ℹ [2025-11-13 12:23:05] Checking a list of <Seurat> object...
#> ! [2025-11-13 12:23:05] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 12:23:05] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:23:07] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:23:07] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:23:08] Number of available HVF: 2000
#> ℹ [2025-11-13 12:23:08] Finished check
#> ℹ [2025-11-13 12:23:08] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:23:08] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 12:23:09] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:23:09] Reorder clusters...
#> ℹ [2025-11-13 12:23:09] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:23:09] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:23:09] UMAP will return its model
#> ℹ [2025-11-13 12:23:14] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:23:14] UMAP will return its model
#> ✔ [2025-11-13 12:23:19] Run scop standard workflow done
pancreas_sub <- RunGLMPCA(pancreas_sub)
#> ◌ [2025-11-13 12:23:19] Installing: glmpca...
#>  
#> → Will install 1 package.
#> → The package (0 B) is cached.
#> + glmpca   0.2.0 
#>   
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Installed glmpca 0.2.0  (1s)
#> ✔ 1 pkg + 1 dep: kept 1, added 1 [2s]
#> ✔ [2025-11-13 12:23:21] glmpca installed successfully
#> ℹ GLMPC_ 1 
#> ℹ Positive:  Barx2, Cartpt, Ptger3, Gm3448, Cypt3, Kng2, Gad1, 3930402G23Rik, Mesp1, Il1r2 
#> ℹ      Dusp26, Cdkn2b, Aard, Platr26, Prl, Ucn3, Slc38a5, Pcdh8, Lrrc6, 4930426D05Rik 
#> ℹ      Ky, Ctxn2, Pax6os1, Msx1, A930017K11Rik, 1700001C02Rik, Ifit1bl1, RP23-385E22.2, Gm933, Gm6086 
#> ℹ Negative:  Col6a1, Sparcl1, Col23a1, Hoxb4, Col1a1, Galnt16, Ctgf, Zfp385b, Col1a2, Islr 
#> ℹ      Isg15, Pmp22, P2ry2, Platr22, Plscr2, Gm6878, Gm26633, Kcnj8, Smpx, A730098A19Rik 
#> ℹ      Col3a1, Gm15640, Aif1, Ctsk, Cdkn2c, Peg12, Olfml2a, Krt20, Tnni3, Tnfrsf19 
#> ℹ GLMPC_ 2 
#> ℹ Positive:  Ncf2, Lmx1a, Cmklr1, Tgm7, Pthlh, Nhlh1, 1520401A03Rik, Gm15567, Lipg, Epb42 
#> ℹ      Adgrb1, Sema3g, Gm16140, Slc52a3, Neurog3, 1700128E19Rik, Crlf1, Notum, Acot11, Siglece 
#> ℹ      Eya2, Neurod2, Laptm5, Prom2, Gm8773, Fgf18, Wnt3, Bhlhe22, Megf11, Rasgrp3 
#> ℹ Negative:  Sst, Dkk2, Klhl14, Aif1, RP23-428N8.3, 4930539E08Rik, Ctsk, Fgb, Tnni3, Col1a2 
#> ℹ      Ctgf, Tac1, Tstd1, Col25a1, Col6a1, Col23a1, Lgr5, Crygn, Sparcl1, Nov 
#> ℹ      Ppy, 4930426D05Rik, Pyy, Otc, Hoxb4, Zfp385b, Sp140, Cbln4, Ceacam10, Gm26633 
#> ℹ GLMPC_ 3 
#> ℹ Positive:  Gtf2ird2, Fam198b, Tac1, 4933440M02Rik, Tstd1, 1700015F17Rik, Sst, Cbln1, Pkd2l1, Lmx1a 
#> ℹ      Slc4a10, Ctsk, RP23-58K20.3, 4430402I18Rik, Gm10382, Igfbp3, Klhl14, Gm15895, Fcgr3, 4930539E08Rik 
#> ℹ      Lst1, Cdc25c, Kif2c, Lgr5, Cenpf, Ankrd1, Icosl, Slfn2, Srgn, D7Ertd443e 
#> ℹ Negative:  Kcnj8, Cbln4, Col3a1, Sparcl1, Gm15640, Islr, Col5a1, Col1a1, Ghrl, Tex36 
#> ℹ      Npy, Guca2a, Col1a2, Galnt16, Col23a1, Sapcd1, Foxd3, Col6a1, Tmem119, L1td1 
#> ℹ      Gm38112, Ptpro, Olfml3, Tmtc1, Irs4, Tagln, Lrrtm3, Cxcl16, Lsp1, Nid1 
#> ℹ GLMPC_ 4 
#> ℹ Positive:  Npy, Aif1, Gm38112, Cldn18, Dlgap1, Rac2, Ins2, Gm11789, Ins1, Syndig1l 
#> ℹ      Sp5, Hist1h1a, Gm15640, Col25a1, P2ry14, Slfn2, Tmem215, Sst, Gip, Cbln4 
#> ℹ      1700024G13Rik, Iqgap3, P2ry1, Pif1, Gm933, Adam32, Ifitm1, RP23-58K20.3, Pf4, Hist1h1b 
#> ℹ Negative:  Sp140, Gast, Tstd1, Tnfaip8l3, D7Ertd443e, Fam46d, Lmod3, Rerg, 1500035N22Rik, Guca2a 
#> ℹ      Arhgap22, Pou6f2, Gm29440, RP23-385E22.2, Snai2, Irs4, Plbd1, Calb1, Ctsk, Nrn1 
#> ℹ      Gcg, Anxa1, Nxph1, Smpx, Ngf, Bhlhe23, Gm13375, Oasl2, 1110002O04Rik, Bmp2 
#> ℹ GLMPC_ 5 
#> ℹ Positive:  Cxcl10, Rac2, Srgn, Kcne2, Fcgr3, Sst, Ccl20, Anxa1, Aif1, Krt17 
#> ℹ      Tyrobp, Elovl4, Ltb, Plaur, Tex36, P2ry14, Itgb7, Lst1, Tmem100, Cyp11a1 
#> ℹ      Tnni3, Gm933, Cpa3, Gm17455, Cxcl16, Bmp2, Cd37, Arhgap22, Lrrtm3, Gpr6 
#> ℹ Negative:  Gad2, Sparcl1, 4933440M02Rik, Gcg, Galnt16, Oasl2, Islr, Col6a1, Hoxb4, Col1a2 
#> ℹ      Gsg1l, Nhs, Gm6878, Calb1, Col1a1, Lgr5, Ryr3, Ins1, Rspo1, Cdkn2c 
#> ℹ      4930539E08Rik, Tmem119, Pid1, Pmp22, Lmx1a, Ska3, Pkd2l1, Aspm, Tmem255b, Fam71b 
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "glmpca"
)
```
