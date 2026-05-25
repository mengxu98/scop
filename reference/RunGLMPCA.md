# Run generalized principal components analysis (GLMPCA)

Run generalized principal components analysis (GLMPCA)

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

  An object. Can be a Seurat object, an assay object, or a matrix-like
  object.

- ...:

  Additional arguments to be passed to the
  [glmpca::glmpca](https://rdrr.io/pkg/glmpca/man/glmpca.html) function.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Which layer to use. Default is `data`.

- features:

  A character vector of features to use. Default is `NULL`.

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

  Random seed for reproducibility. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-25 08:17:02] Start standard processing workflow...
#> ℹ [2026-05-25 08:17:02] Checking a list of <Seurat>...
#> ! [2026-05-25 08:17:02] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 08:17:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 08:17:04] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 08:17:04] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 08:17:05] Number of available HVF: 2000
#> ℹ [2026-05-25 08:17:05] Finished check
#> ℹ [2026-05-25 08:17:05] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 08:17:05] Perform pca linear dimension reduction
#> ℹ [2026-05-25 08:17:05] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-25 08:17:06] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 08:17:06] Reorder clusters...
#> ℹ [2026-05-25 08:17:06] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 08:17:06] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 08:17:06] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-25 08:17:11] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-25 08:17:15] Standard processing workflow completed
pancreas_sub <- RunGLMPCA(pancreas_sub)
#> ℹ GLMPC_ 1 
#> ℹ Positive:  Barx2, Cartpt, Gm3448, Ptger3, Kng2, Cypt3, Il1r2, 3930402G23Rik, Mesp1, Dusp26 
#> ℹ      Gad1, Cdkn2b, Platr26, Gm933, Aard, Ky, Prl, Lrrc6, 4930426D05Rik, Ucn3 
#> ℹ      RP23-385E22.2, Pcdh8, Slc38a5, Ctxn2, A930017K11Rik, Gm6086, 1700001C02Rik, Pax6os1, Ifit1bl1, Msx1 
#> ℹ Negative:  Col6a1, Hoxb4, Sparcl1, Col23a1, Ctgf, Zfp385b, Galnt16, Col1a1, Isg15, Col1a2 
#> ℹ      Pmp22, Islr, Platr22, P2ry2, Plscr2, Gm26633, Gm15640, Gm6878, Kcnj8, Col3a1 
#> ℹ      Smpx, A730098A19Rik, Krt20, Ctsk, Tnfrsf19, Olfml2a, Cxcl10, Peg12, Cdkn2c, Pgr 
#> ℹ GLMPC_ 2 
#> ℹ Positive:  Lmx1a, Ncf2, Cmklr1, Slfn2, Tgm7, Pthlh, Sema3g, Gm15567, 1520401A03Rik, Nhlh1 
#> ℹ      Gm16140, Adgrb1, Epb42, Lipg, Neurog3, Prom2, Slc52a3, 1700128E19Rik, Gm8773, Bhlhe22 
#> ℹ      Siglece, Crlf1, Neurod2, Notum, Eya2, Fn3krp, Rasgrp3, Wnt3, Tmem114, Acot11 
#> ℹ Negative:  Sst, Dkk2, Klhl14, RP23-428N8.3, Aif1, 4930539E08Rik, Tnni3, Fgb, Col1a2, Ctsk 
#> ℹ      Ctgf, Col23a1, Tac1, Col25a1, Crygn, Lgr5, Ppy, Col6a1, Gm26633, Hoxb4 
#> ℹ      Otc, Zfp385b, Pyy, Nov, Tstd1, Sparcl1, 4930426D05Rik, Gad2, Ceacam10, Kcne2 
#> ℹ GLMPC_ 3 
#> ℹ Positive:  Gtf2ird2, Fam198b, Tac1, Cbln1, Ctsk, 1700015F17Rik, Sst, Slc4a10, 4933440M02Rik, Lmx1a 
#> ℹ      Tstd1, Gm10382, RP23-58K20.3, Pkd2l1, Ceacam10, Klhl14, Igfbp3, 4430402I18Rik, Lgr5, Cdc25c 
#> ℹ      Kif2c, 4930539E08Rik, Slfn2, Icosl, Gm15895, Cenpf, Ppfibp2, Cpne5, Fcgr3, Aspm 
#> ℹ Negative:  Kcnj8, Cbln4, Col3a1, Islr, Sparcl1, Gm15640, Col1a1, Col23a1, Ghrl, Col5a1 
#> ℹ      Tex36, Galnt16, Guca2a, Col6a1, Col1a2, Sapcd1, Foxd3, L1td1, Npy, Gm38112 
#> ℹ      Irs4, Tagln, Anxa1, Olfml3, Pid1, Cldn18, Ptpro, Tmtc1, Lrrtm3, Gm17455 
#> ℹ GLMPC_ 4 
#> ℹ Positive:  Npy, Gm38112, Aif1, Cldn18, Dlgap1, Rac2, Ins2, Cbln4, Syndig1l, Sparcl1 
#> ℹ      Hist1h1a, Gm11789, Ins1, P2ry14, 1700024G13Rik, Pif1, Gm15640, Tmem215, Slfn2, Galnt16 
#> ℹ      P2ry1, Pf4, Iqgap3, Col25a1, Gip, RP23-172P1.4, Nnat, RP23-58K20.3, Sp5, Sst 
#> ℹ Negative:  Tstd1, Sp140, D7Ertd443e, Gast, Tnfaip8l3, Lmod3, Fam46d, Ctsk, 1500035N22Rik, Rerg 
#> ℹ      Guca2a, Anxa1, Gm29440, Arhgap22, Pou6f2, Calb1, Snai2, Smpx, Plbd1, Gcg 
#> ℹ      Irs4, Bmp2, RP23-385E22.2, Nxph1, Zfp97, Ngf, Gm13375, Tagln, Sycp3, Bhlhe23 
#> ℹ GLMPC_ 5 
#> ℹ Positive:  Srgn, Kcne2, Cxcl10, Sst, Fcgr3, Aif1, Rac2, Ccl20, P2ry14, Tex36 
#> ℹ      Ltb, Elovl4, Krt17, Anxa1, Plaur, Tyrobp, Lst1, Itgb7, Tmem100, Gm933 
#> ℹ      Cpa3, Gm17455, Kcnk10, Arhgap22, Cyp11a1, Lrrtm3, Cd37, Cbln4, Gpr6, Ifitm1 
#> ℹ Negative:  Gad2, 4933440M02Rik, Sparcl1, Hoxb4, Gcg, Col1a2, Oasl2, Nhs, Gsg1l, Sp140 
#> ℹ      Ins1, Islr, Gm6878, Calb1, Galnt16, Col6a1, Rspo1, Lmx1a, Ryr3, Tmem119 
#> ℹ      1110002O04Rik, 4930539E08Rik, Fam71b, Cdkn2c, Col1a1, Pmp22, Lgr5, Pid1, Pkd2l1, Dlgap1 
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "glmpca"
)
```
