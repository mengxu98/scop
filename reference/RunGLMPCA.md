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
#> ℹ [2026-06-28 10:28:56] Start standard processing workflow...
#> ℹ [2026-06-28 10:28:56] Checking a list of <Seurat>...
#> ! [2026-06-28 10:28:56] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-28 10:28:56] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 10:28:56] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-28 10:28:56] Use the separate HVF from `srt_list`
#> ℹ [2026-06-28 10:28:57] Number of available HVF: 2000
#> ℹ [2026-06-28 10:28:57] Finished check
#> ℹ [2026-06-28 10:28:57] Perform `ScaleData()`
#> ℹ [2026-06-28 10:28:57] Perform pca linear dimension reduction
#> ℹ [2026-06-28 10:28:57] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-28 10:28:58] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-28 10:28:58] Reorder clusters...
#> ℹ [2026-06-28 10:28:58] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-28 10:28:58] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-28 10:29:05] Standard processing workflow completed
pancreas_sub <- RunGLMPCA(pancreas_sub)
#> ℹ GLMPC_ 1 
#> ℹ Positive:  Cartpt, Barx2, Cdkn2b, Gip, Aard, Pax6os1, Prl, Ucn3, Ptger3, G6pc2 
#> ℹ      Pcdh8, Pappa2, Slc30a8, Dusp26, 1700001C02Rik, Spock1, Acsbg1, Vgf, Bace2, Pcp4 
#> ℹ      Kctd8, Ceacam10, 3930402G23Rik, 1700015F17Rik, 4930426D05Rik, Fbln5, Ctxn2, Rnf138rt1, Scn9a, 2410021H03Rik 
#> ℹ Negative:  Col1a1, AA986860, Gm6878, P2ry2, Il18, Ctgf, Guca2b, Scx, S100b, Plscr2 
#> ℹ      Pmp22, Sp140, Anxa9, Adamts16, Prickle2, Smoc2, Islr, 1110002O04Rik, Tmem178, Tns1 
#> ℹ      Cryab, Gsta3, Dcdc2a, Adgrg6, Tmem171, Isg15, Grin3a, Cxcl12, Gm20649, Apcs 
#> ℹ GLMPC_ 2 
#> ℹ Positive:  Guca2a, Gsg1l, Ifitm1, Adra2c, Fam71b, Laptm5, P2ry14, Tmprss6, Bhlhe22, Bhmt2 
#> ℹ      Dlgap1, Jakmip3, Bcas1os1, Slc4a1, Tff3, Nhs, Alb, Ins1, Wnt3, Nkpd1 
#> ℹ      Pou3f1, Rhbg, Tspear, Fam124a, Pdcd1, Lrrn2, Krtap17-1, Kcnq4, Slc39a2, Sema3g 
#> ℹ Negative:  Sst, Aif1, Fam198b, Klhl14, Tac1, Edn1, Kcne2, RP23-58K20.3, Gtf2ird2, Platr22 
#> ℹ      Ctsk, Elovl4, Hoxb4, Slit2, Col25a1, Aspm, Nkain4, Igfbp3, Kif2c, Plscr2 
#> ℹ      Pkd2l1, Nlgn1, Gast, 4430402I18Rik, Cdc25c, D7Ertd443e, Bmp2, Pif1, Slc4a10, Tstd1 
#> ℹ GLMPC_ 3 
#> ℹ Positive:  Fgb, Col1a2, Klhl14, Gcg, Lgr5, Doc2a, Calb1, Gad2, Gast, Ryr3 
#> ℹ      Hist1h4a, Gm11744, 4930426D05Rik, Ctsk, Tstd1, 4930539E08Rik, Cbln4, Crygn, Dkk2, Col23a1 
#> ℹ      P2ry2, Tac1, Pou6f2, Tgfb2, Arhgap36, Ctgf, RP23-428N8.3, BC043934, Sp140, Smpx 
#> ℹ Negative:  Pif1, Msx1, Ppp1r17, Gtse1, Mmel1, Pf4, Mxd3, Depdc1a, Gm42984, Ccnb1 
#> ℹ      Nusap1, Hmmr, Igfbp3, Slfn2, Shox2, Sapcd2, Cbln1, Parpbp, Plk1, Cenpf 
#> ℹ      Cnrip1, Aurka, Fgf8, Icosl, Cdx2, RP23-4H17.3, Cdc20, Kif20a, Cdc25c, Nhlh1 
#> ℹ GLMPC_ 4 
#> ℹ Positive:  Gast, Bmp2, Lmod3, Ifit1bl1, 1500035N22Rik, Rerg, Ngf, RP23-385E22.2, Zfp97, Tnfaip8l3 
#> ℹ      Fam46d, Gm29440, Zcchc12, Avp, Gm13375, Mboat4, Tstd1, Arhgap22, Tmem255b, Sp140 
#> ℹ      Pkd2l1, Cypt3, Nkx6-3, Nipal3, Tox2, 1110002O04Rik, Elovl4, Snai2, Rasgrp3, Trp53cor1 
#> ℹ Negative:  Crygn, P2ry1, Il1r2, Npy, Sp5, RP23-58K20.3, Aif1, Tmem215, Adam32, Lgr5 
#> ℹ      Lrrc6, Adgrf5, Klhl14, Pid1, Arhgap36, 1700024G13Rik, Col25a1, Ins2, Gad2, Sytl4 
#> ℹ      Srgn, Slfn2, Iapp, Gm11789, Gm38112, Dkk2, Doc2a, Nnat, Ins1, Bace2 
#> ℹ GLMPC_ 5 
#> ℹ Positive:  Il1r2, Fcgr3, Srgn, Coro1a, Slc4a10, Anxa1, Lrrc6, Cd37, Itgb7, Tyrobp 
#> ℹ      Tmem100, 4933440M02Rik, 4430402I18Rik, Gtf2ird2, Kcnk10, Gm11636, Lmx1a, Sytl2, Lst1, Tssk6 
#> ℹ      Agmo, Ppfibp2, Slfn2, Gm15895, Ncf2, Ltb, Gm6410, Lgr5, Ryr3, Chgb 
#> ℹ Negative:  Sparcl1, Col3a1, Dkk2, Galnt16, Cbln4, Kcnj8, Tex36, Col6a1, Islr, Gpr6 
#> ℹ      Col1a2, Col1a1, Col23a1, Gm933, Irs4, Ptpro, Olfml3, Spock1, Kcne2, Lrrtm3 
#> ℹ      Gm15640, RP23-428N8.3, Colec12, Ghrl, Gm17455, Th, RP23-172P1.4, Ceacam10, Aunip, Clspn 
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "glmpca"
)
```
