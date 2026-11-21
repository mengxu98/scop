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

  Arguments passed to other functions.

- assay:

  Assay to use. `NULL` uses the default assay.

- layer:

  Assay layer to use.

- features:

  Features used instead of a reduction.

- L:

  The number of components to be computed.

- fam:

  The family of the generalized linear model to be used. Currently
  supported values are `"poi"`, `"nb"`, `"nb2"`, `"binom"`, `"mult"`,
  and `"bern"`.

- rev.gmlpca:

  Whether to perform reverse GLMPCA (i.e., transpose the input matrix)
  before running the analysis.

- ndims.print:

  The dimensions (number of components) to print in the output.

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

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-11-21 11:10:44] Start standard processing workflow...
#> ℹ [2026-11-21 11:10:44] Checking a list of <Seurat>...
#> ! [2026-11-21 11:10:44] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-11-21 11:10:44] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-11-21 11:10:44] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-11-21 11:10:44] Use the separate HVF from `srt_list`
#> ℹ [2026-11-21 11:10:44] Number of available HVF: 2000
#> ℹ [2026-11-21 11:10:44] Finished check
#> ℹ [2026-11-21 11:10:44] Perform `ScaleData()`
#> ℹ [2026-11-21 11:10:44] Perform pca linear dimension reduction
#> ℹ [2026-11-21 11:10:45] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-11-21 11:10:45] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-11-21 11:10:45] Reorder clusters...
#> ℹ [2026-11-21 11:10:45] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-11-21 11:10:45] Perform umap nonlinear dimension reduction
#> ✔ [2026-11-21 11:10:54] Standard processing workflow completed
pancreas_sub <- RunGLMPCA(pancreas_sub)
#> ℹ GLMPC_ 1 
#> ℹ Positive:  Gm26633, Ctgf, Anxa1, Col6a1, Sparcl1, Sp140, Col23a1, Isg15, Ctsk, Zfp385b 
#> ℹ      Islr, Col1a1, Galnt16, Tnni3, Hoxb4, Cxcl10, Platr22, Plscr2, P2ry2, Col1a2 
#> ℹ      Smpx, Pmp22, Gm15640, Kcnj8, Olfml2a, Edn1, 1110002O04Rik, Ccl20, Cxcl16, Peg12 
#> ℹ Negative:  Barx2, Gm3448, Fam71b, Gm933, Ky, Ptger3, Mesp1, Cartpt, Gm6086, Cypt3 
#> ℹ      Ngf, Bcas1os1, 3930402G23Rik, Platr26, Dusp26, Srrm3, Prl, Nlgn1, Klk11, Bhlhe23 
#> ℹ      Gm4319, 1700001C02Rik, 2410021H03Rik, Foxd3, Ucn3, Doc2a, Lrrc6, Cdkn2b, Gad1, Il1r2 
#> ℹ GLMPC_ 2 
#> ℹ Positive:  Islr, Kcnj8, Sparcl1, Tmem119, Galnt16, Col1a1, Col3a1, L1td1, Col1a2, Fam71b 
#> ℹ      Col6a1, Klk11, Oasl2, Jakmip3, Pdcd1, Sapcd1, Col5a1, Gsg1l, Adra2c, Pkd2l1 
#> ℹ      Col23a1, Calb1, Wnt3, Tmem114, Acoxl, 1700128E19Rik, Pid1, Gcg, Gm6086, Tfap2c 
#> ℹ Negative:  Sst, Dkk2, RP23-428N8.3, Gm26633, Kcne2, Tnni3, Col25a1, Rac2, Ceacam10, Ctsk 
#> ℹ      Gtf2ird2, Klhl14, Ctgf, Gm6410, D7Ertd443e, Ppy, Afap1l2, Elovl4, RP23-58K20.3, Nov 
#> ℹ      Tstd1, Srgn, Platr27, Cxcl10, Cttnbp2, Fcgr3, Lst1, Spock1, 1700015F17Rik, Itgb7 
#> ℹ GLMPC_ 3 
#> ℹ Positive:  Ncf2, Rac2, Cd37, P2ry14, Fcgr3, Srgn, Ifitm1, Tyrobp, Prom2, Fn3krp 
#> ℹ      Tex36, Fgf8, Gm11636, Slfn2, Gm16140, Laptm5, Snai2, Itgb7, Pthlh, Bhlhe23 
#> ℹ      Tgm7, Cmklr1, Tox2, Fam212b, Gm8773, Siglece, S100a14, Mapkapk3, Sema3g, Gmfg 
#> ℹ Negative:  Sparcl1, Dkk2, 4930539E08Rik, Gcg, Fgb, Gm26633, Col6a1, Lgr5, Galnt16, Col1a2 
#> ℹ      RP23-428N8.3, Gad2, Calb1, Kng2, Ryr3, Ctgf, 4930426D05Rik, Col23a1, Tnni3, Pid1 
#> ℹ      Npy, Hoxb4, Otc, Klhl14, Zfp385b, Crygn, Barx2, Islr, Doc2a, 1700001C02Rik 
#> ℹ GLMPC_ 4 
#> ℹ Positive:  Cbln4, Islr, Sp140, Ghrl, Kcnj8, Irs4, Tex36, Anxa1, Foxd3, Col1a1 
#> ℹ      Col6a1, Lrrtm3, Arhgap22, Col1a2, Guca2a, Col23a1, Sst, Sparcl1, Ngf, Col3a1 
#> ℹ      Gm17455, Kcne2, Avp, L1td1, Oasl2, Ctsk, D7Ertd443e, Prickle2, Fam71b, St8sia2 
#> ℹ Negative:  Kng2, Lmx1a, Spns2, Pif1, 4930539E08Rik, Slc4a10, Ryr3, Noxred1, Il1r2, Cldn18 
#> ℹ      Lrrc6, Cbln1, Pid1, Sapcd2, Ceacam10, Iqgap3, Espl1, Depdc1a, Aurkb, Npy 
#> ℹ      Ckap2l, Mxd3, Ccnb1, Nusap1, Neil3, 1700024G13Rik, Esco2, P2ry1, Slc16a6, Hist1h1a 
#> ℹ GLMPC_ 5 
#> ℹ Positive:  Anxa1, Nepn, Cxcl10, Gast, Col3a1, Gm29440, Fgf8, Cxcl16, Krt17, Col1a2 
#> ℹ      Gm15640, Tagln, A930017K11Rik, Ltb, Arhgap22, Ccl20, Esyt3, RP23-4H17.3, Ifit1bl1, Nrn1 
#> ℹ      Sytl2, Cypt3, Rerg, Colec12, Cttnbp2, Gm38112, Sapcd1, Otc, Col5a1, Ngf 
#> ℹ Negative:  Cbln4, Dlgap1, Sst, Ins1, Sds, Gsg1l, Sparcl1, Npy, Gm6878, Adam32 
#> ℹ      Smpx, BC043934, Ins2, Olfml2a, Snai2, Fam71b, Pkd2l1, Sp5, Sftpd, Rac2 
#> ℹ      Gm11789, Oasl2, Syndig1l, Col25a1, Nhs, Tmem255b, Gm28875, P2ry2, Sostdc1, Atoh8 
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "glmpca"
)
```
