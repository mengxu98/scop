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

  Passed to the
  [glmpca::glmpca](https://rdrr.io/pkg/glmpca/man/glmpca.html) function.

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
#> ℹ [2026-09-27 22:10:19] Start standard processing workflow...
#> ℹ [2026-09-27 22:10:19] Checking a list of <Seurat>...
#> ! [2026-09-27 22:10:19] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:10:19] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 22:10:19] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 22:10:19] Use the separate HVF from `srt_list`
#> ℹ [2026-09-27 22:10:19] Number of available HVF: 2000
#> ℹ [2026-09-27 22:10:19] Finished check
#> ℹ [2026-09-27 22:10:19] Perform `ScaleData()`
#> ℹ [2026-09-27 22:10:19] Perform pca linear dimension reduction
#> ℹ [2026-09-27 22:10:19] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-27 22:10:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-27 22:10:20] Reorder clusters...
#> ℹ [2026-09-27 22:10:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-27 22:10:20] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-27 22:10:27] Standard processing workflow completed
pancreas_sub <- RunGLMPCA(pancreas_sub)
#> ℹ GLMPC_ 1 
#> ℹ Positive:  Sparcl1, Col6a1, Anxa1, Col23a1, Galnt16, Col1a1, Sp140, Islr, Ctgf, Isg15 
#> ℹ      Col1a2, Gm26633, Ctsk, Tnni3, Zfp385b, Hoxb4, P2ry2, Kcnj8, Plscr2, Pmp22 
#> ℹ      Platr22, Col3a1, Prickle2, Smpx, Cxcl10, 1110002O04Rik, Edn1, Gdf15, Gm6878, Apcs 
#> ℹ Negative:  Barx2, Gm3448, Ky, Mesp1, Ptger3, Platr26, Gm6086, Bcas1os1, Fam71b, Gm933 
#> ℹ      Cartpt, Srrm3, Cypt3, Dusp26, 3930402G23Rik, Gm4319, Ngf, Nlgn1, Lrrc6, Doc2a 
#> ℹ      1700001C02Rik, Cdkn2b, Gad1, Ucn3, Prl, Klk11, 2410021H03Rik, Acoxl, Il1r2, Aard 
#> ℹ GLMPC_ 2 
#> ℹ Positive:  Dkk2, Fgb, Sparcl1, RP23-428N8.3, Gm26633, Tnni3, Lgr5, Klhl14, 4930539E08Rik, Crygn 
#> ℹ      Otc, Ctgf, Gdf15, Gcg, 4930426D05Rik, Galnt16, Sst, Col25a1, Col6a1, Gad2 
#> ℹ      Doc2a, Barx2, Zfp385b, Col1a2, Pyy, Col23a1, Ppy, 1700001C02Rik, Calb1, Tac1 
#> ℹ Negative:  Ncf2, P2ry14, Cd37, Ifitm1, Prom2, Tex36, Slfn2, Tgm7, Bhlhe23, Fn3krp 
#> ℹ      Gm16140, Krt15, Cmklr1, Snai2, Tyrobp, Fgf8, Laptm5, Siglece, Bhlhe22, Tox2 
#> ℹ      Fam212b, Gm15567, Sema3g, Fcgr3, Pthlh, Gm8773, Neurog3, Neurod2, Epb42, S100a14 
#> ℹ GLMPC_ 3 
#> ℹ Positive:  Sst, Ceacam10, Gtf2ird2, Srgn, Gm6410, Rac2, Kcne2, RP23-58K20.3, RP23-428N8.3, Col25a1 
#> ℹ      Afap1l2, Fcgr3, Platr27, Rasl11a, Gm26633, Itgb7, Tmem100, Lst1, Spock1, Elovl4 
#> ℹ      Tstd1, Coro1a, D7Ertd443e, Tnni3, Klhl14, Igfbp3, Tyrobp, Lypd6b, Pik3c2b, Slc4a10 
#> ℹ Negative:  Islr, Kcnj8, Sparcl1, Col1a1, Col1a2, Galnt16, Col6a1, Col3a1, Tmem119, Col23a1 
#> ℹ      Fam71b, Guca2a, Ghrl, L1td1, Gcg, Pkd2l1, Col5a1, Sapcd1, Oasl2, Calb1 
#> ℹ      Adra2c, Bcas1os1, Irs4, Ngf, Tfap2c, Gsg1l, Tmem255b, Foxd3, Gm6878, Klk11 
#> ℹ GLMPC_ 4 
#> ℹ Positive:  Kng2, Pid1, Pif1, Nusap1, Sapcd2, Iqgap3, Ryr3, Mxd3, Ccnb1, Espl1 
#> ℹ      Spns2, Jakmip3, Cldn18, Neil3, Esco2, Casc5, Hmmr, Ckap2l, Aurkb, Plk1 
#> ℹ      Depdc1a, Npy, Noxred1, Cenpf, Parpbp, Prc1, Cdc25c, Troap, Kif2c, Hist1h1b 
#> ℹ Negative:  Rac2, Sst, Cbln4, Ghrl, Sp140, Kcne2, Anxa1, Irs4, Tex36, Foxd3 
#> ℹ      Lrrtm3, Islr, Arhgap22, D7Ertd443e, Gm17455, Ngf, Col1a1, Kcnj8, Ctsk, Elovl4 
#> ℹ      Itgb7, Gpr6, Col6a1, Col23a1, Lst1, Nlgn1, Tstd1, Mesp1, Gm933, Avp 
#> ℹ GLMPC_ 5 
#> ℹ Positive:  Anxa1, Gast, Nepn, Cxcl10, Fgf8, Col1a2, Col3a1, Gm29440, Cxcl16, Esyt3 
#> ℹ      Krt17, A930017K11Rik, Rerg, Cttnbp2, Sytl2, Ifit1bl1, Tagln, Gm15640, Ccl20, Cypt3 
#> ℹ      Pole, RP23-4H17.3, Ryr3, Ltb, Calb1, Pou6f2, Zfp97, Colec12, Nrn1, Sapcd1 
#> ℹ Negative:  Cbln4, Dlgap1, Ins1, Sst, Gsg1l, Sds, Gm6878, Adam32, Npy, Sparcl1 
#> ℹ      Aif1, Gm11789, Oasl2, Guca2a, Rac2, Olfml2a, Aspm, Fam71b, BC043934, Ins2 
#> ℹ      Sftpd, Snai2, Gm28875, Sp5, Syndig1l, Smpx, P2ry2, Guca2b, Rsad2, Tmem255b 
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "glmpca"
)
```
