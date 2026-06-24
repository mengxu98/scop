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
#> ℹ [2026-06-24 04:21:03] Start standard processing workflow...
#> ℹ [2026-06-24 04:21:03] Checking a list of <Seurat>...
#> ! [2026-06-24 04:21:03] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:21:03] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:21:03] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:21:04] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 04:21:04] Number of available HVF: 2000
#> ℹ [2026-06-24 04:21:04] Finished check
#> ℹ [2026-06-24 04:21:04] Perform `ScaleData()`
#> ℹ [2026-06-24 04:21:04] Perform pca linear dimension reduction
#> ℹ [2026-06-24 04:21:05] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-24 04:21:05] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-24 04:21:05] Reorder clusters...
#> ℹ [2026-06-24 04:21:05] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-24 04:21:05] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-24 04:21:12] Standard processing workflow completed
pancreas_sub <- RunGLMPCA(pancreas_sub)
#> ℹ GLMPC_ 1 
#> ℹ Positive:  Sparcl1, Anxa1, Col6a1, Galnt16, Col23a1, Col1a1, Islr, Col1a2, Sp140, Isg15 
#> ℹ      Ctgf, Zfp385b, Gm26633, Tnni3, Ctsk, Kcnj8, Col3a1, Gdf15, P2ry2, Hoxb4 
#> ℹ      Platr22, Prickle2, Plscr2, Pmp22, 1110002O04Rik, Gm15640, Smpx, D630011A20Rik, Krt20, Pgr 
#> ℹ Negative:  Barx2, Gm3448, Ky, Ptger3, Platr26, Srrm3, Mesp1, Gm6086, Gm933, Bcas1os1 
#> ℹ      Dusp26, Gm4319, 3930402G23Rik, Fam71b, Lrrc6, Cypt3, Cartpt, Nlgn1, Cdkn2b, Gad1 
#> ℹ      Il1r2, A130057D12Rik, Ngf, Aard, Ucn3, 1700001C02Rik, 4930426D05Rik, Avp, Doc2a, 2410021H03Rik 
#> ℹ GLMPC_ 2 
#> ℹ Positive:  Dkk2, Fgb, Sparcl1, Crygn, RP23-428N8.3, Klhl14, Otc, 4930539E08Rik, Lgr5, Gcg 
#> ℹ      Gdf15, Gad2, Galnt16, Doc2a, 4930426D05Rik, Gm26633, Tnni3, Col1a2, Aif1, 1700001C02Rik 
#> ℹ      Ctgf, Col6a1, Tac1, Pif1, Rgs4, Col23a1, Pou6f2, Col25a1, Pyy, Aspm 
#> ℹ Negative:  Ncf2, Cd37, P2ry14, Ifitm1, Slfn2, Prom2, Fcgr3, Laptm5, Fn3krp, Tyrobp 
#> ℹ      Bhlhe22, Srgn, Tex36, Snai2, Tgm7, Gm16140, Cmklr1, Gm15567, Lmx1a, Siglece 
#> ℹ      Rac2, Neurog3, Pthlh, Coro1a, Tox2, Krt15, Sema3g, Gmfg, Gm13373, Rasgrp3 
#> ℹ GLMPC_ 3 
#> ℹ Positive:  Sst, Ceacam10, RP23-58K20.3, Gm6410, Srgn, Gtf2ird2, RP23-428N8.3, Afap1l2, Tmem100, Tnni3 
#> ℹ      Col25a1, Kcne2, Platr27, Rac2, Fcgr3, Cbln1, Gm10382, Dkk2, Ctgf, Spock1 
#> ℹ      Pik3c2b, Rasl11a, Tyrobp, Lst1, Gm26633, Coro1a, Tstd1, Itgb7, Igfbp3, Slc4a10 
#> ℹ Negative:  Islr, Kcnj8, Col1a1, Sparcl1, Galnt16, Col6a1, Col1a2, Col3a1, Tmem119, Col23a1 
#> ℹ      Guca2a, Ghrl, Fam71b, Sapcd1, Col5a1, L1td1, Gcg, Irs4, Klk11, Calb1 
#> ℹ      Pkd2l1, Ngf, Adra2c, Bcas1os1, Tfap2c, Bhlhe23, Gsg1l, Oasl2, Anxa1, Lsp1 
#> ℹ GLMPC_ 4 
#> ℹ Positive:  Pid1, Kng2, Nusap1, Sapcd2, Pif1, Neil3, Ccnb1, Aurkb, Mxd3, Iqgap3 
#> ℹ      Ckap2l, Prc1, Ska1, Plk1, Depdc1a, Aurka, Espl1, Hmmr, Esco2, Ube2c 
#> ℹ      Spns2, Tpx2, Cenpf, Nuf2, Ccnf, Sgol1, Kif20a, Troap, Kif11, Sgol2a 
#> ℹ Negative:  Rac2, Sst, Kcne2, Cbln4, Irs4, Ghrl, Sp140, Tex36, Lrrtm3, Anxa1 
#> ℹ      Itgb7, Foxd3, Elovl4, Gm933, Islr, Gpr6, Col1a1, Lst1, Ngf, Arhgap22 
#> ℹ      Gm17455, Col6a1, D7Ertd443e, Avp, Tstd1, Smpx, Ctsk, Col23a1, Bmp2, Nlgn1 
#> ℹ GLMPC_ 5 
#> ℹ Positive:  Gast, Anxa1, Fgf8, Nepn, Cxcl10, Gm15640, Cxcl16, Gm29440, Tagln, A930017K11Rik 
#> ℹ      Sytl2, Rerg, Col3a1, Esyt3, Krt17, Col1a2, Zfp97, Isg15, Cttnbp2, Ifit1bl1 
#> ℹ      Calb1, Nrn1, Pou6f2, Ccl20, Gm38112, Slfn9, Gcg, Cypt3, Pole, Ngf 
#> ℹ Negative:  Cbln4, Dlgap1, Ins1, Sparcl1, Sds, Gm6878, Aif1, Guca2a, Npy, Sst 
#> ℹ      Gsg1l, Aspm, Adam32, Rac2, Oasl2, Smpx, Fam71b, BC043934, Gm28875, Gm11789 
#> ℹ      Ins2, Crygn, Olfml2a, Nhs, Guca2b, Sftpd, Pmp22, Sp5, Pkd2l1, Sostdc1 
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "glmpca"
)
```
