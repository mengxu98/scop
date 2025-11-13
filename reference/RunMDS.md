# Run MDS (multi-dimensional scaling)

Run MDS (multi-dimensional scaling)

## Usage

``` r
RunMDS(object, ...)

# S3 method for class 'Seurat'
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.name = "mds",
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay'
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# S3 method for class 'Assay5'
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.key = "MDS_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunMDS(
  object,
  assay = NULL,
  layer = "data",
  nmds = 50,
  dist.method = "euclidean",
  mds.method = "cmdscale",
  rev.mds = FALSE,
  reduction.key = "MDS_",
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

  Additional arguments to be passed to
  [stats::cmdscale](https://rdrr.io/r/stats/cmdscale.html),
  [MASS::isoMDS](https://rdrr.io/pkg/MASS/man/isoMDS.html) or
  [MASS::sammon](https://rdrr.io/pkg/MASS/man/sammon.html).

- assay:

  The assay to be used for the analysis. Default is `NULL`.

- layer:

  The layer to be used for the analysis. Default is `"data"`.

- features:

  The features to be used for the analysis. Default is `NULL`, which
  uses all variable features.

- nmds:

  The number of dimensions to be computed. Default is `50`.

- dist.method:

  The distance metric to be used. Currently supported values are
  `"euclidean"`, `"chisquared"`, `"kullback"`, `"jeffreys"`, `"jensen"`,
  `"manhattan"`, `"maximum"`, `"canberra"`, `"minkowski"`, and
  `"hamming"`. Default is `"euclidean"`.

- mds.method:

  The MDS algorithm to be used. Currently supported values are
  `"cmdscale"`, `"isoMDS"`, and `"sammon"`. Default is `"cmdscale"`.

- rev.mds:

  Whether to perform reverse MDS (i.e., transpose the input matrix)
  before running the analysis. Default is `FALSE`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"mds"`.

- reduction.key:

  The prefix for the column names of the basis vectors. Default is
  `"MDS_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  The random seed to be used. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 12:28:37] Start standard scop workflow...
#> ℹ [2025-11-13 12:28:38] Checking a list of <Seurat> object...
#> ! [2025-11-13 12:28:38] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 12:28:38] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:28:40] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:28:40] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:28:41] Number of available HVF: 2000
#> ℹ [2025-11-13 12:28:41] Finished check
#> ℹ [2025-11-13 12:28:41] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:28:41] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 12:28:42] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:28:43] Reorder clusters...
#> ℹ [2025-11-13 12:28:43] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:28:43] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:28:43] UMAP will return its model
#> ℹ [2025-11-13 12:28:47] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:28:47] UMAP will return its model
#> ✔ [2025-11-13 12:28:52] Run scop standard workflow done
pancreas_sub <- RunMDS(pancreas_sub)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "mds"
)
```
