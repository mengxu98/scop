# Run DM (diffusion map)

Run DM (diffusion map)

## Usage

``` r
RunDM(object, ...)

# S3 method for class 'Seurat'
RunDM(
  object,
  reduction = "pca",
  dims = 1:30,
  features = NULL,
  assay = NULL,
  layer = "data",
  ndcs = 2,
  sigma = "local",
  k = 30,
  dist.method = "euclidean",
  reduction.name = "dm",
  reduction.key = "DM_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunDM(
  object,
  assay = NULL,
  layer = "data",
  ndcs = 2,
  sigma = "local",
  k = 30,
  dist.method = "euclidean",
  reduction.key = "DM_",
  verbose = TRUE,
  seed.use = 11,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object or a matrix-like object.

- ...:

  Additional arguments to be passed to
  [destiny::DiffusionMap](https://rdrr.io/pkg/destiny/man/DiffusionMap-class.html).

- reduction:

  The reduction to be used. Default is `"pca"`.

- dims:

  The dimensions to be used. Default is `1:30`.

- features:

  The features to be used. Default is `NULL`.

- assay:

  The assay to be used. Default is `NULL`.

- layer:

  The layer to be used. Default is `"data"`.

- ndcs:

  A number of diffusion components (dimensions) to be computed. Default
  is `2`.

- sigma:

  The diffusion scale parameter of the Gaussian kernel. Currently
  supported values are `"local"` (default) and `"global"`.

- k:

  A number of nearest neighbors to be used for the construction of the
  graph. Default is `30`.

- dist.method:

  The distance metric to be used for the construction of the knn graph.
  Currently supported values are `"euclidean"` and `"cosine"`. Default
  is `"euclidean"`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"dm"`.

- reduction.key:

  The prefix for the column names of the basis vectors. Default is
  `"DM_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  An integer specifying the random seed to be used. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
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
pancreas_sub <- RunDM(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#>  
#> → Will install 18 packages.
#> → All 18 packages (0 B) are cached.
#> + DEoptimR            1.1-4  
#> + RcppHNSW            0.6.0  
#> + TTR                 0.24.4 
#> + VIM                 6.2.6  
#> + destiny             3.24.0 [bld][cmp]
#> + ggplot.multistats   1.0.1  
#> + ggthemes            5.2.0  
#> + knn.covertree       1.0    
#> + laeken              0.5.3  
#> + lmtest              0.9-40 
#> + pcaMethods          2.2.0  [bld][cmp]
#> + ranger              0.17.0 
#> + robustbase          0.99-6 
#> + smoother            1.3    
#> + sp                  2.2-0  
#> + vcd                 1.4-13 
#> + xts                 0.14.1 
#> + zoo                 1.8-14 
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 18 pkgs are cached
#> ✔ Got DEoptimR 1.1-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (74.64 kB)
#> ✔ Got RcppHNSW 0.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (217.19 kB)
#> ✔ Got lmtest 0.9-40 (x86_64-pc-linux-gnu-ubuntu-24.04) (403.49 kB)
#> ✔ Got pcaMethods 2.2.0 (source) (1.05 MB)
#> ✔ Got TTR 0.24.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (524.49 kB)
#> ✔ Got knn.covertree 1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (91.51 kB)
#> ✔ Got destiny 3.24.0 (source) (900.68 kB)
#> ✔ Got ranger 0.17.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (473.78 kB)
#> ✔ Got ggplot.multistats 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (31.64 kB)
#> ✔ Got xts 0.14.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.22 MB)
#> ✔ Got vcd 1.4-13 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.30 MB)
#> ✔ Got smoother 1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (22.69 kB)
#> ✔ Got VIM 6.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.77 MB)
#> ✔ Got zoo 1.8-14 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.03 MB)
#> ✔ Got laeken 0.5.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.16 MB)
#> ✔ Got robustbase 0.99-6 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.13 MB)
#> ✔ Got sp 2.2-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.31 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libcurl4-openssl-dev libssl-dev make cmake libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> cmake is already the newest version (3.28.3-1build7).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 62 not upgraded.
#> ℹ Building pcaMethods 2.2.0
#> ✔ Installed DEoptimR 1.1-4  (66ms)
#> ✔ Installed ggplot.multistats 1.0.1  (90ms)
#> ✔ Installed ggthemes 5.2.0  (165ms)
#> ✔ Installed knn.covertree 1.0  (68ms)
#> ✔ Installed laeken 0.5.3  (69ms)
#> ✔ Installed lmtest 0.9-40  (66ms)
#> ✔ Installed ranger 0.17.0  (68ms)
#> ✔ Installed RcppHNSW 0.6.0  (66ms)
#> ✔ Installed robustbase 0.99-6  (152ms)
#> ✔ Installed smoother 1.3  (151ms)
#> ✔ Installed TTR 0.24.4  (27ms)
#> ✔ Installed sp 2.2-0  (144ms)
#> ✔ Installed vcd 1.4-13  (70ms)
#> ✔ Installed VIM 6.2.6  (75ms)
#> ✔ Installed xts 0.14.1  (71ms)
#> ✔ Installed zoo 1.8-14  (49ms)
#> ✔ Built pcaMethods 2.2.0 (8.5s)
#> ✔ Installed pcaMethods 2.2.0  (1s)
#> ℹ Building destiny 3.24.0
#> ✔ Built destiny 3.24.0 (27.4s)
#> ✔ Installed destiny 3.24.0  (1.1s)
#> ✔ 1 pkg + 103 deps: kept 84, added 18, dld 17 (21.71 MB) [45s]
#> 'as(<dsCMatrix>, "dgTMatrix")' is deprecated.
#> Use 'as(as(., "generalMatrix"), "TsparseMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").

CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "dm"
)
```
