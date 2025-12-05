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
#> + VIM                 6.2.6  
#> + destiny             3.24.0 [bld][cmp]
#> + e1071               1.7-16 
#> + ggplot.multistats   1.0.1  
#> + ggthemes            5.2.0  
#> + knn.covertree       1.0    
#> + laeken              0.5.3  
#> + org.Mm.eg.db        3.22.0 [bld]
#> + pcaMethods          2.2.0  [bld][cmp]
#> + proxy               0.4-27 
#> + ranger              0.17.0 
#> + repr                1.1.7  
#> + rgl                 1.3.31  + ✔ libfreetype6-dev, ✖ libglu1-mesa-dev, ✖ texlive, ✔ libpng-dev, ✖ libgl1-mesa-dev, ✔ pandoc, ✔ zlib1g-dev
#> + robustbase          0.99-6 
#> + scatterplot3d       0.3-44 
#> + smoother            1.3    
#> + vcd                 1.4-13 
#> → Will install 3 system packages:
#> + libgl1-mesa-dev   - rgl
#> + libglu1-mesa-dev  - rgl
#> + texlive           - rgl
#> ℹ No downloads are needed, 18 pkgs are cached
#> ✔ Got e1071 1.7-16 (x86_64-pc-linux-gnu-ubuntu-24.04) (596.57 kB)
#> ✔ Got DEoptimR 1.1-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (74.64 kB)
#> ✔ Got smoother 1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (22.69 kB)
#> ✔ Got destiny 3.24.0 (source) (900.68 kB)
#> ✔ Got proxy 0.4-27 (x86_64-pc-linux-gnu-ubuntu-24.04) (175.47 kB)
#> ✔ Got pcaMethods 2.2.0 (source) (1.05 MB)
#> ✔ Got ranger 0.17.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (473.78 kB)
#> ✔ Got knn.covertree 1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (91.51 kB)
#> ✔ Got ggthemes 5.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (969.47 kB)
#> ✔ Got ggplot.multistats 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (31.64 kB)
#> ✔ Got laeken 0.5.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.16 MB)
#> ✔ Got scatterplot3d 0.3-44 (x86_64-pc-linux-gnu-ubuntu-24.04) (348.15 kB)
#> ✔ Got VIM 6.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.77 MB)
#> ✔ Got vcd 1.4-13 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.30 MB)
#> ✔ Got repr 1.1.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (128.02 kB)
#> ✔ Got robustbase 0.99-6 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.13 MB)
#> ✔ Got rgl 1.3.31 (x86_64-pc-linux-gnu-ubuntu-24.04) (6.45 MB)
#> ✔ Got org.Mm.eg.db 3.22.0 (source) (92.82 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install libfreetype6-dev libglu1-mesa-dev texlive libpng-dev libgl1-mesa-dev pandoc zlib1g-dev make libcurl4-openssl-dev libnode-dev libxml2-dev libx11-dev libssl-dev libglpk-dev cmake libjpeg-dev libtiff-dev libwebp-dev libicu-dev libfontconfig1-dev libfribidi-dev libharfbuzz-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> pandoc is already the newest version (3.1.3+ds-2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> make is already the newest version (4.3-4.1build2).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libnode-dev is already the newest version (18.19.1+dfsg-6ubuntu5).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libx11-dev is already the newest version (2:1.8.7-1build1).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> libglpk-dev is already the newest version (5.0-1build2).
#> cmake is already the newest version (3.28.3-1build7).
#> libjpeg-dev is already the newest version (8c-2ubuntu11).
#> libtiff-dev is already the newest version (4.5.1+git230720-4ubuntu2.4).
#> libwebp-dev is already the newest version (1.3.2-0.4build3).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> libfribidi-dev is already the newest version (1.0.13-3build1).
#> libharfbuzz-dev is already the newest version (8.3.0-2build2).
#> The following additional packages will be installed:
#>   dvisvgm fonts-lmodern fonts-texgyre fonts-texgyre-math libegl-dev
#>   libegl-mesa0 libegl1 libgl-dev libgles-dev libgles1 libgles2 libglu1-mesa
#>   libglvnd-core-dev libglvnd-dev libglx-dev libgumbo2 libkpathsea6 libmujs3
#> libopengl-dev libopengl0 libpotrace0 libptexenc1 libsynctex2 libteckit0
#>   libtexlua53-5 libwoff1 libzzip-0-13t64 lmodern mupdf-tools tex-gyre
#>   texlive-base texlive-binaries texlive-fonts-recommended texlive-latex-base
#>   texlive-latex-recommended tipa
#> Suggested packages:
#>   perl-tk xpdf | pdf-viewer xzdec texlive-binaries-sse2 hintview
#>   texlive-fonts-recommended-doc texlive-latex-base-doc wp2latex
#>   texlive-latex-recommended-doc texlive-luatex texlive-pstricks tipa-doc
#> The following NEW packages will be installed:
#>   dvisvgm fonts-lmodern fonts-texgyre fonts-texgyre-math libegl-dev
#>   libegl-mesa0 libegl1 libgl-dev libgl1-mesa-dev libgles-dev libgles1 libgles2
#>   libglu1-mesa libglu1-mesa-dev libglvnd-core-dev libglvnd-dev libglx-dev
#> libgumbo2 libkpathsea6 libmujs3 libopengl-dev libopengl0 libpotrace0
#>   libptexenc1 libsynctex2 libteckit0 libtexlua53-5 libwoff1 libzzip-0-13t64
#>   lmodern mupdf-tools tex-gyre texlive texlive-base texlive-binaries
#>   texlive-fonts-recommended texlive-latex-base texlive-latex-recommended tipa
#> 0 upgraded, 39 newly installed, 0 to remove and 49 not upgraded.
#> Need to get 132 MB of archives.
#> After this operation, 327 MB of additional disk space will be used.
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Get:2 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libkpathsea6 amd64 2023.20230311.66589-9build3 [63.0 kB]
#> Get:3 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libpotrace0 amd64 1.16-2build1 [17.7 kB]
#> Get:4 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libwoff1 amd64 1.0.2-2build1 [45.3 kB]
#> Get:5 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 dvisvgm amd64 3.2.1+ds-1build1 [1042 kB]
#> Get:6 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 fonts-lmodern all 2.005-1 [4799 kB]
#> Get:7 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 fonts-texgyre all 20180621-6 [8350 kB]
#> Get:8 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 fonts-texgyre-math all 20180621-6 [2373 kB]
#> Get:9 http://azure.archive.ubuntu.com/ubuntu noble-updates/main amd64 libegl-mesa0 amd64 25.0.7-0ubuntu0.24.04.2 [124 kB]
#> Get:10 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libegl1 amd64 1.7.0-1build1 [28.7 kB]
#> Get:11 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libglx-dev amd64 1.7.0-1build1 [14.2 kB]
#> Get:12 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libgl-dev amd64 1.7.0-1build1 [102 kB]
#> Get:13 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libegl-dev amd64 1.7.0-1build1 [18.2 kB]
#> Get:14 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libgles1 amd64 1.7.0-1build1 [11.6 kB]
#> Get:15 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libgles2 amd64 1.7.0-1build1 [17.1 kB]
#> Get:16 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libgles-dev amd64 1.7.0-1build1 [50.5 kB]
#> Get:17 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libopengl0 amd64 1.7.0-1build1 [32.8 kB]
#> Get:18 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libglu1-mesa amd64 9.0.2-1.1build1 [152 kB]
#> Get:19 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libopengl-dev amd64 1.7.0-1build1 [3454 B]
#> Get:20 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libglu1-mesa-dev amd64 9.0.2-1.1build1 [237 kB]
#> Get:21 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libgumbo2 amd64 0.12.0+dfsg-2build1 [126 kB]
#> Get:22 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libmujs3 amd64 1.3.3-3build2 [134 kB]
#> Get:23 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libptexenc1 amd64 2023.20230311.66589-9build3 [40.4 kB]
#> Get:24 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libsynctex2 amd64 2023.20230311.66589-9build3 [59.6 kB]
#> Get:25 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libteckit0 amd64 2.5.12+ds1-1 [411 kB]
#> Get:26 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libtexlua53-5 amd64 2023.20230311.66589-9build3 [123 kB]
#> Get:27 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libzzip-0-13t64 amd64 0.13.72+dfsg.1-1.2build1 [28.1 kB]
#> Get:28 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 lmodern all 2.005-1 [9542 kB]
#> Get:29 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 mupdf-tools amd64 1.23.10+ds1-1build3 [49.3 MB]
#> Get:30 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 tex-gyre all 20180621-6 [6396 kB]
#> Get:31 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 texlive-binaries amd64 2023.20230311.66589-9build3 [8529 kB]
#> Get:32 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 texlive-base all 2023.20240207-1 [21.7 MB]
#> Get:33 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 texlive-fonts-recommended all 2023.20240207-1 [4973 kB]
#> Get:34 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 texlive-latex-base all 2023.20240207-1 [1238 kB]
#> Get:35 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 texlive-latex-recommended all 2023.20240207-1 [8826 kB]
#> Get:36 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 texlive all 2023.20240207-1 [14.0 kB]
#> Get:37 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 tipa all 2:1.3-21 [2967 kB]
#> Get:38 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libglvnd-core-dev amd64 1.7.0-1build1 [13.6 kB]
#> Get:39 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libglvnd-dev amd64 1.7.0-1build1 [3198 B]
#> Get:40 http://azure.archive.ubuntu.com/ubuntu noble-updates/main amd64 libgl1-mesa-dev amd64 25.0.7-0ubuntu0.24.04.2 [20.4 kB]
#> Preconfiguring packages ...
#> Fetched 132 MB in 1s (131 MB/s)
#> Selecting previously unselected package libkpathsea6:amd64.
#> (Reading database ...
#> (Reading database ... 5%(Reading database ... 10%(Reading database ... 15%(Reading database ... 20%(Reading database ... 25%(Reading database ... 30%(Reading database ... 35%(Reading database ... 40%(Reading database ... 45%(Reading database ... 50%(Reading database ... 55%
#> (Reading database ... 60%
#> (Reading database ... 65%
#> (Reading database ... 70%
#> (Reading database ... 75%
#> (Reading database ... 80%
#> (Reading database ... 85%
#> (Reading database ... 90%
#> (Reading database ... 95%
#> (Reading database ... 100%(Reading database ... 256281 files and directories currently installed.)
#> Preparing to unpack .../00-libkpathsea6_2023.20230311.66589-9build3_amd64.deb ...
#> Unpacking libkpathsea6:amd64 (2023.20230311.66589-9build3) ...
#> Selecting previously unselected package libpotrace0:amd64.
#> Preparing to unpack .../01-libpotrace0_1.16-2build1_amd64.deb ...
#> Unpacking libpotrace0:amd64 (1.16-2build1) ...
#> Selecting previously unselected package libwoff1:amd64.
#> Preparing to unpack .../02-libwoff1_1.0.2-2build1_amd64.deb ...
#> Unpacking libwoff1:amd64 (1.0.2-2build1) ...
#> Selecting previously unselected package dvisvgm.
#> Preparing to unpack .../03-dvisvgm_3.2.1+ds-1build1_amd64.deb ...
#> Unpacking dvisvgm (3.2.1+ds-1build1) ...
#> Selecting previously unselected package fonts-lmodern.
#> Preparing to unpack .../04-fonts-lmodern_2.005-1_all.deb ...
#> Unpacking fonts-lmodern (2.005-1) ...
#> Selecting previously unselected package fonts-texgyre.
#> Preparing to unpack .../05-fonts-texgyre_20180621-6_all.deb ...
#> Unpacking fonts-texgyre (20180621-6) ...
#> Selecting previously unselected package fonts-texgyre-math.
#> Preparing to unpack .../06-fonts-texgyre-math_20180621-6_all.deb ...
#> Unpacking fonts-texgyre-math (20180621-6) ...
#> Selecting previously unselected package libegl-mesa0:amd64.
#> Preparing to unpack .../07-libegl-mesa0_25.0.7-0ubuntu0.24.04.2_amd64.deb ...
#> Unpacking libegl-mesa0:amd64 (25.0.7-0ubuntu0.24.04.2) ...
#> Selecting previously unselected package libegl1:amd64.
#> Preparing to unpack .../08-libegl1_1.7.0-1build1_amd64.deb ...
#> Unpacking libegl1:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libglx-dev:amd64.
#> Preparing to unpack .../09-libglx-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libglx-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libgl-dev:amd64.
#> Preparing to unpack .../10-libgl-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libgl-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libegl-dev:amd64.
#> Preparing to unpack .../11-libegl-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libegl-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libgles1:amd64.
#> Preparing to unpack .../12-libgles1_1.7.0-1build1_amd64.deb ...
#> Unpacking libgles1:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libgles2:amd64.
#> Preparing to unpack .../13-libgles2_1.7.0-1build1_amd64.deb ...
#> Unpacking libgles2:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libgles-dev:amd64.
#> Preparing to unpack .../14-libgles-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libgles-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libopengl0:amd64.
#> Preparing to unpack .../15-libopengl0_1.7.0-1build1_amd64.deb ...
#> Unpacking libopengl0:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libglu1-mesa:amd64.
#> Preparing to unpack .../16-libglu1-mesa_9.0.2-1.1build1_amd64.deb ...
#> Unpacking libglu1-mesa:amd64 (9.0.2-1.1build1) ...
#> Selecting previously unselected package libopengl-dev:amd64.
#> Preparing to unpack .../17-libopengl-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libopengl-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libglu1-mesa-dev:amd64.
#> Preparing to unpack .../18-libglu1-mesa-dev_9.0.2-1.1build1_amd64.deb ...
#> Unpacking libglu1-mesa-dev:amd64 (9.0.2-1.1build1) ...
#> Selecting previously unselected package libgumbo2:amd64.
#> Preparing to unpack .../19-libgumbo2_0.12.0+dfsg-2build1_amd64.deb ...
#> Unpacking libgumbo2:amd64 (0.12.0+dfsg-2build1) ...
#> Selecting previously unselected package libmujs3:amd64.
#> Preparing to unpack .../20-libmujs3_1.3.3-3build2_amd64.deb ...
#> Unpacking libmujs3:amd64 (1.3.3-3build2) ...
#> Selecting previously unselected package libptexenc1:amd64.
#> Preparing to unpack .../21-libptexenc1_2023.20230311.66589-9build3_amd64.deb ...
#> Unpacking libptexenc1:amd64 (2023.20230311.66589-9build3) ...
#> Selecting previously unselected package libsynctex2:amd64.
#> Preparing to unpack .../22-libsynctex2_2023.20230311.66589-9build3_amd64.deb ...
#> Unpacking libsynctex2:amd64 (2023.20230311.66589-9build3) ...
#> Selecting previously unselected package libteckit0:amd64.
#> Preparing to unpack .../23-libteckit0_2.5.12+ds1-1_amd64.deb ...
#> Unpacking libteckit0:amd64 (2.5.12+ds1-1) ...
#> Selecting previously unselected package libtexlua53-5:amd64.
#> Preparing to unpack .../24-libtexlua53-5_2023.20230311.66589-9build3_amd64.deb ...
#> Unpacking libtexlua53-5:amd64 (2023.20230311.66589-9build3) ...
#> Selecting previously unselected package libzzip-0-13t64:amd64.
#> Preparing to unpack .../25-libzzip-0-13t64_0.13.72+dfsg.1-1.2build1_amd64.deb ...
#> Unpacking libzzip-0-13t64:amd64 (0.13.72+dfsg.1-1.2build1) ...
#> Selecting previously unselected package lmodern.
#> Preparing to unpack .../26-lmodern_2.005-1_all.deb ...
#> Unpacking lmodern (2.005-1) ...
#> Selecting previously unselected package mupdf-tools.
#> Preparing to unpack .../27-mupdf-tools_1.23.10+ds1-1build3_amd64.deb ...
#> Unpacking mupdf-tools (1.23.10+ds1-1build3) ...
#> Selecting previously unselected package tex-gyre.
#> Preparing to unpack .../28-tex-gyre_20180621-6_all.deb ...
#> Unpacking tex-gyre (20180621-6) ...
#> Selecting previously unselected package texlive-binaries.
#> Preparing to unpack .../29-texlive-binaries_2023.20230311.66589-9build3_amd64.deb ...
#> Unpacking texlive-binaries (2023.20230311.66589-9build3) ...
#> Selecting previously unselected package texlive-base.
#> Preparing to unpack .../30-texlive-base_2023.20240207-1_all.deb ...
#> Unpacking texlive-base (2023.20240207-1) ...
#> Selecting previously unselected package texlive-fonts-recommended.
#> Preparing to unpack .../31-texlive-fonts-recommended_2023.20240207-1_all.deb ...
#> Unpacking texlive-fonts-recommended (2023.20240207-1) ...
#> Selecting previously unselected package texlive-latex-base.
#> Preparing to unpack .../32-texlive-latex-base_2023.20240207-1_all.deb ...
#> Unpacking texlive-latex-base (2023.20240207-1) ...
#> Selecting previously unselected package texlive-latex-recommended.
#> Preparing to unpack .../33-texlive-latex-recommended_2023.20240207-1_all.deb ...
#> Unpacking texlive-latex-recommended (2023.20240207-1) ...
#> Selecting previously unselected package texlive.
#> Preparing to unpack .../34-texlive_2023.20240207-1_all.deb ...
#> Unpacking texlive (2023.20240207-1) ...
#> Selecting previously unselected package tipa.
#> Preparing to unpack .../35-tipa_2%3a1.3-21_all.deb ...
#> Unpacking tipa (2:1.3-21) ...
#> Selecting previously unselected package libglvnd-core-dev:amd64.
#> Preparing to unpack .../36-libglvnd-core-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libglvnd-core-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libglvnd-dev:amd64.
#> Preparing to unpack .../37-libglvnd-dev_1.7.0-1build1_amd64.deb ...
#> Unpacking libglvnd-dev:amd64 (1.7.0-1build1) ...
#> Selecting previously unselected package libgl1-mesa-dev:amd64.
#> Preparing to unpack .../38-libgl1-mesa-dev_25.0.7-0ubuntu0.24.04.2_amd64.deb ...
#> Unpacking libgl1-mesa-dev:amd64 (25.0.7-0ubuntu0.24.04.2) ...
#> Setting up libglvnd-core-dev:amd64 (1.7.0-1build1) ...
#> Setting up fonts-texgyre-math (20180621-6) ...
#> Setting up libwoff1:amd64 (1.0.2-2build1) ...
#> Setting up libmujs3:amd64 (1.3.3-3build2) ...
#> Setting up libopengl0:amd64 (1.7.0-1build1) ...
#> Setting up libegl-mesa0:amd64 (25.0.7-0ubuntu0.24.04.2) ...
#> Setting up libgles2:amd64 (1.7.0-1build1) ...
#> Setting up libzzip-0-13t64:amd64 (0.13.72+dfsg.1-1.2build1) ...
#> Setting up libgumbo2:amd64 (0.12.0+dfsg-2build1) ...
#> Setting up libteckit0:amd64 (2.5.12+ds1-1) ...
#> Setting up libgles1:amd64 (1.7.0-1build1) ...
#> Setting up libtexlua53-5:amd64 (2023.20230311.66589-9build3) ...
#> Setting up fonts-texgyre (20180621-6) ...
#> Setting up libkpathsea6:amd64 (2023.20230311.66589-9build3) ...
#> Setting up libegl1:amd64 (1.7.0-1build1) ...
#> Setting up fonts-lmodern (2.005-1) ...
#> Setting up libglx-dev:amd64 (1.7.0-1build1) ...
#> Setting up libglu1-mesa:amd64 (9.0.2-1.1build1) ...
#> Setting up libopengl-dev:amd64 (1.7.0-1build1) ...
#> Setting up tex-gyre (20180621-6) ...
#> Setting up libgl-dev:amd64 (1.7.0-1build1) ...
#> Setting up libsynctex2:amd64 (2023.20230311.66589-9build3) ...
#> Setting up libpotrace0:amd64 (1.16-2build1) ...
#> Setting up libegl-dev:amd64 (1.7.0-1build1) ...
#> Setting up dvisvgm (3.2.1+ds-1build1) ...
#> Setting up mupdf-tools (1.23.10+ds1-1build3) ...
#> Setting up libptexenc1:amd64 (2023.20230311.66589-9build3) ...
#> Setting up libglu1-mesa-dev:amd64 (9.0.2-1.1build1) ...
#> Setting up texlive-binaries (2023.20230311.66589-9build3) ...
#> update-alternatives: using /usr/bin/xdvi-xaw to provide /usr/bin/xdvi.bin (xdvi.bin) in auto mode
#> update-alternatives: using /usr/bin/bibtex.original to provide /usr/bin/bibtex (bibtex) in auto mode
#> Setting up lmodern (2.005-1) ...
#> Setting up libgles-dev:amd64 (1.7.0-1build1) ...
#> Setting up texlive-base (2023.20240207-1) ...
#> tl-paper: setting paper size for dvips to a4: /var/lib/texmf/dvips/config/config-paper.ps
#> tl-paper: setting paper size for dvipdfmx to a4: /var/lib/texmf/dvipdfmx/dvipdfmx-paper.cfg
#> tl-paper: setting paper size for xdvi to a4: /var/lib/texmf/xdvi/XDvi-paper
#> tl-paper: setting paper size for pdftex to a4: /var/lib/texmf/tex/generic/tex-ini-files/pdftexconfig.tex
#> Setting up libglvnd-dev:amd64 (1.7.0-1build1) ...
#> Setting up texlive-latex-base (2023.20240207-1) ...
#> Setting up texlive-latex-recommended (2023.20240207-1) ...
#> Setting up libgl1-mesa-dev:amd64 (25.0.7-0ubuntu0.24.04.2) ...
#> Setting up texlive-fonts-recommended (2023.20240207-1) ...
#> Setting up tipa (2:1.3-21) ...
#> Setting up texlive (2023.20240207-1) ...
#> Processing triggers for man-db (2.12.0-4build2) ...
#> Not building database; man-db/auto-update is not 'true'.
#> Processing triggers for tex-common (6.18) ...
#> Running mktexlsr. This may take some time...
#> done.
#> Running updmap-sys. This may take some time...
#> done.
#> Running mktexlsr /var/lib/texmf ...
#> done.
#> Building format(s) --all.
#>  This may take some time...
#> done.
#> Processing triggers for install-info (7.1-3build2) ...
#> Processing triggers for fontconfig (2.15.0-1.1ubuntu2) ...
#> Processing triggers for libc-bin (2.39-0ubuntu8.6) ...
#> Running kernel seems to be up-to-date.
#> 
#> Restarting services...
#> Service restarts being deferred:
#>  systemctl restart networkd-dispatcher.service
#> 
#> No containers need to be restarted.
#> 
#> No user sessions are running outdated binaries.
#> 
#> No VM guests are running outdated hypervisor (qemu) binaries on this host.
#> ℹ Building pcaMethods 2.2.0
#> ℹ Building org.Mm.eg.db 3.22.0
#> ✔ Installed DEoptimR 1.1-4  (100ms)
#> ✔ Installed e1071 1.7-16  (89ms)
#> ✔ Installed ggplot.multistats 1.0.1  (93ms)
#> ✔ Installed ggthemes 5.2.0  (88ms)
#> ✔ Installed knn.covertree 1.0  (148ms)
#> ✔ Installed laeken 0.5.3  (93ms)
#> ✔ Installed proxy 0.4-27  (83ms)
#> ✔ Installed ranger 0.17.0  (84ms)
#> ✔ Installed repr 1.1.7  (89ms)
#> ✔ Installed robustbase 0.99-6  (109ms)
#> ✔ Installed rgl 1.3.31  (279ms)
#> ✔ Installed scatterplot3d 0.3-44  (83ms)
#> ✔ Installed smoother 1.3  (81ms)
#> ✔ Installed vcd 1.4-13  (88ms)
#> ✔ Installed VIM 6.2.6  (73ms)
#> ✔ Built pcaMethods 2.2.0 (9.6s)
#> ✔ Installed pcaMethods 2.2.0  (41ms)
#> ℹ Building destiny 3.24.0
#> ✔ Built destiny 3.24.0 (28.6s)
#> ✔ Installed destiny 3.24.0  (60ms)
#> ✔ Built org.Mm.eg.db 3.22.0 (3m 31.8s)
#> ✔ Installed org.Mm.eg.db 3.22.0  (1.9s)
#> ✔ 1 pkg + 256 deps: kept 239, added 18, dld 18 (114.49 MB) [4m 30.7s]
#> Registered S3 methods overwritten by 'proxy':
#>   method               from    
#>   print.registry_field registry
#>   print.registry_entry registry
#> 'as(<dsCMatrix>, "dgTMatrix")' is deprecated.
#> Use 'as(as(., "generalMatrix"), "TsparseMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").

CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "dm"
)
```
