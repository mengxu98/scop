# Run CellChat

RunCellChat performs CellChat analysis on a Seurat object to investigate
cell-to-cell communication. The results are stored in the Seurat object
and can be visualized using
[CellChatPlot](https://mengxu98.github.io/scop/reference/CellChatPlot.md).

## Usage

``` r
RunCellChat(
  srt,
  group.by,
  species = c("human", "mouse", "zebrafish"),
  split.by = NULL,
  annotation_selected = NULL,
  group_column = NULL,
  group_cmp = NULL,
  thresh = 0.05,
  min.cells = 10,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of the metadata column in the Seurat object that contains cell
  annotations.

- species:

  The species of the data, either 'human', 'mouse' or 'zebrafish'.

- split.by:

  Name of the metadata column in the Seurat object that contains sample
  information.

- annotation_selected:

  A vector of cell annotations of interest for running the CellChat
  analysis. If not provided, all cell types will be considered.

- group_column:

  Name of the metadata column in the Seurat object that defines
  conditions or groups.

- group_cmp:

  A list of pairwise condition comparisons for differential CellChat
  analysis.

- thresh:

  The threshold for computing centrality scores.

- min.cells:

  the minmum number of expressed cells required for the genes that are
  considered for cell-cell communication analysis.

- verbose:

  Whether to print the message. Default is `TRUE`.

## References

[CellChat](https://github.com/jinworks/CellChat),
[scDown::run_cellchatV2](https://htmlpreview.github.io/?https://raw.githubusercontent.com/BCH-RC/scDown/main/vignettes/scDown_CellChatV2.html)

## See also

[CellChatPlot](https://mengxu98.github.io/scop/reference/CellChatPlot.md)

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
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "mouse"
)
#>  
#> → Will install 62 packages.
#> → Will download 2 CRAN packages (1.72 MB), cached: 60 (0 B).
#> + BiocNeighbors    2.4.0    [bld][cmp]
#> + CellChat         2.2.0    [bld][cmp] (GitHub: 623f48f)
#> + ComplexHeatmap   2.26.0   [bld]
#> + DBI              1.2.3    
#> + GetoptLong       1.1.0     + ✔ perl
#> + GlobalOptions    0.1.3    
#> + NMF              0.28     [bld][cmp][dl] (1.68 MB)
#> + assorthead       1.4.0    [bld]
#> + bit              4.6.0    
#> + bit64            4.6.0-1  
#> + blob             1.2.4    
#> + brew             1.0-10   
#> + bsicons          0.1.2    
#> + cellranger       1.1.0    
#> + circlize         0.4.17   
#> + clipr            0.8.0     + ✔ libx11-dev
#> + clue             0.3-66   
#> + coda             0.19-4.1 
#> + conflicted       1.2.0    
#> + dbplyr           2.5.1    
#> + doParallel       1.0.17   
#> + dtplyr           1.3.2    
#> + foreach          1.5.2    
#> + gargle           1.6.0    
#> + gg.gap           1.3      
#> + ggalluvial       0.12.5   
#> + ggnetwork        0.5.14   
#> + googledrive      2.1.2    
#> + googlesheets4    1.1.2    
#> + gridBase         0.4-7    
#> + haven            2.5.5     + ✔ make, ✔ zlib1g-dev
#> + hms              1.1.4    
#> + ids              1.0.1    
#> + iterators        1.0.14   
#> + lubridate        1.9.4    
#> + network          1.19.0   
#> + ragg             1.5.0     + ✔ libfreetype6-dev, ✔ libjpeg-dev, ✔ libpng-dev, ✔ libtiff-dev, ✔ libwebp-dev
#> + readr            2.1.6    
#> + readxl           1.4.5    
#> + registry         0.5-1    
#> + rematch          2.0.0    
#> + rematch2         2.1.2    
#> + reprex           2.1.1     + ✔ pandoc
#> + rjson            0.2.23   
#> + rngtools         1.5.2    
#> + roxygen2         7.3.3    
#> + rstudioapi       0.17.1   
#> + rvest            1.0.5    
#> + selectr          0.5-0    
#> + shape            1.4.6.1  
#> + sna              2.8      
#> + statnet.common   4.12.0   
#> + svglite          2.2.2     + ✔ libpng-dev
#> + systemfonts      1.3.1     + ✔ libfontconfig1-dev, ✔ libfreetype6-dev
#> + textshaping      1.0.4     + ✔ libfreetype6-dev, ✔ libfribidi-dev, ✔ libharfbuzz-dev
#> + tidyverse        2.0.0    
#> + timechange       0.3.0    
#> + tzdb             0.5.0    
#> + uuid             1.2-1    
#> + vroom            1.6.7    
#> + wordcloud        2.6      [bld][cmp][dl] (42.45 kB)
#> + xml2             1.5.1     + ✔ libxml2-dev
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 2 pkgs (1.72 MB), 60 cached
#> ✔ Got BiocNeighbors 2.4.0 (source) (338.98 kB)
#> ✔ Got blob 1.2.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (47.47 kB)
#> ✔ Got bsicons 0.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (253.96 kB)
#> ✔ Got brew 1.0-10 (x86_64-pc-linux-gnu-ubuntu-24.04) (76.35 kB)
#> ✔ Got bit64 4.6.0-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (492.54 kB)
#> ✔ Got cellranger 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (103.66 kB)
#> ✔ Got bit 4.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (628.10 kB)
#> ✔ Got rematch2 2.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (46.02 kB)
#> ✔ Got conflicted 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (55.04 kB)
#> ✔ Got ComplexHeatmap 2.26.0 (source) (1.47 MB)
#> ✔ Got readxl 1.4.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (412.32 kB)
#> ✔ Got GlobalOptions 0.1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (483.55 kB)
#> ✔ Got rstudioapi 0.17.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (317.69 kB)
#> ✔ Got uuid 1.2-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (48.94 kB)
#> ✔ Got assorthead 1.4.0 (source) (1.91 MB)
#> ✔ Got statnet.common 4.12.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (313.83 kB)
#> ✔ Got iterators 1.0.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (346.57 kB)
#> ✔ Got clipr 0.8.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (51.24 kB)
#> ✔ Got foreach 1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (139.65 kB)
#> ✔ Got dtplyr 1.3.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (409.84 kB)
#> ✔ Got gridBase 0.4-7 (x86_64-pc-linux-gnu-ubuntu-24.04) (162.00 kB)
#> ✔ Got hms 1.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (103.38 kB)
#> ✔ Got wordcloud 2.6 (source) (41.56 kB)
#> ✔ Got vroom 1.6.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (966.05 kB)
#> ✔ Got googlesheets4 1.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (515.22 kB)
#> ✔ Got registry 0.5-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (194.47 kB)
#> ✔ Got selectr 0.5-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (570.89 kB)
#> ✔ Got NMF 0.28 (source) (1.68 MB)
#> ✔ Got lubridate 1.9.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (990.67 kB)
#> ✔ Got ggnetwork 0.5.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.97 MB)
#> ✔ Got timechange 0.3.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (174.23 kB)
#> ✔ Got rematch 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (16.20 kB)
#> ✔ Got network 1.19.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (823.12 kB)
#> ✔ Got textshaping 1.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (189.59 kB)
#> ✔ Got rjson 0.2.23 (x86_64-pc-linux-gnu-ubuntu-24.04) (113.12 kB)
#> ✔ Got reprex 2.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (495.27 kB)
#> ✔ Got gargle 1.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (618.93 kB)
#> ✔ Got tzdb 0.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (666.60 kB)
#> ✔ Got sna 2.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.28 MB)
#> ✔ Got dbplyr 2.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.24 MB)
#> ✔ Got readr 2.1.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (868.61 kB)
#> ✔ Got ragg 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (743.50 kB)
#> ✔ Got tidyverse 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (425.82 kB)
#> ✔ Got DBI 1.2.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (916.93 kB)
#> ✔ Got systemfonts 1.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (816.45 kB)
#> ✔ Got ggalluvial 0.12.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.66 MB)
#> ✔ Got googledrive 2.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.92 MB)
#> ✔ Got rvest 1.0.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (301.36 kB)
#> ✔ Got haven 2.5.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (393.30 kB)
#> ✔ Got rngtools 1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (78.03 kB)
#> ✔ Got doParallel 1.0.17 (x86_64-pc-linux-gnu-ubuntu-24.04) (188.11 kB)
#> ✔ Got clue 0.3-66 (x86_64-pc-linux-gnu-ubuntu-24.04) (998.14 kB)
#> ✔ Got shape 1.4.6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (749.97 kB)
#> ✔ Got gg.gap 1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (20.48 kB)
#> ✔ Got roxygen2 7.3.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (729.96 kB)
#> ✔ Got svglite 2.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (239.34 kB)
#> ✔ Got ids 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (122.46 kB)
#> ✔ Got coda 0.19-4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (327.75 kB)
#> ✔ Got GetoptLong 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.41 MB)
#> ✔ Got CellChat 2.2.0 (source) (29.82 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libx11-dev perl make zlib1g-dev libfreetype6-dev libjpeg-dev libpng-dev libtiff-dev libwebp-dev pandoc libfontconfig1-dev libfribidi-dev libharfbuzz-dev libxml2-dev libcurl4-openssl-dev libssl-dev libglpk-dev cmake python3 libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libx11-dev is already the newest version (2:1.8.7-1build1).
#> libx11-dev set to manually installed.
#> perl is already the newest version (5.38.2-3.2ubuntu0.2).
#> make is already the newest version (4.3-4.1build2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libjpeg-dev is already the newest version (8c-2ubuntu11).
#> libpng-dev is already the newest version (1.6.43-5ubuntu0.1).
#> libtiff-dev is already the newest version (4.5.1+git230720-4ubuntu2.4).
#> libwebp-dev is already the newest version (1.3.2-0.4build3).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> libfribidi-dev is already the newest version (1.0.13-3build1).
#> libharfbuzz-dev is already the newest version (8.3.0-2build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> libglpk-dev is already the newest version (5.0-1build2).
#> cmake is already the newest version (3.28.3-1build7).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 49 not upgraded.
#> ℹ Building wordcloud 2.6
#> ℹ Building assorthead 1.4.0
#> ✔ Installed bit 4.6.0  (114ms)
#> ✔ Installed bit64 4.6.0-1  (121ms)
#> ✔ Installed blob 1.2.4  (118ms)
#> ✔ Installed brew 1.0-10  (111ms)
#> ✔ Installed bsicons 0.1.2  (120ms)
#> ✔ Installed cellranger 1.1.0  (183ms)
#> ✔ Installed circlize 0.4.17  (208ms)
#> ✔ Installed clipr 0.8.0  (125ms)
#> ✔ Installed clue 0.3-66  (123ms)
#> ✔ Installed coda 0.19-4.1  (119ms)
#> ✔ Installed conflicted 1.2.0  (121ms)
#> ✔ Installed DBI 1.2.3  (214ms)
#> ✔ Installed dbplyr 2.5.1  (217ms)
#> ✔ Installed doParallel 1.0.17  (148ms)
#> ✔ Installed dtplyr 1.3.2  (95ms)
#> ✔ Installed foreach 1.5.2  (116ms)
#> ✔ Installed gargle 1.6.0  (117ms)
#> ✔ Installed GetoptLong 1.1.0  (207ms)
#> ✔ Installed gg.gap 1.3  (120ms)
#> ✔ Installed ggalluvial 0.12.5  (113ms)
#> ✔ Installed ggnetwork 0.5.14  (116ms)
#> ✔ Installed GlobalOptions 0.1.3  (118ms)
#> ✔ Installed googledrive 2.1.2  (187ms)
#> ✔ Installed googlesheets4 1.1.2  (211ms)
#> ✔ Installed gridBase 0.4-7  (120ms)
#> ✔ Installed haven 2.5.5  (120ms)
#> ✔ Built assorthead 1.4.0 (2.9s)
#> ! Failed to add assorthead 1.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed hms 1.1.4  (139ms)
#> ✔ Installed ids 1.0.1  (155ms)
#> ✔ Installed iterators 1.0.14  (139ms)
#> ✔ Installed lubridate 1.9.4  (85ms)
#> ✔ Installed network 1.19.0  (88ms)
#> ✔ Installed ragg 1.5.0  (87ms)
#> ✔ Installed readr 2.1.6  (87ms)
#> ✔ Installed readxl 1.4.5  (149ms)
#> ✔ Installed registry 0.5-1  (99ms)
#> ✔ Installed rematch 2.0.0  (76ms)
#> ✔ Installed rematch2 2.1.2  (79ms)
#> ✔ Installed reprex 2.1.1  (83ms)
#> ✔ Installed rjson 0.2.23  (78ms)
#> ✔ Installed rngtools 1.5.2  (155ms)
#> ℹ Building NMF 0.28
#> ✔ Installed roxygen2 7.3.3  (123ms)
#> ✔ Installed rstudioapi 0.17.1  (49ms)
#> ✔ Installed rvest 1.0.5  (33ms)
#> ✔ Installed selectr 0.5-0  (57ms)
#> ✔ Installed shape 1.4.6.1  (48ms)
#> ℹ Building ComplexHeatmap 2.26.0
#> ✔ Installed sna 2.8  (81ms)
#> ✔ Installed statnet.common 4.12.0  (36ms)
#> ✔ Installed svglite 2.2.2  (1s)
#> ✔ Installed systemfonts 1.3.1  (76ms)
#> ✔ Installed textshaping 1.0.4  (47ms)
#> ✔ Installed tidyverse 2.0.0  (39ms)
#> ✔ Installed timechange 0.3.0  (1s)
#> ✔ Built wordcloud 2.6 (7.6s)
#> ! Failed to add wordcloud 2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed tzdb 0.5.0  (104ms)
#> ✔ Installed uuid 1.2-1  (109ms)
#> ✔ Installed vroom 1.6.7  (116ms)
#> ✔ Installed wordcloud 2.6  (200ms)
#> ✔ Installed xml2 1.5.1  (212ms)
#> ✔ Installed assorthead 1.4.0  (125ms)
#> ℹ Building BiocNeighbors 2.4.0
#> ✔ Built ComplexHeatmap 2.26.0 (16.2s)
#> ! Failed to add ComplexHeatmap 2.26.0 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed ComplexHeatmap 2.26.0  (47ms)
#> ✔ Built NMF 0.28 (17s)
#> ! Failed to add NMF 0.28 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed NMF 0.28  (49ms)
#> ✔ Built BiocNeighbors 2.4.0 (31.6s)
#> ! Failed to add BiocNeighbors 2.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed BiocNeighbors 2.4.0  (89ms)
#> ℹ Packaging CellChat 2.2.0
#> ✔ Packaged CellChat 2.2.0 (4s)
#> ℹ Building CellChat 2.2.0
#> ✔ Built CellChat 2.2.0 (23.8s)
#> ✔ Installed CellChat 2.2.0 (github::jinworks/CellChat@623f48f) (63ms)
#> ✔ 1 pkg + 259 deps: kept 195, added 62, dld 60 (NA B) [1m 24.2s]
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Ngn3-high-EP, Endocrine, Ngn3-low-EP, Pre-endocrine 
#> The number of highly variable ligand-receptor pairs used for signaling inference is 841 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-12-16 07:26:19.277345]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-12-16 07:27:17.195971]"

CellChatPlot(pancreas_sub)

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

```
