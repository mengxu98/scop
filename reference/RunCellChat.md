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
#> → Will install 51 packages.
#> → Will download 1 CRAN package (1.68 MB), cached: 50 (0 B).
#> + BiocNeighbors    2.4.0     [bld][cmp]
#> + CellChat         2.2.0     [bld][cmp] (GitHub: 623f48f)
#> + ComplexHeatmap   2.26.0    [bld]
#> + GetoptLong       1.1.0      + ✔ perl
#> + GlobalOptions    0.1.3     
#> + NMF              0.28      [bld][cmp][dl] (1.68 MB)
#> + RcppEigen        0.3.4.0.2 
#> + assorthead       1.4.0     [bld]
#> + brew             1.0-10    
#> + bsicons          0.1.2     
#> + cellranger       1.1.0     
#> + circlize         0.4.16    
#> + clue             0.3-66    
#> + coda             0.19-4.1  
#> + conflicted       1.2.0     
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
#> + haven            2.5.5      + ✔ make, ✔ zlib1g-dev
#> + ids              1.0.1     
#> + iterators        1.0.14    
#> + lubridate        1.9.4     
#> + network          1.19.0    
#> + ragg             1.5.0      + ✔ libfreetype6-dev, ✔ libjpeg-dev, ✔ libpng-dev, ✔ libtiff-dev, ✔ libwebp-dev
#> + readxl           1.4.5     
#> + registry         0.5-1     
#> + rematch          2.0.0     
#> + rematch2         2.1.2     
#> + reprex           2.1.1      + ✔ pandoc
#> + rjson            0.2.23    
#> + rngtools         1.5.2     
#> + roxygen2         7.3.3     
#> + rvest            1.0.5     
#> + selectr          0.5-0     
#> + shape            1.4.6.1   
#> + sna              2.8       
#> + statnet.common   4.12.0    
#> + svglite          2.2.2      + ✔ libpng-dev
#> + textshaping      1.0.4      + ✔ libfreetype6-dev, ✔ libfribidi-dev, ✔ libharfbuzz-dev
#> + tidyverse        2.0.0     
#> + timechange       0.3.0     
#> + uuid             1.2-1     
#> + wordcloud        2.6       
#> + xml2             1.5.1      + ✔ libxml2-dev
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 1 pkg (1.68 MB), 50 cached
#> ✔ Got brew 1.0-10 (x86_64-pc-linux-gnu-ubuntu-24.04) (76.35 kB)
#> ✔ Got cellranger 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (103.66 kB)
#> ✔ Got BiocNeighbors 2.4.0 (source) (338.98 kB)
#> ✔ Got foreach 1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (139.65 kB)
#> ✔ Got bsicons 0.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (253.96 kB)
#> ✔ Got selectr 0.5-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (570.89 kB)
#> ✔ Got coda 0.19-4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (327.75 kB)
#> ✔ Got googlesheets4 1.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (515.22 kB)
#> ✔ Got rematch 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (16.20 kB)
#> ✔ Got rjson 0.2.23 (x86_64-pc-linux-gnu-ubuntu-24.04) (113.12 kB)
#> ✔ Got rvest 1.0.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (301.36 kB)
#> ✔ Got ComplexHeatmap 2.26.0 (source) (1.47 MB)
#> ✔ Got textshaping 1.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (189.59 kB)
#> ✔ Got statnet.common 4.12.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (313.83 kB)
#> ✔ Got uuid 1.2-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (48.94 kB)
#> ✔ Got roxygen2 7.3.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (729.96 kB)
#> ✔ Got assorthead 1.4.0 (source) (1.91 MB)
#> ✔ Got wordcloud 2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (107.74 kB)
#> ✔ Got haven 2.5.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (393.30 kB)
#> ✔ Got clue 0.3-66 (x86_64-pc-linux-gnu-ubuntu-24.04) (998.14 kB)
#> ✔ Got GetoptLong 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.41 MB)
#> ✔ Got registry 0.5-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (194.47 kB)
#> ✔ Got rematch2 2.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (46.02 kB)
#> ✔ Got gridBase 0.4-7 (x86_64-pc-linux-gnu-ubuntu-24.04) (162.00 kB)
#> ✔ Got reprex 2.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (495.27 kB)
#> ✔ Got tidyverse 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (425.82 kB)
#> ✔ Got ids 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (122.46 kB)
#> ✔ Got googledrive 2.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.92 MB)
#> ✔ Got doParallel 1.0.17 (x86_64-pc-linux-gnu-ubuntu-24.04) (188.11 kB)
#> ✔ Got lubridate 1.9.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (990.67 kB)
#> ✔ Got gg.gap 1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (20.48 kB)
#> ✔ Got circlize 0.4.16 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.32 MB)
#> ✔ Got rngtools 1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (78.03 kB)
#> ✔ Got gargle 1.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (618.93 kB)
#> ✔ Got sna 2.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.28 MB)
#> ✔ Got timechange 0.3.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (174.23 kB)
#> ✔ Got conflicted 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (55.04 kB)
#> ✔ Got iterators 1.0.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (346.57 kB)
#> ✔ Got ggalluvial 0.12.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.66 MB)
#> ✔ Got RcppEigen 0.3.4.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.85 MB)
#> ✔ Got svglite 2.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (239.34 kB)
#> ✔ Got readxl 1.4.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (412.32 kB)
#> ✔ Got xml2 1.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (282.14 kB)
#> ✔ Got shape 1.4.6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (749.97 kB)
#> ✔ Got dtplyr 1.3.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (409.84 kB)
#> ✔ Got network 1.19.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (823.12 kB)
#> ✔ Got ragg 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (743.50 kB)
#> ✔ Got GlobalOptions 0.1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (483.55 kB)
#> ✔ Got NMF 0.28 (source) (1.68 MB)
#> ✔ Got ggnetwork 0.5.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.97 MB)
#> ✔ Got CellChat 2.2.0 (source) (29.82 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install perl make zlib1g-dev libfreetype6-dev libjpeg-dev libpng-dev libtiff-dev libwebp-dev pandoc libfribidi-dev libharfbuzz-dev libxml2-dev libx11-dev libcurl4-openssl-dev libssl-dev libglpk-dev cmake python3 libicu-dev libfontconfig1-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> perl is already the newest version (5.38.2-3.2ubuntu0.2).
#> make is already the newest version (4.3-4.1build2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libjpeg-dev is already the newest version (8c-2ubuntu11).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> libtiff-dev is already the newest version (4.5.1+git230720-4ubuntu2.4).
#> libwebp-dev is already the newest version (1.3.2-0.4build3).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libfribidi-dev is already the newest version (1.0.13-3build1).
#> libharfbuzz-dev is already the newest version (8.3.0-2build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libx11-dev is already the newest version (2:1.8.7-1build1).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> libglpk-dev is already the newest version (5.0-1build2).
#> cmake is already the newest version (3.28.3-1build7).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> 0 upgraded, 0 newly installed, 0 to remove and 49 not upgraded.
#> ℹ Building assorthead 1.4.0
#> ✔ Installed brew 1.0-10  (85ms)
#> ✔ Installed bsicons 0.1.2  (111ms)
#> ✔ Installed cellranger 1.1.0  (201ms)
#> ✔ Installed circlize 0.4.16  (103ms)
#> ✔ Installed clue 0.3-66  (77ms)
#> ✔ Installed coda 0.19-4.1  (74ms)
#> ✔ Installed conflicted 1.2.0  (75ms)
#> ✔ Installed doParallel 1.0.17  (78ms)
#> ✔ Installed dtplyr 1.3.2  (139ms)
#> ✔ Installed foreach 1.5.2  (97ms)
#> ✔ Installed gargle 1.6.0  (84ms)
#> ✔ Installed GetoptLong 1.1.0  (82ms)
#> ✔ Installed gg.gap 1.3  (77ms)
#> ✔ Installed ggalluvial 0.12.5  (84ms)
#> ✔ Installed ggnetwork 0.5.14  (145ms)
#> ✔ Installed GlobalOptions 0.1.3  (112ms)
#> ✔ Installed googledrive 2.1.2  (92ms)
#> ✔ Installed googlesheets4 1.1.2  (90ms)
#> ✔ Installed gridBase 0.4-7  (116ms)
#> ✔ Installed haven 2.5.5  (143ms)
#> ✔ Installed ids 1.0.1  (104ms)
#> ✔ Installed iterators 1.0.14  (79ms)
#> ✔ Installed lubridate 1.9.4  (82ms)
#> ✔ Installed network 1.19.0  (84ms)
#> ✔ Installed ragg 1.5.0  (92ms)
#> ✔ Installed RcppEigen 0.3.4.0.2  (167ms)
#> ✔ Installed readxl 1.4.5  (95ms)
#> ✔ Installed registry 0.5-1  (78ms)
#> ✔ Installed rematch 2.0.0  (72ms)
#> ✔ Installed rematch2 2.1.2  (75ms)
#> ✔ Installed reprex 2.1.1  (115ms)
#> ✔ Installed rjson 0.2.23  (131ms)
#> ✔ Installed rngtools 1.5.2  (78ms)
#> ℹ Building NMF 0.28
#> ✔ Installed roxygen2 7.3.3  (94ms)
#> ✔ Installed rvest 1.0.5  (1s)
#> ✔ Built assorthead 1.4.0 (3.7s)
#> ! Failed to add assorthead 1.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed selectr 0.5-0  (74ms)
#> ✔ Installed shape 1.4.6.1  (134ms)
#> ℹ Building ComplexHeatmap 2.26.0
#> ✔ Installed sna 2.8  (100ms)
#> ✔ Installed statnet.common 4.12.0  (59ms)
#> ✔ Installed svglite 2.2.2  (1s)
#> ✔ Installed textshaping 1.0.4  (54ms)
#> ✔ Installed tidyverse 2.0.0  (113ms)
#> ✔ Installed timechange 0.3.0  (128ms)
#> ✔ Installed uuid 1.2-1  (120ms)
#> ✔ Installed wordcloud 2.6  (195ms)
#> ✔ Installed xml2 1.5.1  (135ms)
#> ✔ Installed assorthead 1.4.0  (125ms)
#> ℹ Building BiocNeighbors 2.4.0
#> ✔ Built ComplexHeatmap 2.26.0 (19.3s)
#> ! Failed to add ComplexHeatmap 2.26.0 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed ComplexHeatmap 2.26.0  (61ms)
#> ✔ Built NMF 0.28 (21.2s)
#> ! Failed to add NMF 0.28 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed NMF 0.28  (1s)
#> ✔ Built BiocNeighbors 2.4.0 (30.4s)
#> ! Failed to add BiocNeighbors 2.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) to the cache
#> ✔ Installed BiocNeighbors 2.4.0  (80ms)
#> ℹ Packaging CellChat 2.2.0
#> ✔ Packaged CellChat 2.2.0 (4.1s)
#> ℹ Building CellChat 2.2.0
#> ✔ Built CellChat 2.2.0 (24.3s)
#> ✔ Installed CellChat 2.2.0 (github::jinworks/CellChat@623f48f) (60ms)
#> ✔ 1 pkg + 260 deps: kept 210, added 51, dld 51 (NA B) [1m 21.9s]
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Ngn3-high-EP, Endocrine, Ngn3-low-EP, Pre-endocrine 
#> The number of highly variable ligand-receptor pairs used for signaling inference is 841 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-12-05 08:42:31.257924]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-12-05 08:43:39.994071]"

CellChatPlot(pancreas_sub)

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

```
