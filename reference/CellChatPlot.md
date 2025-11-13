# Plot CellChat analysis results

CellChatPlot creates various visualizations for CellChat analysis
results stored in a Seurat object.

## Usage

``` r
CellChatPlot(
  srt,
  plot_type = "aggregate",
  condition = NULL,
  pathway = NULL,
  dirpath = NULL,
  output_format = "pdf",
  top_n = 10,
  base_height = 1,
  base_width = 1,
  res = 300,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object that has been processed with
  [RunCellChat](https://mengxu98.github.io/scop/reference/RunCellChat.md).

- plot_type:

  Type of plot to create. Options: `"aggregate"`, `"pathway"`,
  `"comparison"`, `"heatmap"`, `"circle"`, `"bubble"`, `"gene"`.

- condition:

  Condition to plot (if multiple conditions exist).

- pathway:

  Specific pathway to visualize (for pathway, bubble, and gene plots).
  If `NULL`, uses top pathways.

- dirpath:

  Directory to save plots.

- output_format:

  Format of output figure: `"png"` or `"pdf"`. Default: `"png"`.

- top_n:

  Number of top pathways to use for plotting. Default: 10.

- base_height:

  Base height multiplier for all plots. Default: 1.

- base_width:

  Base width multiplier for all plots. Default: 1.

- res:

  Resolution for PNG output. Default: 300.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[RunCellChat](https://mengxu98.github.io/scop/reference/RunCellChat.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 11:35:03] Start standard scop workflow...
#> ℹ [2025-11-13 11:35:05] Checking a list of <Seurat> object...
#> ! [2025-11-13 11:35:05] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 11:35:05] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 11:35:07] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 11:35:08] Use the separate HVF from srt_list
#> ℹ [2025-11-13 11:35:08] Number of available HVF: 2000
#> ℹ [2025-11-13 11:35:08] Finished check
#> ℹ [2025-11-13 11:35:09] Perform `Seurat::ScaleData()`
#> Warning: Different features in new layer data than already exists for scale.data
#> ℹ [2025-11-13 11:35:10] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 11:35:11] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 11:35:11] Reorder clusters...
#> First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#> This message is displayed once every 8 hours.
#> ℹ [2025-11-13 11:35:11] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 11:35:11] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 11:35:11] UMAP will return its model
#> ℹ [2025-11-13 11:35:13] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 11:35:13] UMAP will return its model
#> ✔ [2025-11-13 11:35:15] Run scop standard workflow done
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "mouse"
)
#> ℹ [2025-11-13 11:35:15] Start CellChat analysis
#> ◌ [2025-11-13 11:35:15] Installing: CellChat...
#> ℹ Loading metadata database
#> ✔ Loading metadata database ... done
#> 
#>  
#> → Will install 151 packages.
#> → Will download 1 CRAN package (1.68 MB), cached: 150 (0 B).
#> + Biobase          2.70.0     [bld][cmp]
#> + BiocGenerics     0.56.0     [bld]
#> + BiocManager      1.30.26    
#> + BiocNeighbors    2.4.0      [bld][cmp]
#> + CellChat         2.2.0      [bld][cmp] (GitHub: 623f48f)
#> + ComplexHeatmap   2.26.0     [bld]
#> + Deriv            4.2.0      
#> + FNN              1.1.4.1    
#> + Formula          1.2-5      
#> + GetoptLong       1.0.5       + ✔ perl
#> + GlobalOptions    0.1.2      
#> + IRanges          2.44.0     [bld][cmp]
#> + MatrixModels     0.5-4      
#> + NMF              0.28       [bld][cmp][dl] (1.68 MB)
#> + R6               2.6.1      
#> + RColorBrewer     1.1-3      
#> + RSpectra         0.16-2     
#> + Rcpp             1.1.0      
#> + RcppEigen        0.3.4.0.2  
#> + RcppTOML         0.2.3      
#> + Rdpack           2.6.4      
#> + S4Vectors        0.48.0     [bld][cmp]
#> + S7               0.2.0      
#> + SparseM          1.84-2     
#> + abind            1.4-8      
#> + askpass          1.2.1      
#> + assorthead       1.4.0      [bld]
#> + backports        1.5.0      
#> + base64enc        0.1-3      
#> + broom            1.0.10     
#> + bslib            0.9.0      
#> + cachem           1.1.0      
#> + car              3.1-3      
#> + carData          3.0-5      
#> + circlize         0.4.16     
#> + cli              3.6.5      
#> + clue             0.3-66     
#> + coda             0.19-4.1   
#> + colorspace       2.1-2      
#> + commonmark       2.0.0      
#> + corrplot         0.95       
#> + cowplot          1.2.0      
#> + crayon           1.5.3      
#> + crosstalk        1.2.2      
#> + curl             7.0.0       + ✔ libcurl4-openssl-dev, ✔ libssl-dev
#> + data.table       1.17.8     
#> + digest           0.6.37     
#> + doBy             4.7.0      
#> + doParallel       1.0.17     
#> + dplyr            1.1.4      
#> + evaluate         1.0.5      
#> + farver           2.1.2      
#> + fastmap          1.2.0      
#> + fontawesome      0.5.3      
#> + foreach          1.5.2      
#> + fs               1.6.6       + ✔ make
#> + future           1.67.0     
#> + future.apply     1.20.0     
#> + generics         0.1.4      
#> + ggalluvial       0.12.5     
#> + ggnetwork        0.5.14     
#> + ggplot2          4.0.0      
#> + ggpubr           0.6.2      
#> + ggrepel          0.9.6      
#> + ggsci            4.1.0      
#> + ggsignif         0.6.4      
#> + globals          0.18.0     
#> + glue             1.8.0      
#> + gridBase         0.4-7      
#> + gridExtra        2.3        
#> + gtable           0.3.6      
#> + here             1.0.2      
#> + highr            0.11       
#> + htmltools        0.5.8.1    
#> + htmlwidgets      1.6.4      
#> + httpuv           1.6.16      + ✔ make, ✔ zlib1g-dev
#> + httr             1.4.7      
#> + igraph           2.2.1       + ✔ libglpk-dev, ✔ libxml2-dev
#> + irlba            2.3.5.1    
#> + isoband          0.2.7      
#> + iterators        1.0.14     
#> + jquerylib        0.1.4      
#> + jsonlite         2.0.0      
#> + knitr            1.50        + ✔ pandoc
#> + labeling         0.4.3      
#> + later            1.4.4      
#> + lazyeval         0.2.2      
#> + lifecycle        1.0.4      
#> + listenv          0.10.0     
#> + lme4             1.1-37     
#> + magrittr         2.0.4      
#> + matrixStats      1.5.0      
#> + memoise          2.0.1      
#> + microbenchmark   1.5.0      
#> + mime             0.13       
#> + minqa            1.2.8       + ✔ make
#> + modelr           0.1.11     
#> + network          1.19.0     
#> + nloptr           2.2.1       + ✔ cmake
#> + numDeriv         2016.8-1.1 
#> + openssl          2.3.4       + ✔ libssl-dev
#> + otel             0.2.0      
#> + parallelly       1.45.1     
#> + patchwork        1.3.2      
#> + pbapply          1.7-4      
#> + pbkrtest         0.5.5      
#> + pillar           1.11.1     
#> + pkgconfig        2.0.3      
#> + plotly           4.11.0     
#> + plyr             1.8.9      
#> + png              0.1-8       + ✔ libpng-dev
#> + polynom          1.4-1      
#> + promises         1.5.0      
#> + purrr            1.2.0      
#> + quantreg         6.1        
#> + rappdirs         0.3.3      
#> + rbibutils        2.4        
#> + reformulas       0.4.2      
#> + registry         0.5-1      
#> + reshape2         1.4.4      
#> + reticulate       1.44.0      + ✔ python3
#> + rjson            0.2.23     
#> + rlang            1.1.6      
#> + rmarkdown        2.30        + ✔ pandoc
#> + rngtools         1.5.2      
#> + rprojroot        2.1.1      
#> + rstatix          0.7.3      
#> + sass             0.4.10      + ✔ make
#> + scales           1.4.0      
#> + shape            1.4.6.1    
#> + shiny            1.11.1     
#> + sna              2.8        
#> + sourcetools      0.1.7-1    
#> + statnet.common   4.12.0     
#> + stringi          1.8.7       + ✔ libicu-dev
#> + stringr          1.6.0      
#> + svglite          2.2.2       + ✔ libpng-dev
#> + sys              3.4.3      
#> + systemfonts      1.3.1       + ✔ libfontconfig1-dev, ✔ libfreetype6-dev
#> + textshaping      1.0.4       + ✔ libfreetype6-dev, ✔ libfribidi-dev, ✔ libharfbuzz-dev
#> + tibble           3.3.0      
#> + tidyr            1.3.1      
#> + tidyselect       1.2.1      
#> + tinytex          0.57       
#> + utf8             1.2.6      
#> + vctrs            0.6.5      
#> + viridisLite      0.4.2      
#> + withr            3.0.2      
#> + xfun             0.54       
#> + xtable           1.8-4      
#> + yaml             2.3.10     
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 1 pkg (1.68 MB), 150 cached
#> ✔ Cached copy of NMF 0.28 (source) is the latest build
#> ✔ Got CellChat 2.2.0 (source) (29.82 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libcurl4-openssl-dev libssl-dev make perl zlib1g-dev libglpk-dev libxml2-dev pandoc cmake libpng-dev python3 libicu-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> perl is already the newest version (5.38.2-3.2ubuntu0.2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> pandoc is already the newest version (3.1.3+ds-2).
#> cmake is already the newest version (3.28.3-1build7).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libfribidi-dev is already the newest version (1.0.13-3build1).
#> libharfbuzz-dev is already the newest version (8.3.0-2build2).
#> 0 upgraded, 0 newly installed, 0 to remove and 23 not upgraded.
#> ℹ Building assorthead 1.4.0
#> ✔ Installed abind 1.4-8  (79ms)
#> ✔ Installed askpass 1.2.1  (132ms)
#> ✔ Installed backports 1.5.0  (180ms)
#> ✔ Installed base64enc 0.1-3  (72ms)
#> ✔ Installed BiocManager 1.30.26  (72ms)
#> ✔ Installed broom 1.0.10  (103ms)
#> ✔ Installed bslib 0.9.0  (1.1s)
#> ✔ Installed cachem 1.1.0  (1.1s)
#> ✔ Installed car 3.1-3  (76ms)
#> ✔ Installed carData 3.0-5  (76ms)
#> ✔ Installed circlize 0.4.16  (77ms)
#> ✔ Installed cli 3.6.5  (74ms)
#> ✔ Installed clue 0.3-66  (110ms)
#> ✔ Installed coda 0.19-4.1  (121ms)
#> ✔ Installed colorspace 2.1-2  (81ms)
#> ✔ Installed commonmark 2.0.0  (73ms)
#> ✔ Installed corrplot 0.95  (77ms)
#> ✔ Installed cowplot 1.2.0  (73ms)
#> ✔ Installed crayon 1.5.3  (71ms)
#> ✔ Installed crosstalk 1.2.2  (119ms)
#> ✔ Installed curl 7.0.0  (98ms)
#> ✔ Installed data.table 1.17.8  (89ms)
#> ✔ Installed Deriv 4.2.0  (77ms)
#> ✔ Installed digest 0.6.37  (75ms)
#> ✔ Installed doBy 4.7.0  (91ms)
#> ✔ Built assorthead 1.4.0 (2.7s)
#> ✔ Installed doParallel 1.0.17  (155ms)
#> ✔ Installed dplyr 1.1.4  (157ms)
#> ✔ Installed evaluate 1.0.5  (70ms)
#> ✔ Installed farver 2.1.2  (66ms)
#> ✔ Installed fastmap 1.2.0  (64ms)
#> ✔ Installed FNN 1.1.4.1  (64ms)
#> ✔ Installed fontawesome 0.5.3  (65ms)
#> ✔ Installed foreach 1.5.2  (122ms)
#> ✔ Installed Formula 1.2-5  (79ms)
#> ✔ Installed fs 1.6.6  (67ms)
#> ✔ Installed future 1.67.0  (72ms)
#> ✔ Installed future.apply 1.20.0  (69ms)
#> ✔ Installed generics 0.1.4  (66ms)
#> ℹ Building BiocGenerics 0.56.0
#> ✔ Installed GetoptLong 1.0.5  (122ms)
#> ✔ Installed ggalluvial 0.12.5  (38ms)
#> ✔ Installed ggnetwork 0.5.14  (1s)
#> ✔ Installed ggpubr 0.6.2  (1s)
#> ✔ Installed ggplot2 4.0.0  (1.1s)
#> ✔ Installed ggrepel 0.9.6  (82ms)
#> ✔ Installed ggsci 4.1.0  (76ms)
#> ✔ Installed ggsignif 0.6.4  (71ms)
#> ✔ Installed GlobalOptions 0.1.2  (82ms)
#> ✔ Installed globals 0.18.0  (99ms)
#> ✔ Installed glue 1.8.0  (108ms)
#> ✔ Installed gridBase 0.4-7  (68ms)
#> ✔ Installed gridExtra 2.3  (67ms)
#> ✔ Installed gtable 0.3.6  (70ms)
#> ✔ Built BiocGenerics 0.56.0 (3s)
#> ✔ Installed here 1.0.2  (93ms)
#> ✔ Installed highr 0.11  (88ms)
#> ✔ Installed htmltools 0.5.8.1  (64ms)
#> ✔ Installed htmlwidgets 1.6.4  (95ms)
#> ✔ Installed httpuv 1.6.16  (68ms)
#> ✔ Installed httr 1.4.7  (67ms)
#> ✔ Installed irlba 2.3.5.1  (1s)
#> ✔ Installed igraph 2.2.1  (1.1s)
#> ✔ Installed isoband 0.2.7  (69ms)
#> ✔ Installed iterators 1.0.14  (64ms)
#> ✔ Installed jquerylib 0.1.4  (96ms)
#> ✔ Installed jsonlite 2.0.0  (68ms)
#> ✔ Installed knitr 1.50  (77ms)
#> ✔ Installed labeling 0.4.3  (74ms)
#> ✔ Installed later 1.4.4  (65ms)
#> ✔ Installed lazyeval 0.2.2  (65ms)
#> ✔ Installed lifecycle 1.0.4  (65ms)
#> ✔ Installed listenv 0.10.0  (65ms)
#> ✔ Installed lme4 1.1-37  (97ms)
#> ✔ Installed magrittr 2.0.4  (69ms)
#> ✔ Installed MatrixModels 0.5-4  (66ms)
#> ✔ Installed matrixStats 1.5.0  (67ms)
#> ✔ Installed memoise 2.0.1  (66ms)
#> ✔ Installed microbenchmark 1.5.0  (65ms)
#> ✔ Installed mime 0.13  (65ms)
#> ✔ Installed minqa 1.2.8  (93ms)
#> ✔ Installed modelr 0.1.11  (65ms)
#> ✔ Installed network 1.19.0  (67ms)
#> ✔ Installed nloptr 2.2.1  (67ms)
#> ✔ Installed NMF 0.28  (68ms)
#> ✔ Installed numDeriv 2016.8-1.1  (64ms)
#> ✔ Installed openssl 2.3.4  (63ms)
#> ✔ Installed otel 0.2.0  (97ms)
#> ✔ Installed parallelly 1.45.1  (68ms)
#> ✔ Installed patchwork 1.3.2  (69ms)
#> ✔ Installed pbapply 1.7-4  (67ms)
#> ✔ Installed pbkrtest 0.5.5  (66ms)
#> ✔ Installed pillar 1.11.1  (67ms)
#> ✔ Installed pkgconfig 2.0.3  (66ms)
#> ✔ Installed plotly 4.11.0  (1.1s)
#> ✔ Installed plyr 1.8.9  (1.1s)
#> ✔ Installed png 0.1-8  (65ms)
#> ✔ Installed polynom 1.4-1  (63ms)
#> ✔ Installed promises 1.5.0  (66ms)
#> ✔ Installed purrr 1.2.0  (68ms)
#> ✔ Installed quantreg 6.1  (66ms)
#> ✔ Installed R6 2.6.1  (95ms)
#> ✔ Installed rappdirs 0.3.3  (96ms)
#> ✔ Installed rbibutils 2.4  (69ms)
#> ✔ Installed RColorBrewer 1.1-3  (63ms)
#> ✔ Installed Rcpp 1.1.0  (71ms)
#> ✔ Installed RcppEigen 0.3.4.0.2  (86ms)
#> ✔ Installed RcppTOML 0.2.3  (64ms)
#> ✔ Installed Rdpack 2.6.4  (63ms)
#> ✔ Installed reformulas 0.4.2  (120ms)
#> ✔ Installed registry 0.5-1  (60ms)
#> ✔ Installed reshape2 1.4.4  (60ms)
#> ✔ Installed reticulate 1.44.0  (66ms)
#> ✔ Installed rjson 0.2.23  (64ms)
#> ✔ Installed rlang 1.1.6  (63ms)
#> ✔ Installed rmarkdown 2.30  (72ms)
#> ✔ Installed rngtools 1.5.2  (104ms)
#> ✔ Installed rprojroot 2.1.1  (66ms)
#> ✔ Installed RSpectra 0.16-2  (66ms)
#> ✔ Installed rstatix 0.7.3  (69ms)
#> ✔ Installed S7 0.2.0  (69ms)
#> ✔ Installed sass 0.4.10  (80ms)
#> ✔ Installed scales 1.4.0  (67ms)
#> ✔ Installed shape 1.4.6.1  (97ms)
#> ✔ Installed shiny 1.11.1  (103ms)
#> ✔ Installed sna 2.8  (70ms)
#> ✔ Installed sourcetools 0.1.7-1  (63ms)
#> ✔ Installed SparseM 1.84-2  (65ms)
#> ✔ Installed statnet.common 4.12.0  (63ms)
#> ✔ Installed stringr 1.6.0  (1s)
#> ✔ Installed stringi 1.8.7  (1.1s)
#> ✔ Installed svglite 2.2.2  (100ms)
#> ✔ Installed sys 3.4.3  (67ms)
#> ✔ Installed systemfonts 1.3.1  (67ms)
#> ✔ Installed textshaping 1.0.4  (71ms)
#> ✔ Installed tibble 3.3.0  (68ms)
#> ✔ Installed tidyr 1.3.1  (66ms)
#> ✔ Installed tidyselect 1.2.1  (65ms)
#> ✔ Installed tinytex 0.57  (99ms)
#> ✔ Installed utf8 1.2.6  (65ms)
#> ✔ Installed vctrs 0.6.5  (65ms)
#> ✔ Installed viridisLite 0.4.2  (66ms)
#> ✔ Installed withr 3.0.2  (65ms)
#> ✔ Installed xfun 0.54  (64ms)
#> ✔ Installed xtable 1.8-4  (64ms)
#> ✔ Installed yaml 2.3.10  (96ms)
#> ✔ Installed BiocGenerics 0.56.0  (1s)
#> ℹ Building Biobase 2.70.0
#> ✔ Installed assorthead 1.4.0  (1.1s)
#> ℹ Building BiocNeighbors 2.4.0
#> ℹ Building S4Vectors 0.48.0
#> ✔ Built Biobase 2.70.0 (9.8s)
#> ✔ Installed Biobase 2.70.0  (49ms)
#> ✔ Built S4Vectors 0.48.0 (18.4s)
#> ✔ Installed S4Vectors 0.48.0  (50ms)
#> ℹ Building IRanges 2.44.0
#> ✔ Built BiocNeighbors 2.4.0 (32.2s)
#> ✔ Installed BiocNeighbors 2.4.0  (80ms)
#> ✔ Built IRanges 2.44.0 (37s)
#> ✔ Installed IRanges 2.44.0  (1s)
#> ℹ Building ComplexHeatmap 2.26.0
#> ✔ Built ComplexHeatmap 2.26.0 (12.7s)
#> ✔ Installed ComplexHeatmap 2.26.0  (37ms)
#> ℹ Packaging CellChat 2.2.0
#> ✔ Packaged CellChat 2.2.0 (4.2s)
#> ℹ Building CellChat 2.2.0
#> ✔ Built CellChat 2.2.0 (24.6s)
#> ✔ Installed CellChat 2.2.0 (github::jinworks/CellChat@623f48f) (62ms)
#> ✔ 1 pkg + 160 deps: kept 9, added 151, dld 1 (NA B) [2m 15.8s]
#> ✔ [2025-11-13 11:37:31] jinworks/CellChat installed successfully
#> ◌ [2025-11-13 11:37:31] Installing: presto...
#>  
#> → Will install 2 packages.
#> → All 2 packages (0 B) are cached.
#> + RcppArmadillo   15.0.2-2 
#> + presto          1.0.0    [bld][cmp] (GitHub: 7636b3d)
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 2 pkgs are cached
#> ✔ Got presto 1.0.0 (source) (746.05 kB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 23 not upgraded.
#> ✔ Installed RcppArmadillo 15.0.2-2  (91ms)
#> ℹ Packaging presto 1.0.0
#> ✔ Packaged presto 1.0.0 (715ms)
#> ℹ Building presto 1.0.0
#> ✔ Built presto 1.0.0 (12.3s)
#> ✔ Installed presto 1.0.0 (github::immunogenomics/presto@7636b3d) (34ms)
#> ✔ 1 pkg + 24 deps: kept 23, added 2, dld 1 (NA B) [17.7s]
#> ✔ [2025-11-13 11:37:49] immunogenomics/presto installed successfully
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Ngn3-high-EP, Endocrine, Ngn3-low-EP, Pre-endocrine 
#> The number of highly variable ligand-receptor pairs used for signaling inference is 841 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-11-13 11:37:51.449952]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-11-13 11:39:02.205628]"
#> ✔ [2025-11-13 11:39:02] CellChat analysis completed

CellChatPlot(pancreas_sub, plot_type = "aggregate")
#> ℹ [2025-11-13 11:39:02] Creating "aggregate" plot for condition "ALL"

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

#> ✔ [2025-11-13 11:39:02] Plot creation completed


CellChatPlot(pancreas_sub, plot_type = "pathway")
#> ℹ [2025-11-13 11:39:03] Creating "pathway" plot for condition "ALL"
#> ℹ [2025-11-13 11:39:03] Using top 10 pathways
#> Comparing communications on a single object 

#> Comparing communications on a single object 
#> Comparing communications on a single object 

#> Comparing communications on a single object 
#> Comparing communications on a single object 

#> Comparing communications on a single object 
#> Comparing communications on a single object 

#> Comparing communications on a single object 
#> Comparing communications on a single object 

#> Comparing communications on a single object 
#> ✔ [2025-11-13 11:39:07] Plot creation completed
#> NULL

CellChatPlot(pancreas_sub, plot_type = "bubble")
#> ℹ [2025-11-13 11:39:07] Creating "bubble" plot for condition "ALL"
#> ℹ [2025-11-13 11:39:07] No pathway specified. Using top 3 pathways for bubble plots
#> Comparing communications on a single object 
#> Comparing communications on a single object 
#> Comparing communications on a single object 
#> ✔ [2025-11-13 11:39:07] Plot creation completed


CellChatPlot(pancreas_sub, plot_type = "gene")
#> ℹ [2025-11-13 11:39:08] Creating "gene" plot for condition "ALL"
#> ℹ [2025-11-13 11:39:08] No pathway specified. Using top 3 pathways for gene expression plots
#> ✔ [2025-11-13 11:39:08] Plot creation completed


CellChatPlot(pancreas_sub, plot_type = "heatmap")
#> ℹ [2025-11-13 11:39:10] Creating "heatmap" plot for condition "ALL"

#> ✔ [2025-11-13 11:39:12] Plot creation completed
#> NULL
```
