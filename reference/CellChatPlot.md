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

  Format of output figure: `"png"` or `"pdf"`. Default is `"png"`.

- top_n:

  Number of top pathways to use for plotting. Default is `10`.

- base_height:

  Base height multiplier for all plots. Default is `1`.

- base_width:

  Base width multiplier for all plots. Default is `1`.

- res:

  Resolution for PNG output. Default is `300`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[RunCellChat](https://mengxu98.github.io/scop/reference/RunCellChat.md)

## Examples

``` r
options(log_message.verbose = FALSE)
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
#> First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
#> This message is displayed once every 8 hours.
pancreas_sub <- RunCellChat(
  pancreas_sub,
  group.by = "CellType",
  species = "mouse"
)
#> ℹ Loading metadata database
#> ✔ Loading metadata database ... done
#> 
#>  
#> → Will install 151 packages.
#> → Will download 1 CRAN package (1.68 MB), cached: 150 (0 B).
#> + Biobase          2.70.0     [bld][cmp]
#> + BiocGenerics     0.56.0     [bld]
#> + BiocManager      1.30.27    
#> + BiocNeighbors    2.4.0      [bld][cmp]
#> + CellChat         2.2.0      [bld][cmp] (GitHub: 623f48f)
#> + ComplexHeatmap   2.26.0     [bld]
#> + Deriv            4.2.0      
#> + FNN              1.1.4.1    
#> + Formula          1.2-5      
#> + GetoptLong       1.1.0       + ✔ perl
#> + GlobalOptions    0.1.3      
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
#> + S7               0.2.1      
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
#> + digest           0.6.39     
#> + doBy             4.7.0      
#> + doParallel       1.0.17     
#> + dplyr            1.1.4      
#> + evaluate         1.0.5      
#> + farver           2.1.2      
#> + fastmap          1.2.0      
#> + fontawesome      0.5.3      
#> + foreach          1.5.2      
#> + fs               1.6.6       + ✔ make
#> + future           1.68.0     
#> + future.apply     1.20.0     
#> + generics         0.1.4      
#> + ggalluvial       0.12.5     
#> + ggnetwork        0.5.14     
#> + ggplot2          4.0.1      
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
#> + reshape2         1.4.5      
#> + reticulate       1.44.1      + ✔ python3
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
#> + tinytex          0.58       
#> + utf8             1.2.6      
#> + vctrs            0.6.5      
#> + viridisLite      0.4.2      
#> + withr            3.0.2      
#> + xfun             0.54       
#> + xtable           1.8-4      
#> + yaml             2.3.11     
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 1 pkg (1.68 MB), 150 cached
#> ✔ Got BiocGenerics 0.56.0 (source) (61.53 kB)
#> ✔ Got base64enc 0.1-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (26.57 kB)
#> ✔ Got askpass 1.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (21.89 kB)
#> ✔ Got abind 1.4-8 (x86_64-pc-linux-gnu-ubuntu-24.04) (64.92 kB)
#> ✔ Got backports 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (115.90 kB)
#> ✔ Got BiocNeighbors 2.4.0 (source) (338.98 kB)
#> ✔ Got lazyeval 0.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (157.44 kB)
#> ✔ Got IRanges 2.44.0 (source) (496.06 kB)
#> ✔ Got BiocManager 1.30.27 (x86_64-pc-linux-gnu-ubuntu-24.04) (666.57 kB)
#> ✔ Got parallelly 1.45.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (571.86 kB)
#> ✔ Got png 0.1-8 (x86_64-pc-linux-gnu-ubuntu-24.04) (40.57 kB)
#> ✔ Got S4Vectors 0.48.0 (source) (843.90 kB)
#> ✔ Got ComplexHeatmap 2.26.0 (source) (1.47 MB)
#> ✔ Got Rdpack 2.6.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (637.80 kB)
#> ✔ Got RSpectra 0.16-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (529.66 kB)
#> ✔ Got isoband 0.2.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.64 MB)
#> ✔ Got Deriv 4.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (152.58 kB)
#> ✔ Got Biobase 2.70.0 (source) (1.98 MB)
#> ✔ Got broom 1.0.10 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.92 MB)
#> ✔ Got crosstalk 1.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (412.02 kB)
#> ✔ Got assorthead 1.4.0 (source) (1.91 MB)
#> ✔ Got NMF 0.28 (source) (1.68 MB)
#> ✔ Got fs 1.6.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (310.07 kB)
#> ✔ Got htmlwidgets 1.6.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (815.05 kB)
#> ✔ Got dplyr 1.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.49 MB)
#> ✔ Got farver 2.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.47 MB)
#> ✔ Got fontawesome 0.5.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.40 MB)
#> ✔ Got ggrepel 0.9.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (290.23 kB)
#> ✔ Got jquerylib 0.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (526.85 kB)
#> ✔ Got labeling 0.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (60.95 kB)
#> ✔ Got later 1.4.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (155.03 kB)
#> ✔ Got memoise 2.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (48.86 kB)
#> ✔ Got pbapply 1.7-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (100.53 kB)
#> ✔ Got nloptr 2.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (567.65 kB)
#> ✔ Got ggsci 4.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.30 MB)
#> ✔ Got numDeriv 2016.8-1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (114.36 kB)
#> ✔ Got pillar 1.11.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (660.46 kB)
#> ✔ Got otel 0.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (279.62 kB)
#> ✔ Got R6 2.6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (86.81 kB)
#> ✔ Got polynom 1.4-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (406.99 kB)
#> ✔ Got purrr 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (553.55 kB)
#> ✔ Got reformulas 0.4.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (139.23 kB)
#> ✔ Got ggnetwork 0.5.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.97 MB)
#> ✔ Got rprojroot 2.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (113.23 kB)
#> ✔ Got sna 2.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.28 MB)
#> ✔ Got rlang 1.1.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.59 MB)
#> ✔ Got shape 1.4.6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (749.97 kB)
#> ✔ Got Rcpp 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.19 MB)
#> ✔ Got svglite 2.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (239.34 kB)
#> ✔ Got tidyselect 1.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (225.28 kB)
#> ✔ Got lme4 1.1-37 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.28 MB)
#> ✔ Got SparseM 1.84-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (887.98 kB)
#> ✔ Got tibble 3.3.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (680.35 kB)
#> ✔ Got ggplot2 4.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (8.47 MB)
#> ✔ Got vctrs 0.6.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.31 MB)
#> ✔ Got clue 0.3-66 (x86_64-pc-linux-gnu-ubuntu-24.04) (998.14 kB)
#> ✔ Got carData 3.0-5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.70 MB)
#> ✔ Got cachem 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (67.49 kB)
#> ✔ Got Formula 1.2-5 (x86_64-pc-linux-gnu-ubuntu-24.04) (159.13 kB)
#> ✔ Got coda 0.19-4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (327.75 kB)
#> ✔ Got stringi 1.8.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.29 MB)
#> ✔ Got doParallel 1.0.17 (x86_64-pc-linux-gnu-ubuntu-24.04) (188.11 kB)
#> ✔ Got fastmap 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (66.05 kB)
#> ✔ Got foreach 1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (139.65 kB)
#> ✔ Got ggsignif 0.6.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (602.07 kB)
#> ✔ Got globals 0.18.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (158.33 kB)
#> ✔ Got future 1.68.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (966.66 kB)
#> ✔ Got data.table 1.17.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.67 MB)
#> ✔ Got gridBase 0.4-7 (x86_64-pc-linux-gnu-ubuntu-24.04) (162.00 kB)
#> ✔ Got httr 1.4.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (486.52 kB)
#> ✔ Got highr 0.11 (x86_64-pc-linux-gnu-ubuntu-24.04) (37.50 kB)
#> ✔ Got corrplot 0.95 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.82 MB)
#> ✔ Got iterators 1.0.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (346.57 kB)
#> ✔ Got magrittr 2.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (221.44 kB)
#> ✔ Got colorspace 2.1-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.64 MB)
#> ✔ Got MatrixModels 0.5-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (408.50 kB)
#> ✔ Got microbenchmark 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (65.96 kB)
#> ✔ Got modelr 0.1.11 (x86_64-pc-linux-gnu-ubuntu-24.04) (200.70 kB)
#> ✔ Got pkgconfig 2.0.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (18.08 kB)
#> ✔ Got bslib 0.9.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.67 MB)
#> ✔ Got matrixStats 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (474.75 kB)
#> ✔ Got ggalluvial 0.12.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.66 MB)
#> ✔ Got network 1.19.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (823.12 kB)
#> ✔ Got plyr 1.8.9 (x86_64-pc-linux-gnu-ubuntu-24.04) (787.35 kB)
#> ✔ Got S7 0.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (325.15 kB)
#> ✔ Got ggpubr 0.6.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.13 MB)
#> ✔ Got quantreg 6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.46 MB)
#> ✔ Got promises 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.67 MB)
#> ✔ Got rngtools 1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (78.03 kB)
#> ✔ Got stringr 1.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (333.76 kB)
#> ✔ Got RcppEigen 0.3.4.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.85 MB)
#> ✔ Got sys 3.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (40.73 kB)
#> ✔ Got textshaping 1.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (189.59 kB)
#> ✔ Got statnet.common 4.12.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (313.83 kB)
#> ✔ Got scales 1.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (841.36 kB)
#> ✔ Got utf8 1.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (151.81 kB)
#> ✔ Got commonmark 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (148.07 kB)
#> ✔ Got crayon 1.5.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (163.30 kB)
#> ✔ Got plotly 4.11.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.86 MB)
#> ✔ Got reticulate 1.44.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.98 MB)
#> ✔ Got tidyr 1.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.18 MB)
#> ✔ Got sass 0.4.10 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.46 MB)
#> ✔ Got generics 0.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (80.38 kB)
#> ✔ Got here 1.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (53.59 kB)
#> ✔ Got FNN 1.1.4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (128.22 kB)
#> ✔ Got curl 7.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (788.30 kB)
#> ✔ Got listenv 0.10.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (106.64 kB)
#> ✔ Got glue 1.8.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (168.12 kB)
#> ✔ Got minqa 1.2.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (122.60 kB)
#> ✔ Got htmltools 0.5.8.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (354.46 kB)
#> ✔ Got pbkrtest 0.5.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (222.06 kB)
#> ✔ Got RColorBrewer 1.1-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (51.81 kB)
#> ✔ Got irlba 2.3.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (292.19 kB)
#> ✔ Got reshape2 1.4.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (122.33 kB)
#> ✔ Got rjson 0.2.23 (x86_64-pc-linux-gnu-ubuntu-24.04) (113.12 kB)
#> ✔ Got httpuv 1.6.16 (x86_64-pc-linux-gnu-ubuntu-24.04) (656.07 kB)
#> ✔ Got sourcetools 0.1.7-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (46.59 kB)
#> ✔ Got xfun 0.54 (x86_64-pc-linux-gnu-ubuntu-24.04) (583.00 kB)
#> ✔ Got systemfonts 1.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (816.45 kB)
#> ✔ Got rbibutils 2.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.15 MB)
#> ✔ Got evaluate 1.0.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (102.86 kB)
#> ✔ Got cli 3.6.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.34 MB)
#> ✔ Got circlize 0.4.16 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.32 MB)
#> ✔ Got lifecycle 1.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (125.07 kB)
#> ✔ Got gridExtra 2.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.11 MB)
#> ✔ Got rappdirs 0.3.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (45.26 kB)
#> ✔ Got car 3.1-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.54 MB)
#> ✔ Got registry 0.5-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (194.47 kB)
#> ✔ Got knitr 1.50 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.10 MB)
#> ✔ Got xtable 1.8-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (706.31 kB)
#> ✔ Got openssl 2.3.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.31 MB)
#> ✔ Got gtable 0.3.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (222.55 kB)
#> ✔ Got mime 0.13 (x86_64-pc-linux-gnu-ubuntu-24.04) (44.52 kB)
#> ✔ Got RcppTOML 0.2.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (244.01 kB)
#> ✔ Got viridisLite 0.4.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.30 MB)
#> ✔ Got future.apply 1.20.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (191.24 kB)
#> ✔ Got shiny 1.11.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.44 MB)
#> ✔ Got rstatix 0.7.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (615.19 kB)
#> ✔ Got withr 3.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (223.90 kB)
#> ✔ Got cowplot 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.38 MB)
#> ✔ Got jsonlite 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.09 MB)
#> ✔ Got rmarkdown 2.30 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.64 MB)
#> ✔ Got igraph 2.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.84 MB)
#> ✔ Got doBy 4.7.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.91 MB)
#> ✔ Got patchwork 1.3.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.35 MB)
#> ✔ Got CellChat 2.2.0 (source) (29.82 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Get:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease [126 kB]
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Fetched 126 kB in 0s (353 kB/s)
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
#> 0 upgraded, 0 newly installed, 0 to remove and 62 not upgraded.
#> ℹ Building assorthead 1.4.0
#> ✔ Installed abind 1.4-8  (77ms)
#> ✔ Installed askpass 1.2.1  (86ms)
#> ✔ Installed backports 1.5.0  (115ms)
#> ✔ Installed base64enc 0.1-3  (64ms)
#> ✔ Installed BiocManager 1.30.27  (66ms)
#> ✔ Installed broom 1.0.10  (96ms)
#> ✔ Installed bslib 0.9.0  (141ms)
#> ✔ Installed cachem 1.1.0  (70ms)
#> ✔ Installed car 3.1-3  (70ms)
#> ✔ Installed carData 3.0-5  (70ms)
#> ✔ Installed circlize 0.4.16  (73ms)
#> ✔ Installed cli 3.6.5  (70ms)
#> ✔ Installed clue 0.3-66  (66ms)
#> ✔ Installed coda 0.19-4.1  (94ms)
#> ✔ Installed colorspace 2.1-2  (73ms)
#> ✔ Installed commonmark 2.0.0  (69ms)
#> ✔ Installed corrplot 0.95  (72ms)
#> ✔ Installed cowplot 1.2.0  (69ms)
#> ✔ Installed crayon 1.5.3  (68ms)
#> ✔ Installed crosstalk 1.2.2  (87ms)
#> ✔ Installed curl 7.0.0  (129ms)
#> ✔ Installed data.table 1.17.8  (78ms)
#> ✔ Installed Deriv 4.2.0  (98ms)
#> ✔ Installed digest 0.6.39  (99ms)
#> ✔ Installed doBy 4.7.0  (76ms)
#> ✔ Installed doParallel 1.0.17  (73ms)
#> ✔ Installed dplyr 1.1.4  (79ms)
#> ✔ Installed evaluate 1.0.5  (99ms)
#> ✔ Installed farver 2.1.2  (97ms)
#> ✔ Installed fastmap 1.2.0  (64ms)
#> ✔ Installed FNN 1.1.4.1  (64ms)
#> ✔ Installed fontawesome 0.5.3  (68ms)
#> ✔ Installed foreach 1.5.2  (66ms)
#> ✔ Installed Formula 1.2-5  (65ms)
#> ✔ Installed fs 1.6.6  (68ms)
#> ✔ Installed future 1.68.0  (102ms)
#> ✔ Installed future.apply 1.20.0  (69ms)
#> ✔ Installed generics 0.1.4  (67ms)
#> ℹ Building BiocGenerics 0.56.0
#> ✔ Installed GetoptLong 1.1.0  (82ms)
#> ✔ Installed ggalluvial 0.12.5  (1.1s)
#> ✔ Built assorthead 1.4.0 (3.5s)
#> ✔ Installed ggnetwork 0.5.14  (65ms)
#> ✔ Installed ggpubr 0.6.2  (1s)
#> ✔ Installed ggplot2 4.0.1  (1.1s)
#> ✔ Installed ggrepel 0.9.6  (81ms)
#> ✔ Installed ggsci 4.1.0  (101ms)
#> ✔ Installed ggsignif 0.6.4  (90ms)
#> ✔ Installed GlobalOptions 0.1.3  (70ms)
#> ✔ Installed globals 0.18.0  (104ms)
#> ✔ Installed glue 1.8.0  (71ms)
#> ✔ Installed gridBase 0.4-7  (90ms)
#> ✔ Installed gridExtra 2.3  (90ms)
#> ✔ Installed gtable 0.3.6  (70ms)
#> ✔ Installed here 1.0.2  (69ms)
#> ✔ Installed highr 0.11  (68ms)
#> ✔ Installed htmltools 0.5.8.1  (110ms)
#> ✔ Built BiocGenerics 0.56.0 (3.2s)
#> ✔ Installed htmlwidgets 1.6.4  (93ms)
#> ✔ Installed httpuv 1.6.16  (87ms)
#> ✔ Installed httr 1.4.7  (67ms)
#> ✔ Installed irlba 2.3.5.1  (1s)
#> ✔ Installed igraph 2.2.1  (1.1s)
#> ✔ Installed isoband 0.2.7  (101ms)
#> ✔ Installed iterators 1.0.14  (96ms)
#> ✔ Installed jquerylib 0.1.4  (62ms)
#> ✔ Installed jsonlite 2.0.0  (62ms)
#> ✔ Installed knitr 1.50  (73ms)
#> ✔ Installed labeling 0.4.3  (71ms)
#> ✔ Installed later 1.4.4  (62ms)
#> ✔ Installed lazyeval 0.2.2  (62ms)
#> ✔ Installed lifecycle 1.0.4  (97ms)
#> ✔ Installed listenv 0.10.0  (64ms)
#> ✔ Installed lme4 1.1-37  (68ms)
#> ✔ Installed magrittr 2.0.4  (69ms)
#> ✔ Installed MatrixModels 0.5-4  (68ms)
#> ✔ Installed matrixStats 1.5.0  (66ms)
#> ✔ Installed memoise 2.0.1  (62ms)
#> ✔ Installed microbenchmark 1.5.0  (103ms)
#> ✔ Installed mime 0.13  (64ms)
#> ✔ Installed minqa 1.2.8  (62ms)
#> ✔ Installed modelr 0.1.11  (62ms)
#> ✔ Installed network 1.19.0  (65ms)
#> ✔ Installed nloptr 2.2.1  (64ms)
#> ✔ Installed numDeriv 2016.8-1.1  (62ms)
#> ✔ Installed openssl 2.3.4  (99ms)
#> ✔ Installed otel 0.2.0  (67ms)
#> ✔ Installed parallelly 1.45.1  (67ms)
#> ✔ Installed patchwork 1.3.2  (72ms)
#> ✔ Installed pbapply 1.7-4  (67ms)
#> ✔ Installed pbkrtest 0.5.5  (64ms)
#> ✔ Installed pillar 1.11.1  (64ms)
#> ✔ Installed pkgconfig 2.0.3  (97ms)
#> ✔ Installed plotly 4.11.0  (100ms)
#> ✔ Installed plyr 1.8.9  (69ms)
#> ✔ Installed png 0.1-8  (62ms)
#> ✔ Installed polynom 1.4-1  (61ms)
#> ✔ Installed promises 1.5.0  (63ms)
#> ✔ Installed purrr 1.2.0  (65ms)
#> ✔ Installed quantreg 6.1  (64ms)
#> ✔ Installed R6 2.6.1  (98ms)
#> ✔ Installed rappdirs 0.3.3  (69ms)
#> ✔ Installed rbibutils 2.4  (69ms)
#> ✔ Installed RColorBrewer 1.1-3  (62ms)
#> ✔ Installed Rcpp 1.1.0  (65ms)
#> ✔ Installed RcppEigen 0.3.4.0.2  (1.1s)
#> ✔ Installed RcppTOML 0.2.3  (1.1s)
#> ✔ Installed Rdpack 2.6.4  (101ms)
#> ✔ Installed reformulas 0.4.2  (64ms)
#> ✔ Installed registry 0.5-1  (61ms)
#> ✔ Installed reshape2 1.4.5  (62ms)
#> ✔ Installed reticulate 1.44.1  (66ms)
#> ✔ Installed rjson 0.2.23  (62ms)
#> ✔ Installed rlang 1.1.6  (67ms)
#> ✔ Installed rmarkdown 2.30  (115ms)
#> ✔ Installed rngtools 1.5.2  (70ms)
#> ✔ Installed rprojroot 2.1.1  (65ms)
#> ✔ Installed RSpectra 0.16-2  (65ms)
#> ✔ Installed rstatix 0.7.3  (66ms)
#> ✔ Installed S7 0.2.1  (66ms)
#> ✔ Installed sass 0.4.10  (67ms)
#> ✔ Installed scales 1.4.0  (106ms)
#> ✔ Installed shape 1.4.6.1  (69ms)
#> ✔ Installed sna 2.8  (1s)
#> ✔ Installed shiny 1.11.1  (1.1s)
#> ✔ Installed sourcetools 0.1.7-1  (73ms)
#> ✔ Installed SparseM 1.84-2  (62ms)
#> ✔ Installed statnet.common 4.12.0  (62ms)
#> ✔ Installed stringi 1.8.7  (1.1s)
#> ✔ Installed stringr 1.6.0  (1.1s)
#> ✔ Installed svglite 2.2.2  (69ms)
#> ✔ Installed sys 3.4.3  (65ms)
#> ✔ Installed systemfonts 1.3.1  (67ms)
#> ✔ Installed textshaping 1.0.4  (68ms)
#> ✔ Installed tibble 3.3.0  (70ms)
#> ✔ Installed tidyr 1.3.1  (109ms)
#> ✔ Installed tidyselect 1.2.1  (74ms)
#> ✔ Installed tinytex 0.58  (65ms)
#> ✔ Installed utf8 1.2.6  (64ms)
#> ✔ Installed vctrs 0.6.5  (66ms)
#> ✔ Installed viridisLite 0.4.2  (65ms)
#> ✔ Installed withr 3.0.2  (64ms)
#> ✔ Installed xfun 0.54  (104ms)
#> ✔ Installed xtable 1.8-4  (104ms)
#> ✔ Installed yaml 2.3.11  (66ms)
#> ✔ Installed BiocGenerics 0.56.0  (1s)
#> ℹ Building Biobase 2.70.0
#> ✔ Installed assorthead 1.4.0  (1.1s)
#> ℹ Building BiocNeighbors 2.4.0
#> ℹ Building S4Vectors 0.48.0
#> ✔ Built Biobase 2.70.0 (9s)
#> ✔ Installed Biobase 2.70.0  (68ms)
#> ℹ Building NMF 0.28
#> ✔ Built S4Vectors 0.48.0 (22.7s)
#> ✔ Installed S4Vectors 0.48.0  (64ms)
#> ℹ Building IRanges 2.44.0
#> ✔ Built NMF 0.28 (16.7s)
#> ✔ Installed NMF 0.28  (64ms)
#> ✔ Built BiocNeighbors 2.4.0 (35.3s)
#> ✔ Installed BiocNeighbors 2.4.0  (80ms)
#> ✔ Built IRanges 2.44.0 (36.7s)
#> ✔ Installed IRanges 2.44.0  (41ms)
#> ℹ Building ComplexHeatmap 2.26.0
#> ✔ Built ComplexHeatmap 2.26.0 (12.3s)
#> ✔ Installed ComplexHeatmap 2.26.0  (1.1s)
#> ℹ Packaging CellChat 2.2.0
#> ✔ Packaged CellChat 2.2.0 (4s)
#> ℹ Building CellChat 2.2.0
#> ✔ Built CellChat 2.2.0 (23.8s)
#> ✔ Installed CellChat 2.2.0 (github::jinworks/CellChat@623f48f) (58ms)
#> ✔ 1 pkg + 160 deps: kept 9, added 151, dld 146 (NA B) [2m 27.4s]
#>  
#> → Will install 2 packages.
#> → All 2 packages (0 B) are cached.
#> + RcppArmadillo   15.2.2-1 
#> + presto          1.0.0    [bld][cmp] (GitHub: 7636b3d)
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 2 pkgs are cached
#> ✔ Got presto 1.0.0 (source) (746.05 kB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 62 not upgraded.
#> ✔ Installed RcppArmadillo 15.2.2-1  (1.1s)
#> ℹ Packaging presto 1.0.0
#> ✔ Packaged presto 1.0.0 (721ms)
#> ℹ Building presto 1.0.0
#> ✔ Built presto 1.0.0 (11.8s)
#> ✔ Installed presto 1.0.0 (github::immunogenomics/presto@7636b3d) (1s)
#> ✔ 1 pkg + 24 deps: kept 23, added 2, dld 1 (NA B) [18.9s]
#> [1] "Create a CellChat object from a data matrix"
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Ductal, Ngn3-high-EP, Endocrine, Ngn3-low-EP, Pre-endocrine 
#> The number of highly variable ligand-receptor pairs used for signaling inference is 841 
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-12-03 09:22:06.424816]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-12-03 09:23:17.234595]"

CellChatPlot(pancreas_sub, plot_type = "aggregate")

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways



CellChatPlot(pancreas_sub, plot_type = "pathway")
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
#> NULL

CellChatPlot(pancreas_sub, plot_type = "bubble")
#> Comparing communications on a single object 
#> Comparing communications on a single object 
#> Comparing communications on a single object 


CellChatPlot(pancreas_sub, plot_type = "gene")


CellChatPlot(pancreas_sub, plot_type = "heatmap")

#> NULL
```
