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
#> → Will install 173 packages.
#> → Will download 21 CRAN packages (17.65 MB), cached: 152 (0 B).
#> + BH                     1.90.0-1  
#> + Biobase                2.70.0    [bld][cmp]
#> + BiocGenerics           0.56.0    [bld]
#> + BiocManager            1.30.27   
#> + BiocParallel           1.44.0    [bld][cmp]
#> + BiocStyle              2.38.0    [bld]
#> + DESeq2                 1.50.2    [bld][cmp]
#> + DelayedArray           0.36.0    [bld][cmp]
#> + FNN                    1.1.4.1   
#> + GenomicRanges          1.62.1    [bld][cmp]
#> + IRanges                2.44.0    [bld][cmp]
#> + MatrixGenerics         1.22.0    [bld]
#> + R6                     2.6.1     
#> + RANN                   2.6.2     
#> + RColorBrewer           1.1-3     
#> + ROCR                   1.0-11    
#> + RSpectra               0.16-2    [bld][cmp][dl] (128.03 kB)
#> + Rcpp                   1.1.0     [bld][cmp][dl] (3.11 MB)
#> + RcppAnnoy              0.0.22    [bld][cmp][dl] (156.28 kB)
#> + RcppArmadillo          15.2.2-1  [bld][cmp][dl] (2.12 MB)
#> + RcppEigen              0.3.4.0.2 [bld][cmp][dl] (1.76 MB)
#> + RcppHNSW               0.6.0     [bld][cmp][dl] (37.08 kB)
#> + RcppProgress           0.4.2     
#> + RcppTOML               0.2.3     [bld][cmp][dl] (139.89 kB)
#> + Rtsne                  0.17      [bld][cmp][dl] (70.90 kB)
#> + S4Arrays               1.10.1    [bld][cmp]
#> + S4Vectors              0.48.0    [bld][cmp]
#> + S7                     0.2.1     
#> + Seqinfo                1.0.0     [bld]
#> + Seurat                 5.4.0     [bld][cmp][dl] (2.11 MB)
#> + SeuratObject           5.3.0     [bld][cmp][dl] (767.77 kB)
#> + SingleCellExperiment   1.32.0    [bld]
#> + SparseArray            1.10.7    [bld][cmp]
#> + SummarizedExperiment   1.40.0    [bld]
#> + XVector                0.50.0    [bld][cmp]
#> + abind                  1.4-8     
#> + askpass                1.2.1     
#> + backports              1.5.0     
#> + base64enc              0.1-3     
#> + bitops                 1.0-9     
#> + bookdown               0.46       + ✔ pandoc
#> + brio                   1.1.5     
#> + broom                  1.0.11    
#> + bslib                  0.9.0     
#> + caTools                1.18.3    
#> + cachem                 1.1.0     
#> + callr                  3.7.6     
#> + cli                    3.6.5     
#> + commonmark             2.0.0     
#> + cowplot                1.2.0     
#> + cpp11                  0.5.2     
#> + crayon                 1.5.3     
#> + crosstalk              1.2.2     
#> + curl                   7.0.0      + ✔ libcurl4-openssl-dev, ✔ libssl-dev
#> + data.table             1.17.8    
#> + deldir                 2.0-4     
#> + desc                   1.4.3     
#> + diffobj                0.3.6     
#> + digest                 0.6.39    
#> + dotCall64              1.2       
#> + dplyr                  1.1.4     
#> + dqrng                  0.4.1     [bld][cmp][dl] (267.02 kB)
#> + evaluate               1.0.5     
#> + farver                 2.1.2     
#> + fastDummies            1.7.5     
#> + fastmap                1.2.0     
#> + fitdistrplus           1.2-4     
#> + fontawesome            0.5.3     
#> + formatR                1.14      
#> + fs                     1.6.6      + ✔ make
#> + futile.logger          1.4.3     
#> + futile.options         1.0.1     
#> + future                 1.68.0    
#> + future.apply           1.20.1    
#> + generics               0.1.4     
#> + ggplot2                4.0.1     
#> + ggrepel                0.9.6     [bld][cmp][dl] (149.97 kB)
#> + ggridges               0.5.7     
#> + globals                0.18.0    
#> + glue                   1.8.0     
#> + goftest                1.2-3     
#> + gplots                 3.3.0     
#> + gridExtra              2.3       
#> + gtable                 0.3.6     
#> + gtools                 3.9.5     
#> + here                   1.0.2     
#> + highr                  0.11      
#> + htmltools              0.5.9     
#> + htmlwidgets            1.6.4     
#> + httpuv                 1.6.16    [bld][cmp][dl] (1.89 MB) + ✔ make, ✔ zlib1g-dev
#> + httr                   1.4.7     
#> + ica                    1.0-3     
#> + igraph                 2.2.1      + ✔ libglpk-dev, ✔ libxml2-dev
#> + irlba                  2.3.5.1   
#> + isoband                0.3.0     
#> + jquerylib              0.1.4     
#> + jsonlite               2.0.0     
#> + knitr                  1.50       + ✔ pandoc
#> + labeling               0.4.3     
#> + lambda.r               1.2.4     
#> + later                  1.4.4     [bld][cmp][dl] (70.68 kB)
#> + lazyeval               0.2.2     
#> + lifecycle              1.0.4     
#> + listenv                0.10.0    
#> + lmtest                 0.9-40    
#> + locfit                 1.5-9.12  
#> + magrittr               2.0.4     
#> + matrixStats            1.5.0     
#> + memoise                2.0.1     
#> + mime                   0.13      
#> + miniUI                 0.1.2     
#> + openssl                2.3.4      + ✔ libssl-dev
#> + otel                   0.2.0     
#> + parallelly             1.46.0    
#> + patchwork              1.3.2     
#> + pbapply                1.7-4     
#> + pillar                 1.11.1    
#> + pkgbuild               1.4.8     
#> + pkgconfig              2.0.3     
#> + pkgload                1.4.1     
#> + plotly                 4.11.0    
#> + plyr                   1.8.9     [bld][cmp][dl] (401.49 kB)
#> + png                    0.1-8      + ✔ libpng-dev
#> + polyclip               1.10-7    
#> + praise                 1.0.0     
#> + presto                 1.0.0     [bld][cmp] (GitHub: 7636b3d)
#> + processx               3.8.6     
#> + progressr              0.18.0    
#> + promises               1.5.0     
#> + ps                     1.9.1     
#> + purrr                  1.2.0     
#> + rappdirs               0.3.3     
#> + reshape2               1.4.5     [bld][cmp][dl] (38.08 kB)
#> + reticulate             1.44.1    [bld][cmp][dl] (1.69 MB) + ✔ python3
#> + rlang                  1.1.6     
#> + rmarkdown              2.30       + ✔ pandoc
#> + rprojroot              2.1.1     
#> + sass                   0.4.10     + ✔ make
#> + scales                 1.4.0     
#> + scattermore            1.2       
#> + sctransform            0.4.2     [bld][cmp][dl] (186.27 kB)
#> + shiny                  1.12.1    
#> + sitmo                  2.0.2     [bld][cmp][dl] (132.96 kB)
#> + snow                   0.4-4     
#> + sourcetools            0.1.7-1   
#> + sp                     2.2-0     
#> + spam                   2.11-1    [bld][cmp][dl] (1.81 MB)
#> + spatstat.data          3.1-9     
#> + spatstat.explore       3.6-0     
#> + spatstat.geom          3.6-1     
#> + spatstat.random        3.4-3     
#> + spatstat.sparse        3.1-0     
#> + spatstat.univar        3.1-5     
#> + spatstat.utils         3.2-0     
#> + stringi                1.8.7      + ✔ libicu-dev
#> + stringr                1.6.0     
#> + sys                    3.4.3     
#> + tensor                 1.5.1     
#> + testthat               3.3.1     
#> + tibble                 3.3.0     
#> + tidyr                  1.3.1     
#> + tidyselect             1.2.1     
#> + tinytex                0.58      
#> + utf8                   1.2.6     
#> + uwot                   0.2.4     [bld][cmp][dl] (614.67 kB)
#> + vctrs                  0.6.5     
#> + viridisLite            0.4.2     
#> + waldo                  0.6.2     
#> + withr                  3.0.2     
#> + xfun                   0.54      
#> + xtable                 1.8-4     
#> + yaml                   2.3.12    
#> + zoo                    1.8-14    
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 21 pkgs (17.65 MB), 152 cached
#> ✔ Cached copy of Seurat 5.4.0 (source) is the latest build
#> ✔ Cached copy of SeuratObject 5.3.0 (source) is the latest build
#> ✔ Cached copy of sitmo 2.0.2 (source) is the latest build
#> ✔ Got BiocGenerics 0.56.0 (source) (61.53 kB)
#> ✔ Got askpass 1.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (21.89 kB)
#> ✔ Got DelayedArray 0.36.0 (source) (816.30 kB)
#> ✔ Got base64enc 0.1-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (26.57 kB)
#> ✔ Got reshape2 1.4.5 (source) (39.30 kB)
#> ✔ Got BiocStyle 2.38.0 (source) (908.65 kB)
#> ✔ Got abind 1.4-8 (x86_64-pc-linux-gnu-ubuntu-24.04) (64.92 kB)
#> ✔ Got MatrixGenerics 1.22.0 (source) (34.60 kB)
#> ✔ Got backports 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (115.90 kB)
#> ✔ Got BiocParallel 1.44.0 (source) (1.11 MB)
#> ✔ Got rprojroot 2.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (113.23 kB)
#> ✔ Got bitops 1.0-9 (x86_64-pc-linux-gnu-ubuntu-24.04) (26.02 kB)
#> ✔ Got scattermore 1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (339.10 kB)
#> ✔ Got IRanges 2.44.0 (source) (496.06 kB)
#> ✔ Got dqrng 0.4.1 (source) (269.33 kB)
#> ✔ Got Biobase 2.70.0 (source) (1.98 MB)
#> ✔ Got polyclip 1.10-7 (x86_64-pc-linux-gnu-ubuntu-24.04) (120.02 kB)
#> ✔ Got Seqinfo 1.0.0 (source) (254.66 kB)
#> ✔ Got BiocManager 1.30.27 (x86_64-pc-linux-gnu-ubuntu-24.04) (666.57 kB)
#> ✔ Got SummarizedExperiment 1.40.0 (source) (690.87 kB)
#> ✔ Got DESeq2 1.50.2 (source) (2.69 MB)
#> ✔ Got ps 1.9.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (499.64 kB)
#> ✔ Got fastmap 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (66.05 kB)
#> ✔ Got futile.logger 1.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (96.77 kB)
#> ✔ Got FNN 1.1.4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (128.22 kB)
#> ✔ Got S4Vectors 0.48.0 (source) (843.90 kB)
#> ✔ Got SingleCellExperiment 1.32.0 (source) (987.98 kB)
#> ✔ Got xtable 1.8-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (706.31 kB)
#> ✔ Got ggridges 0.5.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.25 MB)
#> ✔ Got future 1.68.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (966.66 kB)
#> ✔ Got gtable 0.3.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (222.55 kB)
#> ✔ Got highr 0.11 (x86_64-pc-linux-gnu-ubuntu-24.04) (37.50 kB)
#> ✔ Got htmlwidgets 1.6.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (815.05 kB)
#> ✔ Got labeling 0.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (60.95 kB)
#> ✔ Got lambda.r 1.2.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (110.01 kB)
#> ✔ Got sp 2.2-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.31 MB)
#> ✔ Got listenv 0.10.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (106.64 kB)
#> ✔ Got bslib 0.9.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.67 MB)
#> ✔ Got lazyeval 0.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (157.44 kB)
#> ✔ Got later 1.4.4 (source) (70.93 kB)
#> ✔ Got matrixStats 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (474.75 kB)
#> ✔ Got pkgload 1.4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (222.06 kB)
#> ✔ Got pkgbuild 1.4.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (208.62 kB)
#> ✔ Got presto 1.0.0 (source) (746.05 kB)
#> ✔ Got pbapply 1.7-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (100.53 kB)
#> ✔ Got openssl 2.3.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.31 MB)
#> ✔ Got purrr 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (553.55 kB)
#> ✔ Got png 0.1-8 (x86_64-pc-linux-gnu-ubuntu-24.04) (40.57 kB)
#> ✔ Got RcppProgress 0.4.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (20.19 kB)
#> ✔ Got tensor 1.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (13.96 kB)
#> ✔ Got ROCR 1.0-11 (x86_64-pc-linux-gnu-ubuntu-24.04) (465.58 kB)
#> ✔ Got RcppTOML 0.2.3 (source) (140.65 kB)
#> ✔ Got gplots 3.3.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.90 MB)
#> ✔ Got spatstat.univar 3.1-5 (x86_64-pc-linux-gnu-ubuntu-24.04) (325.81 kB)
#> ✔ Got RcppArmadillo 15.2.2-1 (source) (2.12 MB)
#> ✔ Got tibble 3.3.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (680.35 kB)
#> ✔ Got utf8 1.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (151.81 kB)
#> ✔ Got ggplot2 4.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (8.47 MB)
#> ✔ Got sctransform 0.4.2 (source) (186.38 kB)
#> ✔ Got igraph 2.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.84 MB)
#> ✔ Got XVector 0.50.0 (source) (71.42 kB)
#> ✔ Got withr 3.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (223.90 kB)
#> ✔ Got stringi 1.8.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.29 MB)
#> ✔ Got deldir 2.0-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (278.91 kB)
#> ✔ Got cpp11 0.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (289.40 kB)
#> ✔ Got crosstalk 1.2.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (412.02 kB)
#> ✔ Got caTools 1.18.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (205.11 kB)
#> ✔ Got digest 0.6.39 (x86_64-pc-linux-gnu-ubuntu-24.04) (230.38 kB)
#> ✔ Got vctrs 0.6.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.31 MB)
#> ✔ Got callr 3.7.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (449.24 kB)
#> ✔ Got evaluate 1.0.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (102.86 kB)
#> ✔ Got formatR 1.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (151.65 kB)
#> ✔ Got generics 0.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (80.38 kB)
#> ✔ Got fs 1.6.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (310.07 kB)
#> ✔ Got bookdown 0.46 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.09 MB)
#> ✔ Got ica 1.0-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (85.54 kB)
#> ✔ Got globals 0.18.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (158.33 kB)
#> ✔ Got testthat 3.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.11 MB)
#> ✔ Got irlba 2.3.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (292.19 kB)
#> ✔ Got farver 2.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.47 MB)
#> ✔ Got jsonlite 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.09 MB)
#> ✔ Got gridExtra 2.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.11 MB)
#> ✔ Got jquerylib 0.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (526.85 kB)
#> ✔ Got lmtest 0.9-40 (x86_64-pc-linux-gnu-ubuntu-24.04) (403.49 kB)
#> ✔ Got memoise 2.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (48.86 kB)
#> ✔ Got miniUI 0.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (35.20 kB)
#> ✔ Got httpuv 1.6.16 (source) (1.90 MB)
#> ✔ Got otel 0.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (279.62 kB)
#> ✔ Got pillar 1.11.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (660.46 kB)
#> ✔ Got fitdistrplus 1.2-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.96 MB)
#> ✔ Got RcppAnnoy 0.0.22 (source) (157.62 kB)
#> ✔ Got RColorBrewer 1.1-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (51.81 kB)
#> ✔ Got rlang 1.1.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.59 MB)
#> ✔ Got RcppHNSW 0.6.0 (source) (37.28 kB)
#> ✔ Got RcppEigen 0.3.4.0.2 (source) (1.76 MB)
#> ✔ Got spatstat.sparse 3.1-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (213.31 kB)
#> ✔ Got Rcpp 1.1.0 (source) (3.11 MB)
#> ✔ Got rmarkdown 2.30 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.64 MB)
#> ✔ Got tinytex 0.58 (x86_64-pc-linux-gnu-ubuntu-24.04) (143.77 kB)
#> ✔ Got sass 0.4.10 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.46 MB)
#> ✔ Got sys 3.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (40.73 kB)
#> ✔ Got plotly 4.11.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.86 MB)
#> ✔ Got stringr 1.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (333.76 kB)
#> ✔ Got waldo 0.6.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (135.16 kB)
#> ✔ Got xfun 0.54 (x86_64-pc-linux-gnu-ubuntu-24.04) (583.00 kB)
#> ✔ Got spatstat.random 3.4-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.23 MB)
#> ✔ Got uwot 0.2.4 (source) (609.90 kB)
#> ✔ Got commonmark 2.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (148.07 kB)
#> ✔ Got tidyr 1.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.18 MB)
#> ✔ Got zoo 1.8-14 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.03 MB)
#> ✔ Got curl 7.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (788.30 kB)
#> ✔ Got spatstat.data 3.1-9 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.17 MB)
#> ✔ Got desc 1.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (334.52 kB)
#> ✔ Got futile.options 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (18.42 kB)
#> ✔ Got cowplot 1.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.38 MB)
#> ✔ Got glue 1.8.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (168.12 kB)
#> ✔ Got goftest 1.2-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (58.64 kB)
#> ✔ Got diffobj 0.3.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.00 MB)
#> ✔ Got spatstat.geom 3.6-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.18 MB)
#> ✔ Got magrittr 2.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (221.44 kB)
#> ✔ Got RANN 2.6.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (43.84 kB)
#> ✔ Got httr 1.4.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (486.52 kB)
#> ✔ Got plyr 1.8.9 (source) (402.04 kB)
#> ✔ Got knitr 1.50 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.10 MB)
#> ✔ Got processx 3.8.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (337.39 kB)
#> ✔ Got progressr 0.18.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (412.00 kB)
#> ✔ Got S7 0.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (325.15 kB)
#> ✔ Got snow 0.4-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (97.07 kB)
#> ✔ Got tidyselect 1.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (225.28 kB)
#> ✔ Got fontawesome 0.5.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.40 MB)
#> ✔ Got brio 1.1.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (35.30 kB)
#> ✔ Got crayon 1.5.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (163.30 kB)
#> ✔ Got cli 3.6.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.34 MB)
#> ✔ Got spam 2.11-1 (source) (1.81 MB)
#> ✔ Got fastDummies 1.7.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (41.76 kB)
#> ✔ Got gtools 3.9.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (361.11 kB)
#> ✔ Got reticulate 1.44.1 (source) (1.69 MB)
#> ✔ Got R6 2.6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (86.81 kB)
#> ✔ Got patchwork 1.3.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.35 MB)
#> ✔ Got spatstat.explore 3.6-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.56 MB)
#> ✔ Got mime 0.13 (x86_64-pc-linux-gnu-ubuntu-24.04) (44.52 kB)
#> ✔ Got rappdirs 0.3.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (45.26 kB)
#> ✔ Got sourcetools 0.1.7-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (46.59 kB)
#> ✔ Got Rtsne 0.17 (source) (72.50 kB)
#> ✔ Got data.table 1.17.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.67 MB)
#> ✔ Got pkgconfig 2.0.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (18.08 kB)
#> ✔ Got cachem 1.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (67.49 kB)
#> ✔ Got lifecycle 1.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (125.07 kB)
#> ✔ Got locfit 1.5-9.12 (x86_64-pc-linux-gnu-ubuntu-24.04) (539.82 kB)
#> ✔ Got RSpectra 0.16-2 (source) (127.62 kB)
#> ✔ Got ggrepel 0.9.6 (source) (150.12 kB)
#> ✔ Got dotCall64 1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (31.26 kB)
#> ✔ Got here 1.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (53.59 kB)
#> ✔ Got praise 1.0.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (16.39 kB)
#> ✔ Got spatstat.utils 3.2-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (401.56 kB)
#> ✔ Got scales 1.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (841.36 kB)
#> ✔ Got viridisLite 0.4.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.30 MB)
#> ✔ Got dplyr 1.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.49 MB)
#> ✔ Got promises 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.67 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install pandoc libcurl4-openssl-dev libssl-dev make zlib1g-dev libglpk-dev libxml2-dev libpng-dev python3 libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> pandoc is already the newest version (3.1.3+ds-2).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libpng-dev is already the newest version (1.6.43-5ubuntu0.1).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 49 not upgraded.
#> ℹ Building Rcpp 1.1.0
#> ✔ Installed abind 1.4-8  (80ms)
#> ✔ Installed askpass 1.2.1  (91ms)
#> ✔ Installed backports 1.5.0  (126ms)
#> ✔ Installed base64enc 0.1-3  (76ms)
#> ✔ Installed BiocManager 1.30.27  (40ms)
#> ✔ Installed bitops 1.0-9  (30ms)
#> ✔ Installed bookdown 0.46  (56ms)
#> ✔ Installed brio 1.1.5  (29ms)
#> ✔ Installed broom 1.0.11  (83ms)
#> ✔ Installed cachem 1.1.0  (30ms)
#> ✔ Installed bslib 0.9.0  (232ms)
#> ✔ Installed BH 1.90.0-1  (932ms)
#> ✔ Installed callr 3.7.6  (228ms)
#> ✔ Installed caTools 1.18.3  (114ms)
#> ✔ Installed cli 3.6.5  (75ms)
#> ✔ Installed commonmark 2.0.0  (75ms)
#> ✔ Installed cowplot 1.2.0  (77ms)
#> ✔ Installed cpp11 0.5.2  (136ms)
#> ✔ Installed crayon 1.5.3  (101ms)
#> ✔ Installed crosstalk 1.2.2  (88ms)
#> ✔ Installed curl 7.0.0  (85ms)
#> ✔ Installed data.table 1.17.8  (87ms)
#> ✔ Installed deldir 2.0-4  (81ms)
#> ✔ Installed desc 1.4.3  (125ms)
#> ✔ Installed diffobj 0.3.6  (134ms)
#> ✔ Installed digest 0.6.39  (80ms)
#> ✔ Installed dotCall64 1.2  (80ms)
#> ✔ Installed dplyr 1.1.4  (78ms)
#> ✔ Installed evaluate 1.0.5  (78ms)
#> ✔ Installed farver 2.1.2  (127ms)
#> ✔ Installed fastDummies 1.7.5  (137ms)
#> ✔ Installed fastmap 1.2.0  (75ms)
#> ✔ Installed fitdistrplus 1.2-4  (78ms)
#> ✔ Installed FNN 1.1.4.1  (75ms)
#> ✔ Installed fontawesome 0.5.3  (77ms)
#> ✔ Installed formatR 1.14  (75ms)
#> ✔ Installed fs 1.6.6  (132ms)
#> ✔ Installed futile.logger 1.4.3  (77ms)
#> ✔ Installed futile.options 1.0.1  (74ms)
#> ✔ Installed future 1.68.0  (78ms)
#> ✔ Installed future.apply 1.20.1  (79ms)
#> ✔ Installed generics 0.1.4  (76ms)
#> ℹ Building BiocGenerics 0.56.0
#> ✔ Installed ggplot2 4.0.1  (158ms)
#> ✔ Installed ggridges 0.5.7  (64ms)
#> ✔ Installed globals 0.18.0  (30ms)
#> ✔ Installed glue 1.8.0  (36ms)
#> ✔ Installed goftest 1.2-3  (66ms)
#> ✔ Installed gplots 3.3.0  (228ms)
#> ✔ Installed gridExtra 2.3  (178ms)
#> ✔ Installed gtable 0.3.6  (96ms)
#> ✔ Installed gtools 3.9.5  (79ms)
#> ✔ Installed here 1.0.2  (76ms)
#> ✔ Installed highr 0.11  (77ms)
#> ✔ Installed htmltools 0.5.9  (132ms)
#> ✔ Installed htmlwidgets 1.6.4  (141ms)
#> ✔ Installed httr 1.4.7  (79ms)
#> ✔ Installed ica 1.0-3  (75ms)
#> ✔ Installed irlba 2.3.5.1  (28ms)
#> ✔ Installed igraph 2.2.1  (169ms)
#> ✔ Installed isoband 0.3.0  (139ms)
#> ✔ Installed jquerylib 0.1.4  (141ms)
#> ✔ Installed jsonlite 2.0.0  (86ms)
#> ✔ Installed knitr 1.50  (109ms)
#> ✔ Installed labeling 0.4.3  (125ms)
#> ✔ Installed lambda.r 1.2.4  (102ms)
#> ✔ Installed lazyeval 0.2.2  (154ms)
#> ✔ Installed lifecycle 1.0.4  (186ms)
#> ✔ Installed listenv 0.10.0  (111ms)
#> ✔ Installed lmtest 0.9-40  (106ms)
#> ✔ Installed locfit 1.5-9.12  (104ms)
#> ✔ Installed magrittr 2.0.4  (108ms)
#> ✔ Installed matrixStats 1.5.0  (107ms)
#> ℹ Building MatrixGenerics 1.22.0
#> ✔ Installed memoise 2.0.1  (191ms)
#> ✔ Installed mime 0.13  (49ms)
#> ✔ Installed miniUI 0.1.2  (40ms)
#> ✔ Installed openssl 2.3.4  (62ms)
#> ✔ Installed otel 0.2.0  (1s)
#> ✔ Built BiocGenerics 0.56.0 (4.2s)
#> ✔ Installed parallelly 1.46.0  (132ms)
#> ✔ Installed patchwork 1.3.2  (92ms)
#> ✔ Installed pbapply 1.7-4  (80ms)
#> ✔ Installed pillar 1.11.1  (82ms)
#> ✔ Installed pkgbuild 1.4.8  (77ms)
#> ✔ Installed pkgconfig 2.0.3  (75ms)
#> ✔ Installed pkgload 1.4.1  (164ms)
#> ✔ Installed plotly 4.11.0  (108ms)
#> ✔ Installed png 0.1-8  (77ms)
#> ✔ Installed polyclip 1.10-7  (73ms)
#> ✔ Installed praise 1.0.0  (75ms)
#> ✔ Installed processx 3.8.6  (79ms)
#> ✔ Installed progressr 0.18.0  (155ms)
#> ✔ Installed promises 1.5.0  (103ms)
#> ✔ Installed ps 1.9.1  (83ms)
#> ✔ Installed purrr 1.2.0  (80ms)
#> ✔ Installed R6 2.6.1  (77ms)
#> ✔ Installed RANN 2.6.2  (78ms)
#> ✔ Installed rappdirs 0.3.3  (143ms)
#> ✔ Installed RColorBrewer 1.1-3  (95ms)
#> ✔ Installed RcppProgress 0.4.2  (79ms)
#> ✔ Installed rlang 1.1.6  (91ms)
#> ✔ Installed rmarkdown 2.30  (1.1s)
#> ✔ Built MatrixGenerics 1.22.0 (4.2s)
#> ✔ Installed ROCR 1.0-11  (1.1s)
#> ✔ Installed rprojroot 2.1.1  (99ms)
#> ✔ Installed S7 0.2.1  (78ms)
#> ✔ Installed sass 0.4.10  (83ms)
#> ✔ Installed scales 1.4.0  (115ms)
#> ✔ Installed scattermore 1.2  (135ms)
#> ✔ Installed SeuratObject 5.3.0  (63ms)
#> ✔ Installed Seurat 5.4.0  (189ms)
#> ✔ Installed sitmo 2.0.2  (37ms)
#> ✔ Installed shiny 1.12.1  (169ms)
#> ✔ Installed snow 0.4-4  (78ms)
#> ℹ Building BiocParallel 1.44.0
#> ✔ Installed sourcetools 0.1.7-1  (142ms)
#> ✔ Installed sp 2.2-0  (153ms)
#> ✔ Installed spatstat.data 3.1-9  (120ms)
#> ✔ Installed spatstat.explore 3.6-0  (122ms)
#> ✔ Installed spatstat.geom 3.6-1  (115ms)
#> ✔ Installed spatstat.random 3.4-3  (170ms)
#> ✔ Installed spatstat.sparse 3.1-0  (160ms)
#> ✔ Installed spatstat.univar 3.1-5  (91ms)
#> ✔ Installed spatstat.utils 3.2-0  (90ms)
#> ✔ Installed stringi 1.8.7  (105ms)
#> ✔ Installed stringr 1.6.0  (98ms)
#> ✔ Installed sys 3.4.3  (154ms)
#> ✔ Installed tensor 1.5.1  (164ms)
#> ✔ Installed testthat 3.3.1  (92ms)
#> ✔ Installed tibble 3.3.0  (96ms)
#> ✔ Installed tidyr 1.3.1  (88ms)
#> ✔ Installed tidyselect 1.2.1  (93ms)
#> ✔ Installed tinytex 0.58  (151ms)
#> ✔ Installed utf8 1.2.6  (113ms)
#> ✔ Installed vctrs 0.6.5  (101ms)
#> ✔ Installed viridisLite 0.4.2  (95ms)
#> ✔ Installed waldo 0.6.2  (78ms)
#> ✔ Installed withr 3.0.2  (82ms)
#> ✔ Installed xfun 0.54  (161ms)
#> ✔ Installed xtable 1.8-4  (119ms)
#> ✔ Installed yaml 2.3.12  (87ms)
#> ℹ Building BiocStyle 2.38.0
#> ✔ Installed zoo 1.8-14  (99ms)
#> ✔ Installed BiocGenerics 0.56.0  (38ms)
#> ℹ Building Biobase 2.70.0
#> ✔ Built BiocStyle 2.38.0 (2.9s)
#> ℹ Building S4Vectors 0.48.0
#> ✔ Built Biobase 2.70.0 (11.2s)
#> ✔ Installed Biobase 2.70.0  (82ms)
#> ✔ Installed BiocStyle 2.38.0  (50ms)
#> ✔ Installed MatrixGenerics 1.22.0  (39ms)
#> ✔ Built BiocParallel 1.44.0 (15.3s)
#> ✔ Installed BiocParallel 1.44.0  (64ms)
#> ✔ Built Rcpp 1.1.0 (27.1s)
#> ✔ Installed Rcpp 1.1.0  (116ms)
#> ℹ Building dqrng 0.4.1
#> ℹ Building ggrepel 0.9.6
#> ℹ Building later 1.4.4
#> ✔ Built S4Vectors 0.48.0 (24.4s)
#> ℹ Building plyr 1.8.9
#> ✔ Built ggrepel 0.9.6 (14.3s)
#> ℹ Building RcppAnnoy 0.0.22
#> ✔ Built dqrng 0.4.1 (15.1s)
#> ℹ Building RcppArmadillo 15.2.2-1
#> ✔ Built later 1.4.4 (21.4s)
#> ℹ Building RcppEigen 0.3.4.0.2
#> ✔ Built plyr 1.8.9 (9.8s)
#> ℹ Building RcppHNSW 0.6.0
#> ✔ Built RcppArmadillo 15.2.2-1 (27.7s)
#> ℹ Building RcppTOML 0.2.3
#> ✔ Built RcppAnnoy 0.0.22 (28.3s)
#> ℹ Building Rtsne 0.17
#> ✔ Built RcppHNSW 0.6.0 (23.3s)
#> ℹ Building spam 2.11-1
#> ✔ Built RcppTOML 0.2.3 (13s)
#> ✔ Installed dqrng 0.4.1  (66ms)
#> ✔ Installed ggrepel 0.9.6  (72ms)
#> ✔ Installed later 1.4.4  (70ms)
#> ℹ Building httpuv 1.6.16
#> ✔ Built Rtsne 0.17 (19.9s)
#> ✔ Installed plyr 1.8.9  (45ms)
#> ℹ Building reshape2 1.4.5
#> ✔ Built reshape2 1.4.5 (9.6s)
#> ✔ Installed RcppAnnoy 0.0.22  (99ms)
#> ✔ Installed RcppArmadillo 15.2.2-1  (1.2s)
#> ℹ Packaging presto 1.0.0
#> ✔ Built RcppEigen 0.3.4.0.2 (52.7s)
#> ✔ Installed RcppEigen 0.3.4.0.2  (387ms)
#> ℹ Building RSpectra 0.16-2
#> ✔ Packaged presto 1.0.0 (1.1s)
#> ℹ Building presto 1.0.0
#> ✔ Built spam 2.11-1 (30.9s)
#> ✔ Installed RcppHNSW 0.6.0  (118ms)
#> ✔ Installed RcppTOML 0.2.3  (1.1s)
#> ℹ Building reticulate 1.44.1
#> ✔ Built presto 1.0.0 (19.7s)
#> ✔ Installed presto 1.0.0 (github::immunogenomics/presto@7636b3d) (62ms)
#> ✔ Installed reshape2 1.4.5  (45ms)
#> ℹ Building sctransform 0.4.2
#> ✔ Built reticulate 1.44.1 (30.1s)
#> ✔ Installed reticulate 1.44.1  (95ms)
#> ✔ Installed Rtsne 0.17  (1.1s)
#> ✔ Installed spam 2.11-1  (78ms)
#> ✔ Installed S4Vectors 0.48.0  (1.1s)
#> ℹ Building IRanges 2.44.0
#> ✔ Built sctransform 0.4.2 (27.6s)
#> ✔ Installed sctransform 0.4.2  (70ms)
#> ✔ Built IRanges 2.44.0 (48.7s)
#> ✔ Installed IRanges 2.44.0  (92ms)
#> ℹ Building S4Arrays 1.10.1
#> ℹ Building Seqinfo 1.0.0
#> ✔ Built httpuv 1.6.16 (1m 47.6s)
#> ℹ Building XVector 0.50.0
#> ✔ Built Seqinfo 1.0.0 (7.8s)
#> ✔ Installed httpuv 1.6.16  (161ms)
#> ✔ Installed Seqinfo 1.0.0  (1s)
#> ℹ Building GenomicRanges 1.62.1
#> ✔ Built S4Arrays 1.10.1 (14.4s)
#> ✔ Installed S4Arrays 1.10.1  (49ms)
#> ✔ Built XVector 0.50.0 (12.9s)
#> ✔ Installed XVector 0.50.0  (1s)
#> ℹ Building SparseArray 1.10.7
#> ✔ Built GenomicRanges 1.62.1 (16.9s)
#> ✔ Installed GenomicRanges 1.62.1  (61ms)
#> ✔ Built RSpectra 0.16-2 (2m 4.1s)
#> ✔ Installed RSpectra 0.16-2  (402ms)
#> ℹ Building uwot 0.2.4
#> ✔ Built SparseArray 1.10.7 (23.2s)
#> ✔ Installed SparseArray 1.10.7  (53ms)
#> ℹ Building DelayedArray 0.36.0
#> ✔ Built DelayedArray 0.36.0 (15.7s)
#> ✔ Installed DelayedArray 0.36.0  (45ms)
#> ℹ Building SummarizedExperiment 1.40.0
#> ✔ Built SummarizedExperiment 1.40.0 (14.2s)
#> ✔ Installed SummarizedExperiment 1.40.0  (46ms)
#> ℹ Building DESeq2 1.50.2
#> ℹ Building SingleCellExperiment 1.32.0
#> ✔ Built uwot 0.2.4 (47.8s)
#> ✔ Installed uwot 0.2.4  (152ms)
#> ✔ Built SingleCellExperiment 1.32.0 (18.1s)
#> ✔ Installed SingleCellExperiment 1.32.0  (50ms)
#> ✔ Built DESeq2 1.50.2 (32s)
#> ✔ Installed DESeq2 1.50.2  (79ms)
#> ✔ 1 pkg + 180 deps: kept 8, added 173, dld 159 (NA B) [5m 12.2s]
#> Error in loadNamespace(x): there is no package called ‘CellChat’

CellChatPlot(pancreas_sub, plot_type = "aggregate")
#> Error in CellChatPlot(pancreas_sub, plot_type = "aggregate"): No CellChat results found in the Seurat object. Please run
#> `RunCellChat()` first.

CellChatPlot(pancreas_sub, plot_type = "pathway")
#> Error in CellChatPlot(pancreas_sub, plot_type = "pathway"): No CellChat results found in the Seurat object. Please run
#> `RunCellChat()` first.

CellChatPlot(pancreas_sub, plot_type = "bubble")
#> Error in CellChatPlot(pancreas_sub, plot_type = "bubble"): No CellChat results found in the Seurat object. Please run
#> `RunCellChat()` first.

CellChatPlot(pancreas_sub, plot_type = "gene")
#> Error in CellChatPlot(pancreas_sub, plot_type = "gene"): No CellChat results found in the Seurat object. Please run
#> `RunCellChat()` first.

CellChatPlot(pancreas_sub, plot_type = "heatmap")
#> Error in CellChatPlot(pancreas_sub, plot_type = "heatmap"): No CellChat results found in the Seurat object. Please run
#> `RunCellChat()` first.
```
