# Statistical plot of cells

Statistical plot of cells

## Usage

``` r
CellStatPlot(
  srt,
  stat.by,
  group.by = NULL,
  split.by = NULL,
  bg.by = NULL,
  cells = NULL,
  flip = FALSE,
  NA_color = "grey",
  NA_stat = TRUE,
  keep_empty = FALSE,
  individual = FALSE,
  stat_level = NULL,
  plot_type = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord",
    "venn", "upset"),
  stat_type = c("percent", "count"),
  position = c("stack", "dodge"),
  palette = "Paired",
  palcolor = NULL,
  alpha = 1,
  bg_palette = "Paired",
  bg_palcolor = NULL,
  bg_alpha = 0.2,
  label = FALSE,
  label.size = 3.5,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11
)
```

## Arguments

- srt:

  A `Seurat` object.

- stat.by:

  The column name(s) in `meta.data` specifying the variable(s) to be
  plotted.

- group.by:

  The column name in `meta.data` specifying the grouping variable.

- split.by:

  The column name in `meta.data` specifying the splitting variable.

- bg.by:

  The column name in `meta.data` specifying the background variable for
  bar plots.

- cells:

  A character vector specifying the cells to include in the plot.
  Default is `NULL`.

- flip:

  Whether to flip the plot. Default is `FALSE`.

- NA_color:

  The color to use for missing values.

- NA_stat:

  Whether to include missing values in the plot. Default is `TRUE`.

- keep_empty:

  Whether to keep empty groups in the plot. Default is `FALSE`.

- individual:

  Whether to plot individual groups separately. Default is `FALSE`.

- stat_level:

  The level(s) of the variable(s) specified in `stat.by` to include in
  the plot. Default is `NULL`.

- plot_type:

  The type of plot to create. Can be one of `"bar"`, `"rose"`, `"ring"`,
  `"pie"`, `"trend"`, `"area"`, `"dot"`, `"sankey"`, `"chord"`,
  `"venn"`, or `"upset"`.

- stat_type:

  The type of statistic to compute for the plot. Can be one of
  `"percent"` or `"count"`.

- position:

  The position adjustment for the plot. Can be one of `"stack"` or
  `"dodge"`.

- palette:

  The name of the color palette to use for the plot.

- palcolor:

  The color to use in the color palette.

- alpha:

  The transparency level for the plot.

- bg_palette:

  The name of the background color palette to use for bar plots.

- bg_palcolor:

  The color to use in the background color palette.

- bg_alpha:

  The transparency level for the background color in bar plots.

- label:

  Whether to add labels on the plot. Default is `FALSE`.

- label.size:

  The size of the labels.

- label.fg:

  The foreground color of the labels.

- label.bg:

  The background color of the labels.

- label.bg.r:

  The radius of the rounded corners of the label background.

- aspect.ratio:

  The aspect ratio of the plot.

- title:

  The main title of the plot.

- subtitle:

  The subtitle of the plot.

- xlab:

  The x-axis label of the plot.

- ylab:

  The y-axis label of the plot.

- legend.position:

  The position of the legend in the plot. Can be one of `"right"`,
  `"left"`, `"bottom"`, `"top"`, or `"none"`.

- legend.direction:

  The direction of the legend in the plot. Can be one of `"vertical"` or
  `"horizontal"`.

- theme_use:

  The name of the theme to use for the plot. Can be one of the
  predefined themes or a custom theme.

- theme_args:

  A list of arguments to be passed to the theme function.

- combine:

  Whether to combine multiple plots into a single plot. Default is TRUE.

- nrow:

  The number of rows in the combined plot. Default is NULL.

- ncol:

  The number of columns in the combined plot. Default is NULL.

- byrow:

  Whether to fill the plot by row or by column. Default is TRUE.

- force:

  Whether to force the plot even if some variables have more than 100
  levels. Default is FALSE.

- seed:

  The random seed to use for reproducible results. Default is `11`.

## See also

[StatPlot](https://mengxu98.github.io/scop/reference/StatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-19 14:11:55] Start standard scop workflow...
#> ℹ [2025-11-19 14:11:55] Checking a list of <Seurat> object...
#> ! [2025-11-19 14:11:55] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:11:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:11:57] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:11:57] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:11:58] Number of available HVF: 2000
#> ℹ [2025-11-19 14:11:58] Finished check
#> ℹ [2025-11-19 14:11:58] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 14:11:58] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 14:11:59] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 14:11:59] Reorder clusters...
#> ℹ [2025-11-19 14:11:59] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 14:11:59] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:11:59] UMAP will return its model
#> ℹ [2025-11-19 14:12:02] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:12:02] UMAP will return its model
#> ✔ [2025-11-19 14:12:05] Run scop standard workflow done
p1 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  label = TRUE
)
p1


thisplot::panel_fix(p1, height = 2, width = 3)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  stat_type = "count",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  bg.by = "CellType",
  palette = "Set1",
  stat_type = "count",
  position = "dodge"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "bar"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "rose"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "ring"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "pie"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "dot"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "rose"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "ring"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "area"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "dot"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "trend"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar",
  individual = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "bar"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "area"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "dot"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "trend"
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "bar",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "rose",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  plot_type = "ring",
  position = "dodge",
  label = TRUE
)


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "sankey"
)
#> ! [2025-11-19 14:12:12] Stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset
#> ◌ [2025-11-19 14:12:12] Installing: ggsankey...
#>  
#> → Will install 2 packages.
#> → All 2 packages (0 B) are cached.
#> + forcats    1.0.1     
#> + ggsankey   0.0.99999 [bld][cmp] (GitHub: b675d0d)
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 2 pkgs are cached
#> ✔ Got forcats 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (422.23 kB)
#> ✔ Got ggsankey 0.0.99999 (source) (277.80 kB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 31 not upgraded.
#> ✔ Installed forcats 1.0.1  (1.2s)
#> ℹ Packaging ggsankey 0.0.99999
#> ✔ Packaged ggsankey 0.0.99999 (525ms)
#> ℹ Building ggsankey 0.0.99999
#> ✔ Built ggsankey 0.0.99999 (2.2s)
#> ✔ Installed ggsankey 0.0.99999 (github::davidsjoberg/ggsankey@b675d0d) (1s)
#> ✔ 1 pkg + 29 deps: kept 28, added 2, dld 2 (NA B) [9.1s]
#> ✔ [2025-11-19 14:12:21] davidsjoberg/ggsankey installed successfully


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "chord"
)
#> ! [2025-11-19 14:12:22] Stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset


CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "venn",
  stat_level = list(
    CellType = c("Ductal", "Ngn3-low-EP"),
    Phase = "S"
  )
)
#> ! [2025-11-19 14:12:23] Stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset
#> ◌ [2025-11-19 14:12:23] Installing: ggVennDiagram...
#>  
#> → Will install 8 packages.
#> → All 8 packages (0 B) are cached.
#> + admisc          0.39  
#> + aplot           0.2.9 
#> + ggVennDiagram   1.5.4 
#> + ggfun           0.2.0 
#> + ggplotify       0.1.3 
#> + gridGraphics    0.5-1 
#> + venn            1.12  
#> + yulab.utils     0.2.1 
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 8 pkgs are cached
#> ✔ Got aplot 0.2.9 (x86_64-pc-linux-gnu-ubuntu-24.04) (105.37 kB)
#> ✔ Got ggplotify 0.1.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (142.82 kB)
#> ✔ Got yulab.utils 0.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (119.17 kB)
#> ✔ Got gridGraphics 0.5-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (249.33 kB)
#> ✔ Got admisc 0.39 (x86_64-pc-linux-gnu-ubuntu-24.04) (374.81 kB)
#> ✔ Got ggfun 0.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (250.63 kB)
#> ✔ Got venn 1.12 (x86_64-pc-linux-gnu-ubuntu-24.04) (308.02 kB)
#> ✔ Got ggVennDiagram 1.5.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (5.27 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install make`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> make is already the newest version (4.3-4.1build2).
#> 0 upgraded, 0 newly installed, 0 to remove and 31 not upgraded.
#> ✔ Installed admisc 0.39  (60ms)
#> ✔ Installed aplot 0.2.9  (72ms)
#> ✔ Installed ggfun 0.2.0  (93ms)
#> ✔ Installed ggplotify 0.1.3  (125ms)
#> ✔ Installed ggVennDiagram 1.5.4  (102ms)
#> ✔ Installed gridGraphics 0.5-1  (63ms)
#> ✔ Installed venn 1.12  (60ms)
#> ✔ Installed yulab.utils 0.2.1  (44ms)
#> ✔ 1 pkg + 36 deps: kept 29, added 8, dld 8 (6.82 MB) [5s]
#> ✔ [2025-11-19 14:12:28] ggVennDiagram installed successfully


pancreas_sub$Progenitor <- pancreas_sub$CellType %in% c("Ngn3-low-EP", "Ngn3-high-EP")
pancreas_sub$G2M <- pancreas_sub$Phase == "G2M"
pancreas_sub$Fancb_Expressed <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)["Fancb", ] > 0
pancreas_sub$Dlg3_Expressed <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)["Dlg3", ] > 0
CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
  ),
  plot_type = "venn",
  stat_level = "TRUE"
)
#> ! [2025-11-19 14:12:29] Stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset
#> ✔ [2025-11-19 14:12:29] ggVennDiagram installed successfully


CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "Progenitor", "G2M", "Fancb_Expressed", "Dlg3_Expressed"
  ),
  plot_type = "upset",
  stat_level = "TRUE"
)
#> ! [2025-11-19 14:12:29] Stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset
#> ◌ [2025-11-19 14:12:29] Installing: ggupset...
#>  
#> → Will install 1 package.
#> → The package (0 B) is cached.
#> + ggupset   0.4.1 
#>   
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Got ggupset 0.4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.59 MB)
#> ✔ Installed ggupset 0.4.1  (1s)
#> ✔ 1 pkg + 21 deps: kept 21, added 1, dld 1 (2.59 MB) [2.8s]
#> ✔ [2025-11-19 14:12:32] ggupset installed successfully


sum(
  pancreas_sub$Progenitor == "FALSE" &
    pancreas_sub$G2M == "FALSE" &
    pancreas_sub$Fancb_Expressed == "TRUE" &
    pancreas_sub$Dlg3_Expressed == "FALSE"
)
#> [1] 6
```
