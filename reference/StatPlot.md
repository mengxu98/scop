# Statistic Plot

Visualizes data using various plot types such as bar plots, rose plots,
ring plots, pie charts, trend plots, area plots, dot plots, sankey
plots, chord plots, venn diagrams, and upset plots.

## Usage

``` r
StatPlot(
  meta.data,
  stat.by,
  group.by = NULL,
  split.by = NULL,
  bg.by = NULL,
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

- meta.data:

  The data frame containing the data to be plotted.

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

[CellStatPlot](https://mengxu98.github.io/scop/reference/CellStatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-19 15:13:01] Start standard scop workflow...
#> ℹ [2025-11-19 15:13:02] Checking a list of <Seurat> object...
#> ! [2025-11-19 15:13:02] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 15:13:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:13:04] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:13:05] Use the separate HVF from srt_list
#> ℹ [2025-11-19 15:13:05] Number of available HVF: 2000
#> ℹ [2025-11-19 15:13:05] Finished check
#> ℹ [2025-11-19 15:13:05] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 15:13:06] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 15:13:07] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 15:13:07] Reorder clusters...
#> ℹ [2025-11-19 15:13:07] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 15:13:07] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:13:07] UMAP will return its model
#> ℹ [2025-11-19 15:13:12] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:13:12] UMAP will return its model
#> ✔ [2025-11-19 15:13:16] Run scop standard workflow done
head(pancreas_sub@meta.data)
#>                     orig.ident nCount_RNA nFeature_RNA     S_score  G2M_score
#> AAACCTGAGCCTTGAT SeuratProject       7071         2613 -0.01470664 -0.2326104
#> AAACCTGGTAAGTGGC SeuratProject       5026         2211 -0.17998111 -0.1260295
#> AAACGGGAGATATGGT SeuratProject       6194         2486 -0.16794573 -0.1666881
#> AAACGGGCAAAGAATC SeuratProject       6370         2581 -0.17434505 -0.2216024
#> AAACGGGGTACAGTTC SeuratProject       9182         2906 -0.18656224 -0.1571970
#> AAACGGGTCAGCTCTC SeuratProject       5892         2282 -0.12347911 -0.2298337
#>                  nCount_spliced nFeature_spliced nCount_unspliced
#> AAACCTGAGCCTTGAT           7071             2613              953
#> AAACCTGGTAAGTGGC           5026             2211             1000
#> AAACGGGAGATATGGT           6194             2486              721
#> AAACGGGCAAAGAATC           6370             2581             1544
#> AAACGGGGTACAGTTC           9182             2906             2262
#> AAACGGGTCAGCTCTC           5892             2282              935
#>                  nFeature_unspliced     CellType  SubCellType Phase
#> AAACCTGAGCCTTGAT                638       Ductal       Ductal    G1
#> AAACCTGGTAAGTGGC                623 Ngn3-high-EP Ngn3-high-EP    G1
#> AAACGGGAGATATGGT                550       Ductal       Ductal    G1
#> AAACGGGCAAAGAATC               1015    Endocrine         Beta    G1
#> AAACGGGGTACAGTTC               1096    Endocrine         Beta    G1
#> AAACGGGTCAGCTCTC                712  Ngn3-low-EP  Ngn3-low-EP    G1
#>                  Standardpca_SNN_res.0.6 ident Standardpcaclusters
#> AAACCTGAGCCTTGAT                       0  <NA>                <NA>
#> AAACCTGGTAAGTGGC                       3  <NA>                <NA>
#> AAACGGGAGATATGGT                       0  <NA>                <NA>
#> AAACGGGCAAAGAATC                       1  <NA>                <NA>
#> AAACGGGGTACAGTTC                       1  <NA>                <NA>
#> AAACGGGTCAGCTCTC                       0  <NA>                <NA>
#>                  Standardclusters
#> AAACCTGAGCCTTGAT             <NA>
#> AAACCTGGTAAGTGGC             <NA>
#> AAACGGGAGATATGGT             <NA>
#> AAACGGGCAAAGAATC             <NA>
#> AAACGGGGTACAGTTC             <NA>
#> AAACGGGTCAGCTCTC             <NA>
StatPlot(
  pancreas_sub@meta.data,
  stat.by = "Phase",
  group.by = "CellType",
  plot_type = "bar",
  label = TRUE
)


StatPlot(
  pancreas_sub[["RNA"]]@meta.data,
  stat.by = "highly_variable_genes",
  plot_type = "ring",
  label = TRUE,
  NA_stat = FALSE
)


if (FALSE) { # \dontrun{
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  IDtype = "symbol",
  db = c("CSPA", "TF")
)
StatPlot(
  GetFeaturesData(pancreas_sub, "RNA"),
  stat.by = "TF",
  group.by = "CSPA",
  stat_type = "count",
  plot_type = "bar",
  position = "dodge",
  label = TRUE,
  NA_stat = FALSE
)
} # }
```
