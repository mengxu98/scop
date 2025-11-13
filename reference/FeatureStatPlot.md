# Statistical plot of features

This function generates a statistical plot for features.

## Usage

``` r
FeatureStatPlot(
  srt,
  stat.by,
  group.by = NULL,
  split.by = NULL,
  bg.by = NULL,
  plot.by = c("group", "feature"),
  fill.by = c("group", "feature", "expression"),
  cells = NULL,
  layer = "data",
  assay = NULL,
  keep_empty = FALSE,
  individual = FALSE,
  plot_type = c("violin", "box", "bar", "dot", "col"),
  palette = "Paired",
  palcolor = NULL,
  alpha = 1,
  bg_palette = "Paired",
  bg_palcolor = NULL,
  bg_alpha = 0.2,
  add_box = FALSE,
  box_color = "black",
  box_width = 0.1,
  box_ptsize = 2,
  add_point = FALSE,
  pt.color = "grey30",
  pt.size = NULL,
  pt.alpha = 1,
  jitter.width = 0.4,
  jitter.height = 0.1,
  add_trend = FALSE,
  trend_color = "black",
  trend_linewidth = 1,
  trend_ptsize = 2,
  add_stat = c("none", "mean", "median"),
  stat_color = "black",
  stat_size = 1,
  stat_stroke = 1,
  stat_shape = 25,
  add_line = NULL,
  line_color = "red",
  line_size = 1,
  line_type = 1,
  cells.highlight = NULL,
  cols.highlight = "red",
  sizes.highlight = 1,
  alpha.highlight = 1,
  calculate_coexp = FALSE,
  same.y.lims = FALSE,
  y.min = NULL,
  y.max = NULL,
  y.trans = "identity",
  y.nbreaks = 5,
  sort = FALSE,
  stack = FALSE,
  flip = FALSE,
  comparisons = NULL,
  ref_group = NULL,
  pairwise_method = "wilcox.test",
  multiplegroup_comparisons = FALSE,
  multiple_method = "kruskal.test",
  sig_label = c("p.signif", "p.format"),
  sig_labelsize = 3.5,
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = "Expression level",
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

  A Seurat object.

- stat.by:

  A character vector specifying the features to plot.

- group.by:

  A character vector specifying the groups to group by. Default is
  \`NULL\`.

- split.by:

  A character vector specifying the variable to split the plot by.
  Default is \`NULL\`.

- bg.by:

  A character vector specifying the variable to use as the background
  color. Default is \`NULL\`.

- plot.by:

  A character vector specifying how to plot the data, by group or
  feature. Possible values are \`"group"\` or \`"feature"\`. Default is
  \`"group"\`.

- fill.by:

  A string specifying what to fill the plot by. Possible values are
  \`"group"\`, \`"feature"\`, or \`"expression"\`. Default is
  \`"group"\`.

- cells:

  A character vector specifying the cells to include in the plot.
  Default is \`NULL\`.

- layer:

  A string specifying which layer of the Seurat object to use. Default
  is \`"data"\`.

- assay:

  A string specifying which assay to use. Default is \`NULL\`.

- keep_empty:

  Whether to keep empty levels in the plot. Default is \`FALSE\`.

- individual:

  Whether to create individual plots for each group. Default is
  \`FALSE\`.

- plot_type:

  A string specifying the type of plot to create. Possible values are
  \`"violin"\`, \`"box"\`, \`"bar"\`, \`"dot"\`, or \`"col"\`. Default
  is \`"violin"\`.

- palette:

  A string specifying the color palette to use for filling. Default is
  \`"Paired"\`.

- palcolor:

  A character vector specifying specific colors to use for filling.
  Default is \`NULL\`.

- alpha:

  The transparency of the plot. Default is \`1\`.

- bg_palette:

  A string specifying the color palette to use for the background.
  Default is \`"Paired"\`.

- bg_palcolor:

  A character vector specifying specific colors to use for the
  background. Default is \`NULL\`.

- bg_alpha:

  The transparency of the background. Default is \`0.2\`.

- add_box:

  Whether to add a box plot to the plot. Default is \`FALSE\`.

- box_color:

  A string specifying the color of the box plot. Default is \`"black"\`.

- box_width:

  The width of the box plot. Default is \`0.1\`.

- box_ptsize:

  The size of the points of the box plot. Default is \`2\`.

- add_point:

  Whether to add individual data points to the plot. Default is
  \`FALSE\`.

- pt.color:

  A string specifying the color of the data points. Default is
  \`"grey30"\`.

- pt.size:

  The size of the data points. If NULL, the size is automatically
  determined. Default is \`NULL\`.

- pt.alpha:

  The transparency of the data points. Default is \`1\`.

- jitter.width:

  The width of the jitter. Default is \`0.5\`.

- jitter.height:

  The height of the jitter. Default is \`0.1\`.

- add_trend:

  Whether to add a trend line to the plot. Default is \`FALSE\`.

- trend_color:

  A string specifying the color of the trend line. Default is
  \`"black"\`.

- trend_linewidth:

  The width of the trend line. Default is \`1\`.

- trend_ptsize:

  The size of the points of the trend line. Default is \`2\`.

- add_stat:

  A string specifying which statistical summary to add to the plot.
  Possible values are \`"none"\`, \`"mean"\`, or \`"median"\`. Default
  is \`"none"\`.

- stat_color:

  A string specifying the color of the statistical summary. Default is
  \`"black"\`.

- stat_size:

  The size of the statistical summary. Default is \`1\`.

- stat_stroke:

  The stroke width of the statistical summary. Default is \`1\`.

- stat_shape:

  The shape of the statistical summary. Default is \`25\`.

- add_line:

  The y-intercept for adding a horizontal line. Default is \`NULL\`.

- line_color:

  A string specifying the color of the horizontal line. Default is
  \`"red"\`.

- line_size:

  The width of the horizontal line. Default is \`1\`.

- line_type:

  The type of the horizontal line. Default is \`1\`.

- cells.highlight:

  A logical or character vector specifying the cells to highlight in the
  plot. If TRUE, all cells are highlighted. If FALSE, no cells are
  highlighted. Default is \`NULL\`.

- cols.highlight:

  A string specifying the color of the highlighted cells. Default is
  \`"red"\`.

- sizes.highlight:

  The size of the highlighted cells. Default is \`1\`.

- alpha.highlight:

  The transparency of the highlighted cells. Default is \`1\`.

- calculate_coexp:

  Whether to calculate co-expression values. Default is \`FALSE\`.

- same.y.lims:

  Whether to use the same y-axis limits for all plots. Default is
  \`FALSE\`.

- y.min:

  A numeric or character value specifying the minimum y-axis limit. If a
  character value is provided, it must be of the form "qN" where N is a
  number between 0 and 100 (inclusive) representing the quantile to use
  for the limit. Default is \`NULL\`.

- y.max:

  A numeric or character value specifying the maximum y-axis limit. If a
  character value is provided, it must be of the form "qN" where N is a
  number between 0 and 100 (inclusive) representing the quantile to use
  for the limit. Default is \`NULL\`.

- y.trans:

  A string specifying the transformation to apply to the y-axis.
  Possible values are \`"identity"\` or \`"log2"\`. Default is
  \`"identity"\`.

- y.nbreaks:

  A number of breaks to use for the y-axis. Default is \`5\`.

- sort:

  A logical or character value specifying whether to sort the groups on
  the x-axis. If TRUE, groups are sorted in increasing order. If FALSE,
  groups are not sorted. If "increasing", groups are sorted in
  increasing order. If "decreasing", groups are sorted in decreasing
  order. Default is \`FALSE\`.

- stack:

  A logical specifying whether to stack the plots on top of each other.
  Default is \`FALSE\`.

- flip:

  A logical specifying whether to flip the plot vertically. Default is
  \`FALSE\`.

- comparisons:

  A list of length-2 vectors. The entries in the vector are either the
  names of 2 values on the x-axis or the 2 integers that correspond to
  the index of the groups of interest, to be compared.

- ref_group:

  A string specifying the reference group for pairwise comparisons.
  Default is \`NULL\`.

- pairwise_method:

  Method to use for pairwise comparisons. Default is \`"wilcox.test"\`.

- multiplegroup_comparisons:

  Whether to add multiple group comparisons to the plot. Default is
  \`FALSE\`.

- multiple_method:

  Method to use for multiple group comparisons. Default is
  \`"kruskal.test"\`.

- sig_label:

  A string specifying the label to use for significant comparisons.
  Possible values are \`"p.signif"\` or \`"p.format"\`. Default is
  \`"p.format"\`.

- sig_labelsize:

  The size of the significant comparison labels. Default is \`3.5\`.

- aspect.ratio:

  The aspect ratio of the plot. Default is \`NULL\`.

- title:

  A string specifying the title of the plot. Default is \`NULL\`.

- subtitle:

  A string specifying the subtitle of the plot. Default is \`NULL\`.

- xlab:

  A string specifying the label of the x-axis. Default is \`NULL\`.

- ylab:

  A string specifying the label of the y-axis. Default is \`"Expression
  level"\`.

- legend.position:

  A string specifying the position of the legend. Possible values are
  \`"right"\`, \`"left"\`, \`"top"\`, \`"bottom"\`, or \`"none"\`.
  Default is \`"right"\`.

- legend.direction:

  A string specifying the direction of the legend. Possible values are
  \`"vertical"\` or \`"horizontal"\`. Default is \`"vertical"\`.

- theme_use:

  A string specifying the theme to use for the plot. Default is
  \`"theme_scop"\`.

- theme_args:

  A list of arguments to pass to the theme function. Default is an empty
  list.

- combine:

  Whether to combine the individual plots into a single plot. Default is
  \`TRUE\`.

- nrow:

  A number of rows for the combined plot. Default is \`NULL\`.

- ncol:

  A number of columns for the combined plot. Default is \`NULL\`.

- byrow:

  Whether to fill the combined plot by row or by column. Default is
  \`TRUE\`.

- force:

  Whether to force the plot creation even if there are more than 100
  levels in a variable. Default is \`FALSE\`.

- seed:

  An integer specifying the random seed to use for generating jitter.
  Default is \`11\`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 11:58:56] Start standard scop workflow...
#> ℹ [2025-11-13 11:58:56] Checking a list of <Seurat> object...
#> ! [2025-11-13 11:58:56] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 11:58:56] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 11:58:58] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 11:58:59] Use the separate HVF from srt_list
#> ℹ [2025-11-13 11:58:59] Number of available HVF: 2000
#> ℹ [2025-11-13 11:58:59] Finished check
#> ℹ [2025-11-13 11:58:59] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 11:59:00] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 11:59:01] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 11:59:01] Reorder clusters...
#> ℹ [2025-11-13 11:59:01] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 11:59:01] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 11:59:01] UMAP will return its model
#> ℹ [2025-11-13 11:59:05] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 11:59:05] UMAP will return its model
#> ✔ [2025-11-13 11:59:08] Run scop standard workflow done
FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType"
) |> thisplot::panel_fix(height = 1, width = 2)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "box"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "bar"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "dot"
)

FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "col"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_box = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_point = TRUE
)

FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_trend = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_stat = "mean"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_line = 0.2,
  line_type = 2
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  split.by = "Phase"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  split.by = "Phase",
  add_box = TRUE,
  add_trend = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  split.by = "Phase",
  comparisons = TRUE
)
#> ✔ [2025-11-13 11:59:21] ggpubr installed successfully
#> ℹ [2025-11-13 11:59:21] Detected more than 2 groups. Use "kruskal.test" for comparison
#> ℹ [2025-11-13 11:59:21] Detected more than 2 groups. Use "kruskal.test" for comparison


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  fill.by = "expression",
  palette = "Blues",
  same.y.lims = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  multiplegroup_comparisons = TRUE
)
#> ✔ [2025-11-13 11:59:23] ggpubr installed successfully


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta"))
)
#> ✔ [2025-11-13 11:59:25] ggpubr installed successfully


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")),
  sig_label = "p.format"
)
#> ✔ [2025-11-13 11:59:26] ggpubr installed successfully


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  bg.by = "CellType",
  add_box = TRUE, stack = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c(
    "Sox9", "Anxa2", "Bicc1", # Ductal
    "Neurog3", "Hes6", # EPs
    "Fev", "Neurod1", # Pre-endocrine
    "Rbp4", "Pyy", # Endocrine
    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
  ),
  legend.position = "top",
  legend.direction = "horizontal",
  group.by = "SubCellType",
  bg.by = "CellType",
  stack = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c(
    "Sox9", "Anxa2", "Bicc1", # Ductal
    "Neurog3", "Hes6", # EPs
    "Fev", "Neurod1", # Pre-endocrine
    "Rbp4", "Pyy", # Endocrine
    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
  ),
  fill.by = "feature",
  plot_type = "box",
  group.by = "SubCellType",
  bg.by = "CellType", stack = TRUE, flip = TRUE
) |> thisplot::panel_fix_overall(
  width = 8, height = 5
)

# As the plot is created by combining,
# we can adjust the overall height and width directly.

FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "group"
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "feature"
)
#> ℹ [2025-11-13 11:59:39] Setting `group.by` to "Features" as `plot.by` is set to "feature"


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "feature",
  multiplegroup_comparisons = TRUE,
  sig_label = "p.format",
  sig_labelsize = 4
)
#> ℹ [2025-11-13 11:59:40] Setting `group.by` to "Features" as `plot.by` is set to "feature"
#> ✔ [2025-11-13 11:59:40] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:40] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:41] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:41] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:41] ggpubr installed successfully


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "feature",
  comparisons = list(c("Neurog3", "Rbp4"), c("Rbp4", "Ins1")),
  stack = TRUE
)
#> ℹ [2025-11-13 11:59:43] Setting `group.by` to "Features" as `plot.by` is set to "feature"
#> ✔ [2025-11-13 11:59:43] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:43] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:44] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:44] ggpubr installed successfully
#> ✔ [2025-11-13 11:59:44] ggpubr installed successfully


FeatureStatPlot(pancreas_sub,
  stat.by = c(
    "Sox9", "Anxa2", "Bicc1", # Ductal
    "Neurog3", "Hes6", # EPs
    "Fev", "Neurod1", # Pre-endocrine
    "Rbp4", "Pyy", # Endocrine
    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
  ), group.by = "SubCellType",
  plot.by = "feature",
  stack = TRUE
)
#> ℹ [2025-11-13 11:59:46] Setting `group.by` to "Features" as `plot.by` is set to "feature"


data <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "data"
)
pancreas_sub <- SeuratObject::SetAssayData(
  object = pancreas_sub,
  layer = "scale.data",
  assay = "RNA",
  new.data = data / Matrix::rowMeans(data)
)
FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4"),
  group.by = "CellType",
  layer = "scale.data",
  ylab = "FoldChange",
  same.y.lims = TRUE,
  y.max = 4
)
```
