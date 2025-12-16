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
  `NULL`.

- split.by:

  A character vector specifying the variable to split the plot by.
  Default is `NULL`.

- bg.by:

  A character vector specifying the variable to use as the background
  color. Default is `NULL`.

- plot.by:

  A character vector specifying how to plot the data, by group or
  feature. Possible values are `"group"` or `"feature"`. Default is
  `"group"`.

- fill.by:

  A string specifying what to fill the plot by. Possible values are
  `"group"`, `"feature"`, or `"expression"`. Default is `"group"`.

- cells:

  A character vector specifying the cells to include in the plot.
  Default is `NULL`.

- layer:

  A string specifying which layer of the Seurat object to use. Default
  is `"data"`.

- assay:

  A string specifying which assay to use. Default is `NULL`.

- keep_empty:

  Whether to keep empty levels in the plot. Default is `FALSE`.

- individual:

  Whether to create individual plots for each group. Default is `FALSE`.

- plot_type:

  A string specifying the type of plot to create. Possible values are
  `"violin"`, `"box"`, `"bar"`, `"dot"`, or `"col"`. Default is
  `"violin"`.

- palette:

  A string specifying the color palette to use for filling. Default is
  `"Paired"`.

- palcolor:

  A character vector specifying specific colors to use for filling.
  Default is `NULL`.

- alpha:

  The transparency of the plot. Default is `1`.

- bg_palette:

  A string specifying the color palette to use for the background.
  Default is `"Paired"`.

- bg_palcolor:

  A character vector specifying specific colors to use for the
  background. Default is `NULL`.

- bg_alpha:

  The transparency of the background. Default is `0.2`.

- add_box:

  Whether to add a box plot to the plot. Default is `FALSE`.

- box_color:

  A string specifying the color of the box plot. Default is `"black"`.

- box_width:

  The width of the box plot. Default is `0.1`.

- box_ptsize:

  The size of the points of the box plot. Default is `2`.

- add_point:

  Whether to add individual data points to the plot. Default is `FALSE`.

- pt.color:

  A string specifying the color of the data points. Default is
  `"grey30"`.

- pt.size:

  The size of the data points. If NULL, the size is automatically
  determined. Default is `NULL`.

- pt.alpha:

  The transparency of the data points. Default is `1`.

- jitter.width:

  The width of the jitter. Default is `0.5`.

- jitter.height:

  The height of the jitter. Default is `0.1`.

- add_trend:

  Whether to add a trend line to the plot. Default is `FALSE`.

- trend_color:

  A string specifying the color of the trend line. Default is `"black"`.

- trend_linewidth:

  The width of the trend line. Default is `1`.

- trend_ptsize:

  The size of the points of the trend line. Default is `2`.

- add_stat:

  A string specifying which statistical summary to add to the plot.
  Possible values are `"none"`, `"mean"`, or `"median"`. Default is
  `"none"`.

- stat_color:

  A string specifying the color of the statistical summary. Default is
  `"black"`.

- stat_size:

  The size of the statistical summary. Default is `1`.

- stat_stroke:

  The stroke width of the statistical summary. Default is `1`.

- stat_shape:

  The shape of the statistical summary. Default is `25`.

- add_line:

  The y-intercept for adding a horizontal line. Default is `NULL`.

- line_color:

  A string specifying the color of the horizontal line. Default is
  `"red"`.

- line_size:

  The width of the horizontal line. Default is `1`.

- line_type:

  The type of the horizontal line. Default is `1`.

- cells.highlight:

  A logical or character vector specifying the cells to highlight in the
  plot. If TRUE, all cells are highlighted. If FALSE, no cells are
  highlighted. Default is `NULL`.

- cols.highlight:

  A string specifying the color of the highlighted cells. Default is
  `"red"`.

- sizes.highlight:

  The size of the highlighted cells. Default is `1`.

- alpha.highlight:

  The transparency of the highlighted cells. Default is `1`.

- calculate_coexp:

  Whether to calculate co-expression values. Default is `FALSE`.

- same.y.lims:

  Whether to use the same y-axis limits for all plots. Default is
  `FALSE`.

- y.min:

  A numeric or character value specifying the minimum y-axis limit. If a
  character value is provided, it must be of the form "qN" where N is a
  number between 0 and 100 (inclusive) representing the quantile to use
  for the limit. Default is `NULL`.

- y.max:

  A numeric or character value specifying the maximum y-axis limit. If a
  character value is provided, it must be of the form "qN" where N is a
  number between 0 and 100 (inclusive) representing the quantile to use
  for the limit. Default is `NULL`.

- y.trans:

  A string specifying the transformation to apply to the y-axis.
  Possible values are `"identity"` or `"log2"`. Default is `"identity"`.

- y.nbreaks:

  A number of breaks to use for the y-axis. Default is `5`.

- sort:

  A logical or character value specifying whether to sort the groups on
  the x-axis. If TRUE, groups are sorted in increasing order. If FALSE,
  groups are not sorted. If "increasing", groups are sorted in
  increasing order. If "decreasing", groups are sorted in decreasing
  order. Default is `FALSE`.

- stack:

  A logical specifying whether to stack the plots on top of each other.
  Default is `FALSE`.

- flip:

  A logical specifying whether to flip the plot vertically. Default is
  `FALSE`.

- comparisons:

  A list of length-2 vectors. The entries in the vector are either the
  names of 2 values on the x-axis or the 2 integers that correspond to
  the index of the groups of interest, to be compared.

- ref_group:

  A string specifying the reference group for pairwise comparisons.
  Default is `NULL`.

- pairwise_method:

  Method to use for pairwise comparisons. Default is `"wilcox.test"`.

- multiplegroup_comparisons:

  Whether to add multiple group comparisons to the plot. Default is
  `FALSE`.

- multiple_method:

  Method to use for multiple group comparisons. Default is
  `"kruskal.test"`.

- sig_label:

  A string specifying the label to use for significant comparisons.
  Possible values are `"p.signif"` or `"p.format"`. Default is
  `"p.format"`.

- sig_labelsize:

  The size of the significant comparison labels. Default is `3.5`.

- aspect.ratio:

  The aspect ratio of the plot. Default is `NULL`.

- title:

  A string specifying the title of the plot. Default is `NULL`.

- subtitle:

  A string specifying the subtitle of the plot. Default is `NULL`.

- xlab:

  A string specifying the label of the x-axis. Default is `NULL`.

- ylab:

  A string specifying the label of the y-axis. Default is
  `"Expression level"`.

- legend.position:

  A string specifying the position of the legend. Possible values are
  `"right"`, `"left"`, `"top"`, `"bottom"`, or `"none"`. Default is
  `"right"`.

- legend.direction:

  A string specifying the direction of the legend. Possible values are
  `"vertical"` or `"horizontal"`. Default is `"vertical"`.

- theme_use:

  A string specifying the theme to use for the plot. Default is
  `"theme_scop"`.

- theme_args:

  A list of arguments to pass to the theme function. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Whether to combine the individual plots into a single plot. Default is
  `TRUE`.

- nrow:

  A number of rows for the combined plot. Default is `NULL`.

- ncol:

  A number of columns for the combined plot. Default is `NULL`.

- byrow:

  Whether to fill the combined plot by row or by column. Default is
  `TRUE`.

- force:

  Whether to force the plot creation even if there are more than 100
  levels in a variable. Default is `FALSE`.

- seed:

  An integer specifying the random seed to use for generating jitter.
  Default is `11`.

## See also

[CellStatPlot](https://mengxu98.github.io/scop/reference/CellStatPlot.md),
[StatPlot](https://mengxu98.github.io/scop/reference/StatPlot.md)

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
FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType"
) |> thisplot::panel_fix(height = 1, width = 2)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "box"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "bar"
)
#> Warning: Computation failed in `stat_summary()`.
#> Caused by error in `fun.data()`:
#> ! The package "Hmisc" is required.
#> Warning: Computation failed in `stat_summary()`.
#> Caused by error in `fun.data()`:
#> ! The package "Hmisc" is required.
#> Warning: Computation failed in `stat_summary()`.
#> Caused by error in `fun.data()`:
#> ! The package "Hmisc" is required.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Computation failed in `stat_summary()`.
#> Caused by error in `fun.data()`:
#> ! The package "Hmisc" is required.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "dot"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.

FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  plot_type = "col"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_box = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_point = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.

FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_trend = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_stat = "mean"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  add_line = 0.2,
  line_type = 2
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  split.by = "Phase"
)
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  split.by = "Phase",
  add_box = TRUE,
  add_trend = TRUE
)
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score", "Fev"),
  group.by = "SubCellType",
  split.by = "Phase",
  comparisons = TRUE
)
#>  
#> → Will install 34 packages.
#> → Will download 3 CRAN packages (4.41 MB), cached: 31 (0 B).
#> + Deriv            4.2.0      
#> + Formula          1.2-5      
#> + MatrixModels     0.5-4      
#> + Rdpack           2.6.4      
#> + SparseM          1.84-2     
#> + TTR              0.24.4     
#> + car              3.1-3      
#> + carData          3.0-5      
#> + colorspace       2.1-2      
#> + corrplot         0.95       
#> + doBy             4.7.1      
#> + forecast         8.24.0     [bld][cmp][dl] (581.95 kB)
#> + fracdiff         1.5-3      
#> + ggpubr           0.6.2      
#> + ggsci            4.1.0      
#> + ggsignif         0.6.4      
#> + lme4             1.1-38     [bld][cmp][dl] (3.77 MB)
#> + microbenchmark   1.5.0      
#> + minqa            1.2.8      [bld][cmp][dl] (54.64 kB) + ✔ make
#> + modelr           0.1.11     
#> + nloptr           2.2.1       + ✔ cmake
#> + numDeriv         2016.8-1.1 
#> + pbkrtest         0.5.5      
#> + polynom          1.4-1      
#> + quadprog         1.5-8      
#> + quantmod         0.4.28     
#> + quantreg         6.1        
#> + rbibutils        2.4        
#> + reformulas       0.4.2      
#> + rstatix          0.7.3      
#> + timeDate         4051.111   
#> + tseries          0.10-58    
#> + urca             1.3-4      
#> + xts              0.14.1     
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 3 pkgs (4.41 MB), 31 cached
#> ✔ Cached copy of forecast 8.24.0 (source) is the latest build
#> ✔ Cached copy of lme4 1.1-38 (source) is the latest build
#> ✔ Got modelr 0.1.11 (x86_64-pc-linux-gnu-ubuntu-24.04) (200.70 kB)
#> ✔ Got carData 3.0-5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.70 MB)
#> ✔ Got car 3.1-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.54 MB)
#> ✔ Got Deriv 4.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (152.58 kB)
#> ✔ Got TTR 0.24.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (524.49 kB)
#> ✔ Got rbibutils 2.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.15 MB)
#> ✔ Got Formula 1.2-5 (x86_64-pc-linux-gnu-ubuntu-24.04) (159.13 kB)
#> ✔ Got ggpubr 0.6.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.13 MB)
#> ✔ Got numDeriv 2016.8-1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (114.36 kB)
#> ✔ Got colorspace 2.1-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.64 MB)
#> ✔ Got minqa 1.2.8 (source) (55.10 kB)
#> ✔ Got polynom 1.4-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (406.99 kB)
#> ✔ Got ggsignif 0.6.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (602.07 kB)
#> ✔ Got Rdpack 2.6.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (637.80 kB)
#> ✔ Got microbenchmark 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (65.96 kB)
#> ✔ Got corrplot 0.95 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.82 MB)
#> ✔ Got rstatix 0.7.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (615.19 kB)
#> ✔ Got pbkrtest 0.5.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (222.06 kB)
#> ✔ Got SparseM 1.84-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (887.98 kB)
#> ✔ Got ggsci 4.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.30 MB)
#> ✔ Got reformulas 0.4.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (139.23 kB)
#> ✔ Got nloptr 2.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (567.65 kB)
#> ✔ Got quantreg 6.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.46 MB)
#> ✔ Got MatrixModels 0.5-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (408.50 kB)
#> ✔ Got xts 0.14.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.22 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install make cmake libcurl4-openssl-dev libssl-dev pandoc libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> make is already the newest version (4.3-4.1build2).
#> cmake is already the newest version (3.28.3-1build7).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 49 not upgraded.
#> ℹ Building minqa 1.2.8
#> ✔ Installed car 3.1-3  (117ms)
#> ✔ Installed carData 3.0-5  (152ms)
#> ✔ Installed colorspace 2.1-2  (166ms)
#> ✔ Installed corrplot 0.95  (136ms)
#> ✔ Installed Deriv 4.2.0  (114ms)
#> ✔ Installed doBy 4.7.1  (126ms)
#> ✔ Installed forecast 8.24.0  (95ms)
#> ✔ Installed Formula 1.2-5  (71ms)
#> ✔ Installed fracdiff 1.5-3  (68ms)
#> ✔ Installed ggpubr 0.6.2  (79ms)
#> ✔ Installed ggsci 4.1.0  (83ms)
#> ✔ Installed ggsignif 0.6.4  (75ms)
#> ✔ Installed MatrixModels 0.5-4  (78ms)
#> ✔ Installed microbenchmark 1.5.0  (96ms)
#> ✔ Installed modelr 0.1.11  (1s)
#> ✔ Installed lme4 1.1-38  (1.3s)
#> ✔ Installed nloptr 2.2.1  (79ms)
#> ✔ Installed numDeriv 2016.8-1.1  (72ms)
#> ✔ Installed pbkrtest 0.5.5  (73ms)
#> ✔ Installed polynom 1.4-1  (72ms)
#> ✔ Installed quadprog 1.5-8  (70ms)
#> ✔ Installed quantmod 0.4.28  (107ms)
#> ✔ Installed quantreg 6.1  (108ms)
#> ✔ Installed rbibutils 2.4  (77ms)
#> ✔ Installed Rdpack 2.6.4  (73ms)
#> ✔ Installed reformulas 0.4.2  (71ms)
#> ✔ Installed rstatix 0.7.3  (75ms)
#> ✔ Installed SparseM 1.84-2  (74ms)
#> ✔ Installed timeDate 4051.111  (72ms)
#> ✔ Installed tseries 0.10-58  (103ms)
#> ✔ Installed TTR 0.24.4  (71ms)
#> ✔ Installed urca 1.3-4  (72ms)
#> ✔ Installed xts 0.14.1  (51ms)
#> ✔ Built minqa 1.2.8 (7.2s)
#> ✔ Installed minqa 1.2.8  (1s)
#> ✔ 1 pkg + 101 deps: kept 65, added 34, dld 25 (23.73 MB) [15.6s]
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta"))
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Rbp4", "Pyy"),
  group.by = "SubCellType",
  comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")),
  sig_label = "p.format"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "feature"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "feature",
  multiplegroup_comparisons = TRUE,
  sig_label = "p.format",
  sig_labelsize = 4
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4", "Ins1"),
  group.by = "CellType",
  plot.by = "feature",
  comparisons = list(c("Neurog3", "Rbp4"), c("Rbp4", "Ins1")),
  stack = TRUE
)


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
#> Warning: Different features in new layer data than already exists for scale.data
FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Neurog3", "Rbp4"),
  group.by = "CellType",
  layer = "scale.data",
  ylab = "FoldChange",
  same.y.lims = TRUE,
  y.max = 4
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
```
