# Volcano Plot

Generate a volcano plot based on differential expression analysis
results.

## Usage

``` r
VolcanoPlot(
  srt,
  group_by = NULL,
  test.use = "wilcox",
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  x_metric = "diff_pct",
  palette = "RdBu",
  palcolor = NULL,
  pt.size = 1,
  pt.alpha = 1,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  nlabel = 5,
  features_label = NULL,
  label.fg = "black",
  label.bg = "white",
  label.bg.r = 0.1,
  label.size = 4,
  aspect.ratio = NULL,
  xlab = x_metric,
  ylab = "-log10(p-adjust)",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE
)
```

## Arguments

- srt:

  An object of class `Seurat` containing the results of differential
  expression analysis.

- group_by:

  A character vector specifying the column in `srt` to group the samples
  by. Default is `NULL`.

- test.use:

  A character string specifying the type of statistical test to use.
  Default is `"wilcox"`.

- DE_threshold:

  A character string specifying the threshold for differential
  expression. Default is `"avg_log2FC > 0 & p_val_adj < 0.05"`.

- x_metric:

  A character string specifying the metric to use for the x-axis.
  Default is `"diff_pct"`.

- palette:

  A character string specifying the color palette to use for the plot.
  Default is `"RdBu"`.

- palcolor:

  A character string specifying the color for the palette. Default is
  `NULL`.

- pt.size:

  The size of the points. Default is `1`.

- pt.alpha:

  The transparency of the points. Default is `1`.

- cols.highlight:

  A character string specifying the color for highlighted points.
  Default is `"black"`.

- sizes.highlight:

  The size of the highlighted points. Default is `1`.

- alpha.highlight:

  The transparency of the highlighted points. Default is `1`.

- stroke.highlight:

  The stroke width for the highlighted points. Default is `0.5`.

- nlabel:

  An integer value specifying the number of labeled points per group.
  Default is `5`.

- features_label:

  A character vector specifying the feature labels to plot. Default is
  `NULL`.

- label.fg:

  A character string specifying the color for the labels' foreground.
  Default is `"black"`.

- label.bg:

  A character string specifying the color for the labels' background.
  Default is `"white"`.

- label.bg.r:

  The radius of the rounding of the labels' background. Default is 0.1.

- label.size:

  The size of the labels. Default is `4`.

- aspect.ratio:

  The aspect ratio of the plot. Default is `NULL`.

- xlab:

  A character string specifying the x-axis label. Default is the value
  of `x_metric`.

- ylab:

  A character string specifying the y-axis label. Default is
  `"-log10(p-adjust)"`.

- theme_use:

  A character string specifying the theme to use for the plot. Default
  is `"theme_scop"`.

- theme_args:

  A list of theme arguments to pass to the `theme_use` function. Default
  is an empty list.

- combine:

  Whether to combine the plots for each group into a single plot.
  Default is `TRUE`.

- nrow:

  An integer value specifying the number of rows in the combined plot.
  Default is `NULL`.

- ncol:

  An integer value specifying the number of columns in the combined
  plot. Default is `NULL`.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

## See also

[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md)

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
pancreas_sub <- RunDEtest(
 pancreas_sub,
 group_by = "CellType"
)
#> ⠙ [2025-12-02 05:04:02] Running [1/5] Processing: Ductal  ETA:  1s
#> ✔ [2025-12-02 05:04:02] Completed 5 tasks in 723ms
#> 
VolcanoPlot(
  pancreas_sub,
  group_by = "CellType",
  ncol = 2
)


VolcanoPlot(
  pancreas_sub,
  group_by = "CellType",
  DE_threshold = "abs(diff_pct) > 0.3 & p_val_adj < 0.05",
  ncol = 2
)


VolcanoPlot(
  pancreas_sub,
  group_by = "CellType",
  x_metric = "avg_log2FC",
  ncol = 2
)
```
