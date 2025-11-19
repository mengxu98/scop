# Visualize feature values on a 2-dimensional reduction plot

Plotting cell points on a reduced 2D plane and coloring according to the
values of the features.

## Usage

``` r
FeatureDimPlot(
  srt,
  features,
  reduction = NULL,
  dims = c(1, 2),
  split.by = NULL,
  cells = NULL,
  layer = "data",
  assay = NULL,
  show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
  palette = ifelse(isTRUE(compare_features), "Set1", "Spectral"),
  palcolor = NULL,
  pt.size = NULL,
  pt.alpha = 1,
  bg_cutoff = 0,
  bg_color = "grey80",
  keep_scale = "feature",
  lower_quantile = 0,
  upper_quantile = 0.99,
  lower_cutoff = NULL,
  upper_cutoff = NULL,
  add_density = FALSE,
  density_color = "grey80",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  calculate_coexp = FALSE,
  compare_features = FALSE,
  color_blend_mode = c("blend", "average", "screen", "multiply"),
  label = FALSE,
  label.size = 4,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1,
  label_insitu = FALSE,
  label_repel = FALSE,
  label_repulsion = 20,
  label_point_size = 1,
  label_point_color = "black",
  label_segment_color = "black",
  lineages = NULL,
  lineages_trim = c(0.01, 0.99),
  lineages_span = 0.75,
  lineages_palette = "Dark2",
  lineages_palcolor = NULL,
  lineages_arrow = grid::arrow(length = grid::unit(0.1, "inches")),
  lineages_linewidth = 1,
  lineages_line_bg = "white",
  lineages_line_bg_stroke = 0.5,
  lineages_whiskers = FALSE,
  lineages_whiskers_linewidth = 0.5,
  lineages_whiskers_alpha = 0.5,
  graph = NULL,
  edge_size = c(0.05, 0.5),
  edge_alpha = 0.1,
  edge_color = "grey40",
  hex = FALSE,
  hex.linewidth = 0.5,
  hex.color = "grey90",
  hex.bins = 50,
  hex.binwidth = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  aspect.ratio = 1,
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

  A Seurat object.

- features:

  A character vector or a named list of features to plot. Features can
  be gene names in Assay or names of numeric columns in meta.data.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a two-length numeric vector specifying x-
  and y-dimensions.

- split.by:

  Name of a column in meta.data to split plot by.

- cells:

  Subset cells to plot.

- layer:

  Which layer to pull expression data from? Default is `data`.

- assay:

  Which assay to pull expression data from. If `NULL`, will use the
  assay returned by
  [SeuratObject::DefaultAssay](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html).

- show_stat:

  Whether to show statistical information on the plot.

- palette:

  Name of a color palette name collected in scop.

- palcolor:

  Custom colors used to create a color palette.

- pt.size:

  Point size for plotting.

- pt.alpha:

  Point transparency.

- bg_cutoff:

  Background cutoff. Points with feature values lower than the cutoff
  will be considered as background and will be colored with `bg_color`.

- bg_color:

  Color value for background points.

- keep_scale:

  How to handle the color scale across multiple plots. Options are:

  - `NULL (no scaling):` Each individual plot is scaled to the maximum
    expression value of the feature in the condition provided to
    'split.by'. Be aware setting NULL will result in color scales that
    are not comparable between plots.

  - `"feature" (default; by row/feature scaling):` The plots for each
    individual feature are scaled to the maximum expression of the
    feature across the conditions provided to 'split.by'.

  - `"all" (universal scaling):` The plots for all features and
    conditions are scaled to the maximum expression value for the
    feature with the highest overall expression.

- lower_quantile, upper_quantile, lower_cutoff, upper_cutoff:

  Vector of minimum and maximum cutoff values or quantile values for
  each feature.

- add_density:

  Whether to add a density layer on the plot.

- density_color:

  Color of the density contours lines.

- density_filled:

  Whether to add filled contour bands instead of contour lines.

- density_filled_palette:

  Color palette used to fill contour bands.

- density_filled_palcolor:

  Custom colors used to fill contour bands.

- cells.highlight:

  A vector of cell names to highlight.

- cols.highlight:

  Color used to highlight the cells.

- sizes.highlight:

  Size of highlighted cells.

- alpha.highlight:

  Transparency of highlighted cell points.

- stroke.highlight:

  Border width of highlighted cell points.

- calculate_coexp:

  Whether to calculate the co-expression value (geometric mean) of the
  features.

- compare_features:

  Whether to show the values of multiple features on a single plot.

- color_blend_mode:

  Blend mode to use when `compare_features = TRUE`

- label:

  Whether the feature name is labeled in the center of the location of
  cells wieh high expression.

- label.size:

  Size of labels.

- label.fg:

  Foreground color of label.

- label.bg:

  Background color of label.

- label.bg.r:

  Background ratio of label.

- label_insitu:

  Whether the labels is feature names instead of numbers. Valid only
  when `compare_features = TRUE`.

- label_repel:

  Logical value indicating whether the label is repel away from the
  center location.

- label_repulsion:

  Force of repulsion between overlapping text labels. Defaults to 20.

- label_point_size:

  Size of the center points.

- label_point_color:

  Color of the center points

- label_segment_color:

  Color of the line segment for labels.

- lineages:

  Lineages/pseudotime to add to the plot. If specified, curves will be
  fitted using [stats::loess](https://rdrr.io/r/stats/loess.html)
  method.

- lineages_trim:

  Trim the leading and the trailing data in the lineages.

- lineages_span:

  The parameter α which controls the degree of smoothing in
  [stats::loess](https://rdrr.io/r/stats/loess.html) method.

- lineages_palette:

  Color palette used for lineages.

- lineages_palcolor:

  Custom colors used for lineages.

- lineages_arrow:

  Set arrows of the lineages. See
  [grid::arrow](https://rdrr.io/r/grid/arrow.html).

- lineages_linewidth:

  Width of fitted curve lines for lineages.

- lineages_line_bg:

  Background color of curve lines for lineages.

- lineages_line_bg_stroke:

  Border width of curve lines background.

- lineages_whiskers:

  Whether to add whiskers for lineages.

- lineages_whiskers_linewidth:

  Width of whiskers for lineages.

- lineages_whiskers_alpha:

  Transparency of whiskers for lineages.

- graph:

  Specify the graph name to add edges between cell neighbors to the
  plot.

- edge_size:

  Size of edges.

- edge_alpha:

  Transparency of edges.

- edge_color:

  Color of edges.

- hex:

  Whether to chane the plot type from point to the hexagonal bin.

- hex.linewidth:

  Border width of hexagonal bins.

- hex.color:

  Border color of hexagonal bins.

- hex.bins:

  Number of hexagonal bins.

- hex.binwidth:

  Hexagonal bin width.

- raster:

  Convert points to raster format, default is NULL which automatically
  rasterizes if plotting more than 100,000 cells

- raster.dpi:

  Pixel resolution for rasterized plots, passed to geom_scattermore().
  Default is `c(512, 512)`.

- aspect.ratio:

  Aspect ratio of the panel.

- title:

  The text for the title.

- subtitle:

  The text for the subtitle for the plot which will be displayed below
  the title.

- xlab:

  x-axis label.

- ylab:

  y-axis label.

- legend.position:

  The position of legends ("none", "left", "right", "bottom", "top").

- legend.direction:

  Layout of items in legends ("horizontal" or "vertical")

- theme_use:

  Theme used. Can be a character string or a theme function. For
  example, `"theme_blank"` or ggplot2::theme_classic.

- theme_args:

  Other arguments passed to the `theme_use`.

- combine:

  Combine plots into a single `patchwork` object. If `FALSE`, return a
  list of ggplot objects.

- nrow:

  Number of rows in the combined plot.

- ncol:

  Number of columns in the combined plot.

- byrow:

  Logical value indicating if the plots should be arrange by row
  (default) or by column.

- force:

  Whether to force drawing regardless of the number of features greater
  than 100.

- seed:

  Random seed set for reproducibility

## See also

[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-19 14:24:59] Start standard scop workflow...
#> ℹ [2025-11-19 14:25:00] Checking a list of <Seurat> object...
#> ! [2025-11-19 14:25:00] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:25:00] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:25:02] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:25:03] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:25:03] Number of available HVF: 2000
#> ℹ [2025-11-19 14:25:03] Finished check
#> ℹ [2025-11-19 14:25:03] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 14:25:03] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 14:25:04] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 14:25:04] Reorder clusters...
#> ℹ [2025-11-19 14:25:05] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 14:25:05] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:25:05] UMAP will return its model
#> ℹ [2025-11-19 14:25:08] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:25:08] UMAP will return its model
#> ✔ [2025-11-19 14:25:12] Run scop standard workflow done
FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score", reduction = "UMAP"
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP",
  bg_cutoff = -Inf
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP",
  theme_use = "theme_blank"
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP",
  theme_use = ggplot2::theme_classic,
  theme_args = list(base_size = 16)
)


FeatureDimPlot(
  pancreas_sub,
  features = "G2M_score",
  reduction = "UMAP"
) |> thisplot::panel_fix(
  height = 2,
  raster = TRUE,
  dpi = 30
)


# Label and highlight cell points
FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  label = TRUE,
  cells.highlight = colnames(
    pancreas_sub
  )[pancreas_sub$SubCellType == "Delta"]
)


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  split.by = "Phase",
  reduction = "UMAP",
  cells.highlight = TRUE,
  theme_use = "theme_blank"
)


# Add a density layer
FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  label = TRUE,
  add_density = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  label = TRUE,
  add_density = TRUE,
  density_filled = TRUE
)


# Chane the plot type from point to the hexagonal bin
FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  hex = TRUE
)
#> ✔ [2025-11-19 14:25:15] hexbin installed successfully


FeatureDimPlot(
  pancreas_sub,
  features = "Rbp4",
  reduction = "UMAP",
  hex = TRUE,
  hex.bins = 20
)
#> ✔ [2025-11-19 14:25:15] hexbin installed successfully


# Show lineages on the plot based on the pseudotime
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)


FeatureDimPlot(
  pancreas_sub,
  features = "Lineage2",
  reduction = "UMAP",
  lineages = "Lineage2"
)


FeatureDimPlot(
  pancreas_sub,
  features = "Lineage2",
  reduction = "UMAP",
  lineages = "Lineage2",
  lineages_whiskers = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "Lineage2",
  reduction = "UMAP",
  lineages = "Lineage2",
  lineages_span = 0.1
)


# Input a named feature list
markers <- list(
  "Ductal" = c("Sox9", "Anxa2", "Bicc1"),
  "EPs" = c("Neurog3", "Hes6"),
  "Pre-endocrine" = c("Fev", "Neurod1"),
  "Endocrine" = c("Rbp4", "Pyy"),
  "Beta" = "Ins1",
  "Alpha" = "Gcg",
  "Delta" = "Sst",
  "Epsilon" = "Ghrl"
)
FeatureDimPlot(
  pancreas_sub,
  features = markers,
  reduction = "UMAP",
  theme_use = "theme_blank",
  theme_args = list(
    plot.subtitle = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 8)
  )
)


# Plot multiple features with different scales
endocrine_markers <- c(
  "Beta" = "Ins1",
  "Alpha" = "Gcg",
  "Delta" = "Sst",
  "Epsilon" = "Ghrl"
)
FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP"
)


FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP",
  lower_quantile = 0,
  upper_quantile = 0.8
)


FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP",
  lower_cutoff = 1,
  upper_cutoff = 4
)


FeatureDimPlot(
  pancreas_sub,
  endocrine_markers,
  reduction = "UMAP",
  keep_scale = "all"
)


FeatureDimPlot(
  pancreas_sub,
  c("Delta" = "Sst", "Epsilon" = "Ghrl"),
  split.by = "Phase",
  reduction = "UMAP",
  keep_scale = "feature"
)


# Plot multiple features on one picture
FeatureDimPlot(
  pancreas_sub,
  features = endocrine_markers,
  pt.size = 1,
  compare_features = TRUE,
  color_blend_mode = "blend",
  label = TRUE,
  label_insitu = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "blend",
  title = "blend",
  label = TRUE,
  label_insitu = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "average",
  title = "average",
  label = TRUE,
  label_insitu = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "screen",
  title = "screen",
  label = TRUE,
  label_insitu = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = c("S_score", "G2M_score"),
  pt.size = 1,
  palcolor = c("red", "green"),
  compare_features = TRUE,
  color_blend_mode = "multiply",
  title = "multiply",
  label = TRUE,
  label_insitu = TRUE
)
```
