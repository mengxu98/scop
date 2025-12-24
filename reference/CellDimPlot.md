# Cell Dimensional Plot

Visualize cell groups on a 2-dimensional reduction plot. Plotting cell
points on a reduced 2D plane and coloring according to the groups.

## Usage

``` r
CellDimPlot(
  srt,
  group.by,
  reduction = NULL,
  dims = c(1, 2),
  split.by = NULL,
  cells = NULL,
  show_na = FALSE,
  show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Paired",
  palcolor = NULL,
  bg_color = "grey80",
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
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  add_density = FALSE,
  density_color = "grey80",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  add_mark = FALSE,
  mark_type = c("hull", "ellipse", "rect", "circle"),
  mark_expand = grid::unit(3, "mm"),
  mark_alpha = 0.1,
  mark_linetype = 1,
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
  stat.by = NULL,
  stat_type = "percent",
  stat_plot_type = "pie",
  stat_plot_position = c("stack", "dodge"),
  stat_plot_size = 0.15,
  stat_plot_palette = "Set1",
  stat_palcolor = NULL,
  stat_plot_alpha = 1,
  stat_plot_label = FALSE,
  stat_plot_label_size = 3,
  graph = NULL,
  edge_size = c(0.05, 0.5),
  edge_alpha = 0.1,
  edge_color = "grey40",
  paga = NULL,
  paga_type = "connectivities",
  paga_node_size = 4,
  paga_edge_threshold = 0.01,
  paga_edge_size = c(0.2, 1),
  paga_edge_color = "grey40",
  paga_edge_alpha = 0.5,
  paga_transition_threshold = 0.01,
  paga_transition_size = c(0.2, 1),
  paga_transition_color = "black",
  paga_transition_alpha = 1,
  paga_show_transition = FALSE,
  velocity = NULL,
  velocity_plot_type = "raw",
  velocity_n_neighbors = ceiling(ncol(srt@assays[[1]])/50),
  velocity_density = 1,
  velocity_smooth = 0.5,
  velocity_scale = 1,
  velocity_min_mass = 1,
  velocity_cutoff_perc = 5,
  velocity_arrow_color = "black",
  velocity_arrow_angle = 20,
  streamline_L = 5,
  streamline_minL = 1,
  streamline_res = 1,
  streamline_n = 15,
  streamline_width = c(0, 0.8),
  streamline_alpha = 1,
  streamline_color = NULL,
  streamline_palette = "RdYlBu",
  streamline_palcolor = NULL,
  streamline_bg_color = "white",
  streamline_bg_stroke = 0.5,
  hex = FALSE,
  hex.linewidth = 0.5,
  hex.count = TRUE,
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

- group.by:

  Name of one or more meta.data columns to group (color) cells by (for
  example, orig.ident).

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- dims:

  Dimensions to plot, must be a two-length numeric vector specifying x-
  and y-dimensions

- split.by:

  Name of a column in meta.data column to split plot by.

- cells:

  Subset cells to plot.

- show_na:

  Whether to assign a color from the color palette to NA group. If
  `TRUE`, cell points with NA level will be colored by `bg_color`. If
  `FALSE`, cell points with NA level will be removed from the plot.

- show_stat:

  Whether to show statistical information on the plot.

- pt.size:

  Point size.

- pt.alpha:

  Point transparency.

- palette:

  Name of a color palette name collected in scop. Default is `"Paired"`.

- palcolor:

  Custom colors used to create a color palette.

- bg_color:

  Color value for background(NA) points.

- label:

  Whether to label the cell groups.

- label.size:

  Size of labels.

- label.fg:

  Foreground color of label.

- label.bg:

  Background color of label.

- label.bg.r:

  Background ratio of label.

- label_insitu:

  Whether to place the raw labels (group names) in the center of the
  cells with the corresponding group. Default is `FALSE`, which using
  numbers instead of raw labels.

- label_repel:

  Logical value indicating whether the label is repel away from the
  center points.

- label_repulsion:

  Force of repulsion between overlapping text labels. Default is `20`.

- label_point_size:

  Size of the center points.

- label_point_color:

  Color of the center points.

- label_segment_color:

  Color of the line segment for labels.

- cells.highlight:

  A vector of cell names to highlight.

- cols.highlight:

  Color used to highlight the cells.

- sizes.highlight:

  Size of highlighted cell points.

- alpha.highlight:

  Transparency of highlighted cell points.

- stroke.highlight:

  Border width of highlighted cell points.

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

- add_mark:

  Whether to add marks around cell groups. Default is `FALSE`.

- mark_type:

  Type of mark to add around cell groups. One of "hull", "ellipse",
  "rect", or "circle". Default is `"hull"`.

- mark_expand:

  Expansion of the mark around the cell group. Default is
  `grid::unit(3, "mm")`.

- mark_alpha:

  Transparency of the mark. Default is `0.1`.

- mark_linetype:

  Line type of the mark border. Default is `1` (solid line).

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

- stat.by:

  The name of a metadata column to stat.

- stat_type:

  Set stat types ("percent" or "count").

- stat_plot_type:

  Set the statistical plot type.

- stat_plot_position:

  Position adjustment in statistical plot.

- stat_plot_size:

  Set the statistical plot size. Default is `0.1`.

- stat_plot_palette:

  Color palette used in statistical plot.

- stat_palcolor:

  Custom colors used in statistical plot

- stat_plot_alpha:

  Transparency of the statistical plot.

- stat_plot_label:

  Whether to add labels in the statistical plot.

- stat_plot_label_size:

  Label size in the statistical plot.

- graph:

  Specify the graph name to add edges between cell neighbors to the
  plot.

- edge_size:

  Size of edges.

- edge_alpha:

  Transparency of edges.

- edge_color:

  Color of edges.

- paga:

  Specify the calculated paga results to add a PAGA graph layer to the
  plot.

- paga_type:

  PAGA plot type. "connectivities" or "connectivities_tree".

- paga_node_size:

  Size of the nodes in PAGA plot.

- paga_edge_threshold:

  Threshold of edge connectivities in PAGA plot.

- paga_edge_size:

  Size of edges in PAGA plot.

- paga_edge_color:

  Color of edges in PAGA plot.

- paga_edge_alpha:

  Transparency of edges in PAGA plot.

- paga_transition_threshold:

  Threshold of transition edges in PAGA plot.

- paga_transition_size:

  Size of transition edges in PAGA plot.

- paga_transition_color:

  Color of transition edges in PAGA plot.

- paga_transition_alpha:

  Transparency of transition edges in PAGA plot.

- paga_show_transition:

  Whether to show transitions between edges.

- velocity:

  Specify the calculated RNA velocity mode to add a velocity layer to
  the plot.

- velocity_plot_type:

  Set the velocity plot type.

- velocity_n_neighbors:

  Set the number of neighbors used in velocity plot.

- velocity_density:

  Set the density value used in velocity plot.

- velocity_smooth:

  Set the smooth value used in velocity plot.

- velocity_scale:

  Set the scale value used in velocity plot.

- velocity_min_mass:

  Set the min_mass value used in velocity plot.

- velocity_cutoff_perc:

  Set the cutoff_perc value used in velocity plot.

- velocity_arrow_color:

  Color of arrows in velocity plot.

- velocity_arrow_angle:

  Angle of arrows in velocity plot.

- streamline_L:

  Typical length of a streamline in x and y units

- streamline_minL:

  Minimum length of segments to show.

- streamline_res:

  Resolution parameter (higher numbers increases the resolution).

- streamline_n:

  Number of points to draw.

- streamline_width:

  Size of streamline.

- streamline_alpha:

  Transparency of streamline.

- streamline_color:

  Color of streamline.

- streamline_palette:

  Color palette used for streamline.

- streamline_palcolor:

  Custom colors used for streamline.

- streamline_bg_color:

  Background color of streamline.

- streamline_bg_stroke:

  Border width of streamline background.

- hex:

  Whether to chane the plot type from point to the hexagonal bin.

- hex.linewidth:

  Border width of hexagonal bins.

- hex.count:

  Whether show cell counts in each hexagonal bin.

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

  Whether to force drawing regardless of maximum levels in any cell
  group is greater than 100.

- seed:

  Random seed set for reproducibility

## See also

[CellDimPlot3D](https://mengxu98.github.io/scop/reference/CellDimPlot3D.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md),
[FeatureDimPlot3D](https://mengxu98.github.io/scop/reference/FeatureDimPlot3D.md)

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
p1 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
p1
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


thisplot::panel_fix(
  p1,
  height = 2,
  raster = TRUE,
  dpi = 30
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  theme_use = "theme_blank"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  theme_use = ggplot2::theme_classic,
  theme_args = list(base_size = 16)
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


# Highlight cells
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  cells.highlight = colnames(
    pancreas_sub
  )[pancreas_sub$SubCellType == "Epsilon"]
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  split.by = "Phase",
  reduction = "UMAP",
  cells.highlight = TRUE,
  theme_use = "theme_blank",
  legend.position = "none"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


# Add group labels
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label.fg = "orange",
  label.bg = "red",
  label.size = 5
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE,
  label_segment_color = "red"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


# Add various shape of marks
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_expand = grid::unit(1, "mm")
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_alpha = 0.3
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_linetype = 2
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_type = "ellipse"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_type = "rect"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_mark = TRUE,
  mark_type = "circle"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


# Add a density layer
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_density = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_density = TRUE,
  density_filled = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: Removed 396 rows containing missing values or values outside the scale range
#> (`geom_raster()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  add_density = TRUE,
  density_filled = TRUE,
  density_filled_palette = "Blues",
  cells.highlight = TRUE
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: Removed 396 rows containing missing values or values outside the scale range
#> (`geom_raster()`).


# Add statistical charts
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  stat.by = "Phase"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  stat.by = "Phase",
  stat_plot_type = "ring",
  stat_plot_label = TRUE,
  stat_plot_size = 0.15
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text_repel()`).
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  stat.by = "Phase",
  stat_plot_type = "bar",
  stat_type = "count",
  stat_plot_position = "dodge"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


# Chane the plot type from point to the hexagonal bin
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE
)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_hex()`).


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE,
  hex.bins = 20
)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_hex()`).


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE,
  hex.count = FALSE
)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_hex()`).


# Show neighbors graphs on the plot
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  graph = "Standardpca_SNN"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  graph = "Standardpca_SNN",
  edge_color = "grey80"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


# Show lineages based on the pseudotime
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  show_plot = FALSE
)
#> Error in loadNamespace(x): there is no package called ‘slingshot’

FeatureDimPlot(
  pancreas_sub,
  features = paste0("Lineage", 1:2),
  reduction = "UMAP"
)
#> Error in FeatureDimPlot(pancreas_sub, features = paste0("Lineage", 1:2),     reduction = "UMAP"): There are no valid features present.

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2)
)
#> Error in CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP",     lineages = paste0("Lineage", 1:2)): Lineage "Lineage1" is not in the meta.data of srt object

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_whiskers = TRUE
)
#> Error in CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP",     lineages = paste0("Lineage", 1:2), lineages_whiskers = TRUE): Lineage "Lineage1" is not in the meta.data of srt object

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1
)
#> Error in CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP",     lineages = paste0("Lineage", 1:2), lineages_span = 0.1): Lineage "Lineage1" is not in the meta.data of srt object

if (FALSE) { # \dontrun{
# Show PAGA results on the plot
pancreas_sub <- RunPAGA(
  pancreas_sub,
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  return_seurat = TRUE
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@misc$paga
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@misc$paga,
  paga_type = "connectivities_tree"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  label = TRUE,
  label_repel = TRUE,
  label_insitu = TRUE,
  label_segment_color = "transparent",
  paga = pancreas_sub@misc$paga,
  paga_edge_threshold = 0.1,
  paga_edge_color = "black",
  paga_edge_alpha = 1,
  legend.position = "none",
  theme_use = "theme_blank"
)

# Show RNA velocity results on the plot
pancreas_sub <- RunSCVELO(
  pancreas_sub,
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  mode = "stochastic",
  return_seurat = TRUE
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@misc$paga,
  paga_show_transition = TRUE
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = NA,
  velocity = "stochastic"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  velocity = "stochastic",
  velocity_plot_type = "grid"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  velocity = "stochastic",
  velocity_plot_type = "grid",
  velocity_scale = 1.5
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  velocity = "stochastic",
  velocity_plot_type = "stream"
)

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 5,
  pt.alpha = 0.2,
  label = TRUE,
  label_insitu = TRUE,
  velocity = "stochastic",
  velocity_plot_type = "stream",
  velocity_arrow_color = "yellow",
  velocity_density = 2,
  velocity_smooth = 1,
  streamline_n = 20,
  streamline_color = "black",
  legend.position = "none",
  theme_use = "theme_blank"
)
} # }
```
