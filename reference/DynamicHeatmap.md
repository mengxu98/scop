# Heatmap plot for dynamic features along lineages

Heatmap plot for dynamic features along lineages

## Usage

``` r
DynamicHeatmap(
  srt,
  lineages,
  features = NULL,
  use_fitted = FALSE,
  border = TRUE,
  flip = FALSE,
  min_expcells = 20,
  r.sq = 0.2,
  dev.expl = 0.2,
  padjust = 0.05,
  num_intersections = NULL,
  cell_density = 1,
  cell_bins = 100,
  order_by = c("peaktime", "valleytime"),
  layer = "counts",
  assay = NULL,
  exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"),
  exp_legend_title = NULL,
  limits = NULL,
  lib_normalize = identical(layer, "counts"),
  libsize = NULL,
  family = NULL,
  cluster_features_by = NULL,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_side = ifelse(flip, "left", "right"),
  column_names_side = ifelse(flip, "bottom", "top"),
  row_names_rot = 0,
  column_names_rot = 90,
  row_title = NULL,
  column_title = NULL,
  row_title_side = "left",
  column_title_side = "top",
  row_title_rot = 0,
  column_title_rot = ifelse(flip, 90, 0),
  feature_split = NULL,
  feature_split_by = NULL,
  n_split = NULL,
  split_order = NULL,
  split_method = c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime"),
  decreasing = FALSE,
  fuzzification = NULL,
  anno_terms = FALSE,
  anno_keys = FALSE,
  anno_features = FALSE,
  terms_width = grid::unit(4, "in"),
  terms_fontsize = 8,
  keys_width = grid::unit(2, "in"),
  keys_fontsize = c(6, 10),
  features_width = grid::unit(2, "in"),
  features_fontsize = c(6, 10),
  IDtype = "symbol",
  species = "Homo_sapiens",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = FALSE,
  Ensembl_version = NULL,
  mirror = NULL,
  db = "GO_BP",
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  GO_simplify = FALSE,
  GO_simplify_cutoff = "p.adjust < 0.05",
  simplify_method = "Wang",
  simplify_similarityCutoff = 0.7,
  pvalueCutoff = NULL,
  padjustCutoff = 0.05,
  topTerm = 5,
  show_termid = FALSE,
  topWord = 20,
  words_excluded = NULL,
  nlabel = 20,
  features_label = NULL,
  label_size = 10,
  label_color = "black",
  pseudotime_label = NULL,
  pseudotime_label_color = "black",
  pseudotime_label_linetype = 2,
  pseudotime_label_linewidth = 3,
  heatmap_palette = "viridis",
  heatmap_palcolor = NULL,
  pseudotime_palette = "cividis",
  pseudotime_palcolor = NULL,
  feature_split_palette = "simspec",
  feature_split_palcolor = NULL,
  cell_annotation = NULL,
  cell_annotation_palette = "Paired",
  cell_annotation_palcolor = NULL,
  cell_annotation_params = if (flip) {
     list(width = grid::unit(5, "mm"))
 } else {
 
       list(height = grid::unit(5, "mm"))
 },
  feature_annotation = NULL,
  feature_annotation_palette = "Dark2",
  feature_annotation_palcolor = NULL,
  feature_annotation_params = if (flip) {
     list(height = grid::unit(5, "mm"))
 } else
    {
     list(width = grid::unit(5, "mm"))
 },
  separate_annotation = NULL,
  separate_annotation_palette = "Paired",
  separate_annotation_palcolor = NULL,
  separate_annotation_params = if (flip) {
     list(width = grid::unit(10, "mm"))
 }
    else {
     list(height = grid::unit(10, "mm"))
 },
  reverse_ht = NULL,
  use_raster = NULL,
  raster_device = "png",
  raster_by_magick = FALSE,
  height = NULL,
  width = NULL,
  units = "inch",
  seed = 11,
  ht_params = list()
)
```

## Arguments

- srt:

  A Seurat object.

- lineages:

  A character vector specifying the lineages to plot.

- features:

  A character vector specifying the features to plot. By default, this
  parameter is set to NULL, and the dynamic features will be determined
  by the parameters `min_expcells`, `r.sq`, `dev.expl`, `padjust` and
  `num_intersections`.

- use_fitted:

  Whether to use fitted values. Default is `FALSE`.

- border:

  Whether to add a border to the heatmap. Default is `TRUE`.

- flip:

  Whether to flip the heatmap. Default is `FALSE`.

- min_expcells:

  The minimum number of expected cells. Default is `20`.

- r.sq:

  The R-squared threshold. Default is `0.2`.

- dev.expl:

  The deviance explained threshold. Default is `0.2`.

- padjust:

  The p-value adjustment threshold. Default is `0.05`.

- num_intersections:

  This parameter is a numeric vector used to determine the number of
  intersections among lineages. It helps in selecting which dynamic
  features will be used. By default, when this parameter is set to
  `NULL`, all dynamic features that pass the specified threshold will be
  used for each lineage.

- cell_density:

  The cell density within each cell bin. By default, this parameter is
  set to `1`, which means that all cells will be included within each
  cell bin.

- cell_bins:

  The number of cell bins. Default is `100`.

- order_by:

  The order of the heatmap. Default is `"peaktime"`.

- layer:

  A character vector specifying the layer in the Seurat object to use.
  Default is `"counts"`.

- assay:

  A character vector specifying the assay in the Seurat object to use.
  Default is `NULL`.

- exp_method:

  A character vector specifying the method for calculating expression
  values. Options are `"zscore"`, `"raw"`, `"fc"`, `"log2fc"`, or
  `"log1p"`. Default is `"zscore"`.

- exp_legend_title:

  A character vector specifying the title for the legend of expression
  value. Default is `NULL`.

- limits:

  A two-length numeric vector specifying the limits for the color scale.
  Default is `NULL`.

- lib_normalize:

  Whether to normalize the data by library size.

- libsize:

  A numeric vector specifying the library size for each cell. Default is
  `NULL`.

- family:

  The model used to calculate the dynamic features if needed. By
  default, this parameter is set to `NULL`, and the appropriate family
  will be automatically determined.

- cluster_features_by:

  Which lineage to use when clustering features. By default, this
  parameter is set to `NULL`, which means that all lineages will be
  used.

- cluster_rows:

  Whether to cluster rows in the heatmap. Default is `FALSE`.

- cluster_row_slices:

  Whether to cluster row slices in the heatmap. Default is `FALSE`.

- cluster_columns:

  Whether to cluster columns in the heatmap. Default is `FALSE`.

- cluster_column_slices:

  Whether to cluster column slices in the heatmap. Default is `FALSE`.

- show_row_names:

  Whether to show row names in the heatmap. Default is `FALSE`.

- show_column_names:

  Whether to show column names in the heatmap. Default is `FALSE`.

- row_names_side:

  A character vector specifying the side to place row names.

- column_names_side:

  A character vector specifying the side to place column names.

- row_names_rot:

  The rotation angle for row names. Default is `0`.

- column_names_rot:

  The rotation angle for column names. Default is `90`.

- row_title:

  A character vector specifying the title for rows. Default is `NULL`.

- column_title:

  A character vector specifying the title for columns. Default is
  `NULL`.

- row_title_side:

  A character vector specifying the side to place row title. Default is
  `"left"`.

- column_title_side:

  A character vector specifying the side to place column title. Default
  is `"top"`.

- row_title_rot:

  The rotation angle for row title. Default is `0`.

- column_title_rot:

  The rotation angle for column title.

- feature_split:

  A factor specifying how to split the features. Default is `NULL`.

- feature_split_by:

  A character vector specifying which group.by to use when splitting
  features (into n_split feature clusters). Default is `NULL`.

- n_split:

  A number of feature splits (feature clusters) to create. Default is
  `NULL`.

- split_order:

  A numeric vector specifying the order of splits. Default is `NULL`.

- split_method:

  A character vector specifying the method for splitting features.
  Options are `"kmeans"`, `"hclust"`, or `"mfuzz"`. Default is
  `"kmeans"`.

- decreasing:

  Whether to sort feature splits in decreasing order. Default is
  `FALSE`.

- fuzzification:

  The fuzzification coefficient. Default is `NULL`.

- anno_terms:

  Whether to include term annotations. Default is `FALSE`.

- anno_keys:

  Whether to include key annotations. Default is `FALSE`.

- anno_features:

  Whether to include feature annotations. Default is `FALSE`.

- terms_width:

  A unit specifying the width of term annotations. Default is
  `unit(4, "in")`.

- terms_fontsize:

  A numeric vector specifying the font size(s) for term annotations.
  Default is `8`.

- keys_width:

  A unit specifying the width of key annotations. Default is
  `unit(2, "in")`.

- keys_fontsize:

  A two-length numeric vector specifying the minimum and maximum font
  size(s) for key annotations. Default is `c(6, 10)`.

- features_width:

  A unit specifying the width of feature annotations. Default is
  `unit(2, "in")`.

- features_fontsize:

  A two-length numeric vector specifying the minimum and maximum font
  size(s) for feature annotations. Default is `c(6, 10)`.

- IDtype:

  A character vector specifying the type of IDs for features. Default is
  `"symbol"`.

- species:

  A character vector specifying the species for features. Default is
  `"Homo_sapiens"`.

- db_update:

  Whether to update the database. Default is `FALSE`.

- db_version:

  A character vector specifying the version of the database. Default is
  `"latest"`.

- db_combine:

  Whether to use a combined database. Default is `FALSE`.

- convert_species:

  Whether to use a species-converted database if annotation is missing
  for `species`. Default is `FALSE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `103`.

- mirror:

  A character vector specifying the mirror for the Ensembl database.
  Default is `NULL`.

- db:

  A character vector specifying the database to use. Default is
  `"GO_BP"`.

- TERM2GENE:

  A data.frame specifying the TERM2GENE mapping for the database.
  Default is `NULL`.

- TERM2NAME:

  A data.frame specifying the TERM2NAME mapping for the database.
  Default is `NULL`.

- minGSSize:

  An integer specifying the minimum gene set size for the database.
  Default is `10`.

- maxGSSize:

  An integer specifying the maximum gene set size for the database.
  Default is `500`.

- GO_simplify:

  Whether to simplify gene ontology terms. Default is `FALSE`.

- GO_simplify_cutoff:

  A character vector specifying the cutoff for GO simplification.
  Default is `"p.adjust < 0.05"`.

- simplify_method:

  A character vector specifying the method for GO simplification.
  Default is `"Wang"`.

- simplify_similarityCutoff:

  The similarity cutoff for GO simplification. Default is `0.7`.

- pvalueCutoff:

  A numeric vector specifying the p-value cutoff(s) for significance.
  Default is `NULL`.

- padjustCutoff:

  The adjusted p-value cutoff for significance. Default is `0.05`.

- topTerm:

  A number of top terms to include. Default is `5`.

- show_termid:

  Whether to show term IDs. Default is `FALSE`.

- topWord:

  A number of top words to include. Default is `20`.

- words_excluded:

  A character vector specifying the words to exclude. Default is `NULL`.

- nlabel:

  A number of labels to include. Default is `0`.

- features_label:

  A character vector specifying the features to label. Default is
  `NULL`.

- label_size:

  The size of labels. Default is `10`.

- label_color:

  A character vector specifying the color of labels. Default is
  `"black"`.

- pseudotime_label:

  The pseudotime label. Default is `NULL`.

- pseudotime_label_color:

  The pseudotime label color. Default is `"black"`.

- pseudotime_label_linetype:

  The pseudotime label line type. Default is `2`.

- pseudotime_label_linewidth:

  The pseudotime label line width. Default is `3`.

- heatmap_palette:

  A character vector specifying the palette to use for the heatmap.
  Default is `"RdBu"`.

- heatmap_palcolor:

  A character vector specifying the heatmap color to use. Default is
  `NULL`.

- pseudotime_palette:

  The color palette to use for pseudotime. Default is `"cividis"`.

- pseudotime_palcolor:

  The colors to use for the pseudotime in the heatmap. Default is
  `NULL`.

- feature_split_palette:

  A character vector specifying the palette to use for feature splits.
  Default is `"simspec"`.

- feature_split_palcolor:

  A character vector specifying the feature split color to use. Default
  is `NULL`.

- cell_annotation:

  A character vector specifying the cell annotation(s) to include.
  Default is `NULL`.

- cell_annotation_palette:

  A character vector specifying the palette to use for cell annotations.
  The length of the vector should match the number of cell_annotation.
  Default is `"Paired"`.

- cell_annotation_palcolor:

  A list of character vector specifying the cell annotation color(s) to
  use. The length of the list should match the number of
  cell_annotation. Default is `NULL`.

- cell_annotation_params:

  A list specifying additional parameters for cell annotations. Default
  is a list with `width = unit(1, "cm")` if flip is TRUE, else a list
  with `height = unit(1, "cm")`.

- feature_annotation:

  A character vector specifying the feature annotation(s) to include.
  Default is `NULL`.

- feature_annotation_palette:

  A character vector specifying the palette to use for feature
  annotations. The length of the vector should match the number of
  feature_annotation. Default is `"Dark2"`.

- feature_annotation_palcolor:

  A list of character vector specifying the feature annotation color to
  use. The length of the list should match the number of
  feature_annotation. Default is `NULL`.

- feature_annotation_params:

  A list specifying additional parameters for feature annotations.
  Default is an empty list.

- separate_annotation:

  Names of the annotations to be displayed in separate annotation
  blocks. Each name should match a column name in the metadata of the
  `Seurat` object.

- separate_annotation_palette:

  The color palette to use for separate annotations. Default is
  `"Paired"`.

- separate_annotation_palcolor:

  The colors to use for each level of the separate annotations. Default
  is `NULL`.

- separate_annotation_params:

  Other parameters to
  [ComplexHeatmap::HeatmapAnnotation](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)
  when creating a separate annotation blocks.

- reverse_ht:

  Whether to reverse the heatmap. Default is `NULL`.

- use_raster:

  Whether to use a raster device for plotting. Default is `NULL`.

- raster_device:

  A character vector specifying the raster device to use. Default is
  `"png"`.

- raster_by_magick:

  Whether to use the 'magick' package for raster. Default is `FALSE`.

- height:

  A numeric vector specifying the height(s) of the heatmap body. Default
  is `NULL`.

- width:

  A numeric vector specifying the width(s) of the heatmap body. Default
  is `NULL`.

- units:

  A character vector specifying the units for the height and width.
  Default is `"inch"`.

- seed:

  An integer specifying the random seed. Default is `11`.

- ht_params:

  A list specifying additional parameters passed to the
  [ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  function. Default is an empty list.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)

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

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200
)
#>  
#> → Will update 1 package.
#> → The package (0 B) is cached.
#> + mgcv 1.9-3 → 1.9-4 
#>   
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Got mgcv 1.9-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.63 MB)
#> ✔ Installed mgcv 1.9-4  (1s)
#> ✔ 1 pkg + 3 deps: kept 3, upd 1, dld 1 (3.63 MB) [2.6s]
#> ⠙ [2025-12-02 04:05:26] Running [1/225] Processing: 1  ETA: 19s
#> ⠹ [2025-12-02 04:05:26] Running [51/225] Processing: 51  ETA:  6s
#> ⠸ [2025-12-02 04:05:26] Running [146/225] Processing: 146  ETA:  3s
#> ✔ [2025-12-02 04:05:26] Completed 225 tasks in 7.3s
#> 
#> ⠙ [2025-12-02 04:05:34] Running [14/225] Processing: 14  ETA:  5s
#> ⠹ [2025-12-02 04:05:34] Running [115/225] Processing: 115  ETA:  3s
#> ⠸ [2025-12-02 04:05:34] Running [212/225] Processing: 212  ETA:  0s
#> ✔ [2025-12-02 04:05:34] Completed 225 tasks in 6.8s
#> 

ht1 <- DynamicHeatmap(
  pancreas_sub,
  lineages = "Lineage1",
  n_split = 5,
  split_method = "kmeans-peaktime",
  cell_annotation = "SubCellType"
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.

ht1$plot


thisplot::panel_fix(ht1$plot, raster = TRUE, dpi = 50)


ht2 <- DynamicHeatmap(
  pancreas_sub,
  lineages = "Lineage1",
  features = c(
    "Sox9",
    "Neurod2",
    "Isl1",
    "Rbp4",
    "Pyy", "S_score", "G2M_score"
  ),
  cell_annotation = "SubCellType"
)
#> ⠙ [2025-12-02 04:05:45] Running [1/2] Processing: 1  ETA:  0s
#> ✔ [2025-12-02 04:05:45] Completed 2 tasks in 42ms
#> 
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
ht2$plot


thisplot::panel_fix(
  ht2$plot,
  height = 5,
  width = 5,
  raster = TRUE,
  dpi = 50
)


ht3 <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_split = 5,
  split_method = "kmeans",
  cluster_rows = TRUE,
  cell_annotation = "SubCellType"
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.

ht3$plot


ht4 <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  cell_annotation = "SubCellType",
  n_split = 5,
  split_method = "mfuzz",
  species = "Mus_musculus",
  db = "GO_BP",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE
)
#>  
#> → Will install 2 packages.
#> → All 2 packages (0 B) are cached.
#> + e1071   1.7-16 
#> + proxy   0.4-27 
#>   
#> ℹ No downloads are needed, 2 pkgs are cached
#> ✔ Got proxy 0.4-27 (x86_64-pc-linux-gnu-ubuntu-24.04) (175.47 kB)
#> ✔ Got e1071 1.7-16 (x86_64-pc-linux-gnu-ubuntu-24.04) (596.57 kB)
#> ✔ Installed e1071 1.7-16  (31ms)
#> ✔ Installed proxy 0.4-27  (43ms)
#> ✔ 1 pkg + 3 deps: kept 2, added 2, dld 2 (772.04 kB) [1.1s]
#> Registered S3 methods overwritten by 'proxy':
#>   method               from    
#>   print.registry_field registry
#>   print.registry_entry registry
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#>  
#> → Will install 35 packages.
#> → All 35 packages (0 B) are cached.
#> + BH                  1.87.0-1 
#> + BiocParallel        1.44.0   [bld][cmp]
#> + DOSE                4.4.0    [bld]
#> + GO.db               3.22.0   [bld]
#> + GOSemSim            2.36.0   [bld][cmp]
#> + R.methodsS3         1.8.2    
#> + R.oo                1.27.1   
#> + R.utils             2.13.0   
#> + ape                 5.8-1    
#> + clusterProfiler     4.18.2   [bld]
#> + cpp11               0.5.2    
#> + enrichplot          1.30.3   [bld]
#> + fastmatch           1.1-6    
#> + fgsea               1.36.0   [bld][cmp]
#> + fontBitstreamVera   0.1.1    
#> + fontLiberation      0.1.0    
#> + fontquiver          0.2.1    
#> + formatR             1.14     
#> + futile.logger       1.4.3    
#> + futile.options      1.0.1    
#> + gdtools             0.4.4     + ✔ libcairo2-dev, ✔ libfontconfig1-dev, ✔ libfreetype6-dev
#> + ggforce             0.5.0    
#> + ggiraph             0.9.2     + ✔ libpng-dev
#> + ggnewscale          0.5.2    
#> + ggtangle            0.0.9    
#> + ggtree              4.0.1    [bld]
#> + gson                0.1.0    
#> + lambda.r            1.2.4    
#> + polyclip            1.10-7   
#> + qvalue              2.42.0   [bld]
#> + scatterpie          0.2.6    
#> + snow                0.4-4    
#> + tidytree            0.4.6    
#> + treeio              1.34.0   [bld]
#> + tweenr              2.0.3    
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 35 pkgs are cached
#> ✔ Got enrichplot 1.30.3 (source) (100.35 kB)
#> ✔ Got ggtree 4.0.1 (source) (370.24 kB)
#> ✔ Got fastmatch 1.1-6 (x86_64-pc-linux-gnu-ubuntu-24.04) (35.95 kB)
#> ✔ Got GOSemSim 2.36.0 (source) (610.99 kB)
#> ✔ Got cpp11 0.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (289.40 kB)
#> ✔ Got BiocParallel 1.44.0 (source) (1.11 MB)
#> ✔ Got treeio 1.34.0 (source) (701.64 kB)
#> ✔ Got fontBitstreamVera 0.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (699.47 kB)
#> ✔ Got futile.logger 1.4.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (96.77 kB)
#> ✔ Got scatterpie 0.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (150.90 kB)
#> ✔ Got R.utils 2.13.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.45 MB)
#> ✔ Got formatR 1.14 (x86_64-pc-linux-gnu-ubuntu-24.04) (151.65 kB)
#> ✔ Got tweenr 2.0.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (461.68 kB)
#> ✔ Got snow 0.4-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (97.07 kB)
#> ✔ Got ggiraph 0.9.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.67 MB)
#> ✔ Got ggnewscale 0.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (351.15 kB)
#> ✔ Got polyclip 1.10-7 (x86_64-pc-linux-gnu-ubuntu-24.04) (120.02 kB)
#> ✔ Got gdtools 0.4.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (203.79 kB)
#> ✔ Got qvalue 2.42.0 (source) (2.77 MB)
#> ✔ Got gson 0.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (232.54 kB)
#> ✔ Got tidytree 0.4.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (343.02 kB)
#> ✔ Got R.methodsS3 1.8.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (82.67 kB)
#> ✔ Got lambda.r 1.2.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (110.01 kB)
#> ✔ Got futile.options 1.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (18.42 kB)
#> ✔ Got ape 5.8-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.97 MB)
#> ✔ Got fontquiver 0.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.28 MB)
#> ✔ Got R.oo 1.27.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (996.09 kB)
#> ✔ Got fontLiberation 0.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.54 MB)
#> ✔ Got DOSE 4.4.0 (source) (5.75 MB)
#> ✔ Got fgsea 1.36.0 (source) (6.17 MB)
#> ✔ Got ggforce 0.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.95 MB)
#> ✔ Got BH 1.87.0-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (13.35 MB)
#> ✔ Got GO.db 3.22.0 (source) (25.25 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install libcairo2-dev libfontconfig1-dev libfreetype6-dev libpng-dev libcurl4-openssl-dev libssl-dev make libglpk-dev libxml2-dev pandoc libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcairo2-dev is already the newest version (1.18.0-3build1).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 62 not upgraded.
#> ℹ Building qvalue 2.42.0
#> ℹ Building GO.db 3.22.0
#> ✔ Installed ape 5.8-1  (1.1s)
#> ✔ Installed cpp11 0.5.2  (47ms)
#> ✔ Installed fastmatch 1.1-6  (28ms)
#> ✔ Installed BH 1.87.0-1  (1.4s)
#> ✔ Installed fontBitstreamVera 0.1.1  (103ms)
#> ✔ Installed fontLiberation 0.1.0  (111ms)
#> ✔ Installed fontquiver 0.2.1  (159ms)
#> ✔ Installed formatR 1.14  (109ms)
#> ✔ Installed futile.logger 1.4.3  (106ms)
#> ✔ Installed futile.options 1.0.1  (102ms)
#> ✔ Installed gdtools 0.4.4  (102ms)
#> ✔ Installed ggforce 0.5.0  (164ms)
#> ✔ Installed ggiraph 0.9.2  (108ms)
#> ✔ Installed ggnewscale 0.5.2  (172ms)
#> ✔ Installed ggtangle 0.0.9  (112ms)
#> ✔ Installed gson 0.1.0  (107ms)
#> ✔ Installed lambda.r 1.2.4  (104ms)
#> ✔ Installed polyclip 1.10-7  (107ms)
#> ✔ Installed R.methodsS3 1.8.2  (104ms)
#> ✔ Installed R.oo 1.27.1  (106ms)
#> ✔ Installed R.utils 2.13.0  (158ms)
#> ✔ Installed scatterpie 0.2.6  (158ms)
#> ✔ Installed snow 0.4-4  (101ms)
#> ℹ Building BiocParallel 1.44.0
#> ✔ Installed tidytree 0.4.6  (126ms)
#> ℹ Building treeio 1.34.0
#> ✔ Built qvalue 2.42.0 (3.8s)
#> ✔ Installed tweenr 2.0.3  (1s)
#> ✔ Installed qvalue 2.42.0  (53ms)
#> ✔ Built treeio 1.34.0 (6.1s)
#> ✔ Installed treeio 1.34.0  (45ms)
#> ℹ Building ggtree 4.0.1
#> ✔ Built BiocParallel 1.44.0 (13.1s)
#> ✔ Installed BiocParallel 1.44.0  (84ms)
#> ℹ Building fgsea 1.36.0
#> ✔ Built ggtree 4.0.1 (8.4s)
#> ✔ Installed ggtree 4.0.1  (42ms)
#> ✔ Built GO.db 3.22.0 (41.3s)
#> ✔ Built fgsea 1.36.0 (24.9s)
#> ✔ Installed fgsea 1.36.0  (118ms)
#> ✔ Installed GO.db 3.22.0  (1.5s)
#> ℹ Building GOSemSim 2.36.0
#> ✔ Built GOSemSim 2.36.0 (12.3s)
#> ✔ Installed GOSemSim 2.36.0  (1s)
#> ℹ Building DOSE 4.4.0
#> ✔ Built DOSE 4.4.0 (11s)
#> ✔ Installed DOSE 4.4.0  (48ms)
#> ℹ Building enrichplot 1.30.3
#> ✔ Built enrichplot 1.30.3 (10.9s)
#> ✔ Installed enrichplot 1.30.3  (21ms)
#> ℹ Building clusterProfiler 4.18.2
#> ✔ Built clusterProfiler 4.18.2 (11.5s)
#> ✔ Installed clusterProfiler 4.18.2  (1s)
#> ✔ 1 pkg + 125 deps: kept 90, added 35, dld 33 (75.48 MB) [1m 42.3s]
#>  
#> → Will install 1 package.
#> → The package (0 B) is cached.
#> + org.Mm.eg.db   3.22.0 [bld]
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 1 pkg is cached
#> ✔ Got org.Mm.eg.db 3.22.0 (source) (92.82 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install libcurl4-openssl-dev libssl-dev libpng-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> 0 upgraded, 0 newly installed, 0 to remove and 62 not upgraded.
#> ℹ Building org.Mm.eg.db 3.22.0
#> ✔ Built org.Mm.eg.db 3.22.0 (3m 29.7s)
#> ✔ Installed org.Mm.eg.db 3.22.0  (2.9s)
#> ✔ 1 pkg + 34 deps: kept 34, added 1, dld 1 (92.82 MB) [3m 46.8s]
#> Registered S3 method overwritten by 'ggtree':
#>   method         from     
#>   fortify.igraph ggnetwork
#> ⠙ [2025-12-02 04:12:52] Running [1/5] Processing: 1  ETA: 27s
#> ⠹ [2025-12-02 04:12:52] Running [2/5] Processing: 2  ETA: 20s
#> ⠸ [2025-12-02 04:12:52] Running [3/5] Processing: 3  ETA: 12s
#> ⠼ [2025-12-02 04:12:52] Running [4/5] Processing: 4  ETA:  6s
#> ✔ [2025-12-02 04:12:52] Completed 5 tasks in 29.5s
#> 
#>  
#> → Will install 7 packages.
#> → All 7 packages (0 B) are cached.
#> + NLP                  0.3-2  
#> + Polychrome           1.5.4  
#> + scatterplot3d        0.3-44 
#> + simona               1.8.0  [bld][cmp]
#> + simplifyEnrichment   2.4.0  [bld]
#> + slam                 0.1-55 
#> + tm                   0.7-16 
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 7 pkgs are cached
#> ✔ Got slam 0.1-55 (x86_64-pc-linux-gnu-ubuntu-24.04) (187.94 kB)
#> ✔ Got scatterplot3d 0.3-44 (x86_64-pc-linux-gnu-ubuntu-24.04) (348.15 kB)
#> ✔ Got simona 1.8.0 (source) (971.58 kB)
#> ✔ Got simplifyEnrichment 2.4.0 (source) (1.08 MB)
#> ✔ Got NLP 0.3-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (381.71 kB)
#> ✔ Got tm 0.7-16 (x86_64-pc-linux-gnu-ubuntu-24.04) (638.54 kB)
#> ✔ Got Polychrome 1.5.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (684.84 kB)
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
#> ℹ Executing `sudo sh -c apt-get -y install perl libcurl4-openssl-dev libssl-dev make zlib1g-dev libglpk-dev libxml2-dev libpng-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> perl is already the newest version (5.38.2-3.2ubuntu0.2).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> 0 upgraded, 0 newly installed, 0 to remove and 62 not upgraded.
#> ✔ Installed NLP 0.3-2  (59ms)
#> ✔ Installed Polychrome 1.5.4  (75ms)
#> ✔ Installed scatterplot3d 0.3-44  (94ms)
#> ℹ Building simona 1.8.0
#> ✔ Installed slam 0.1-55  (140ms)
#> ✔ Installed tm 0.7-16  (43ms)
#> ✔ Built simona 1.8.0 (47.6s)
#> ✔ Installed simona 1.8.0  (110ms)
#> ℹ Building simplifyEnrichment 2.4.0
#> ✔ Built simplifyEnrichment 2.4.0 (7.7s)
#> ✔ Installed simplifyEnrichment 2.4.0  (1.3s)
#> ✔ 1 pkg + 81 deps: kept 74, added 7, dld 7 (4.29 MB) [1m 1.7s]

ht4$plot


if (FALSE) { # \dontrun{
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = c("CSPA", "TF")
)
ht5 <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  use_fitted = TRUE,
  n_split = 6,
  split_method = "mfuzz",
  heatmap_palette = "viridis",
  cell_annotation = c(
    "SubCellType", "Phase", "G2M_score"
  ),
  cell_annotation_palette = c(
    "Paired", "simspec", "Purples"
  ),
  separate_annotation = list(
    "SubCellType", c("Arxes1", "Ncoa2")
  ),
  separate_annotation_palette = c(
    "Paired", "Set1"
  ),
  separate_annotation_params = list(
    height = grid::unit(10, "mm")
  ),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"),
    c("forestgreen")
  ),
  pseudotime_label = 25,
  pseudotime_label_color = "red"
)
ht5$plot

ht6 <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  use_fitted = TRUE,
  n_split = 6,
  split_method = "mfuzz",
  heatmap_palette = "viridis",
  cell_annotation = c(
    "SubCellType", "Phase", "G2M_score"
  ),
  cell_annotation_palette = c(
    "Paired", "simspec", "Purples"
  ),
  separate_annotation = list(
    "SubCellType", c("Arxes1", "Ncoa2")
  ),
  separate_annotation_palette = c("Paired", "Set1"),
  separate_annotation_params = list(width = grid::unit(10, "mm")),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"),
    c("forestgreen")
  ),
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  flip = TRUE, column_title_rot = 45
)
ht6$plot
} # }
```
