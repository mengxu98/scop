# Feature Heatmap

Feature Heatmap

## Usage

``` r
FeatureHeatmap(
  srt,
  features = NULL,
  cells = NULL,
  group.by = NULL,
  split.by = NULL,
  within_groups = FALSE,
  max_cells = 100,
  cell_order = NULL,
  border = TRUE,
  flip = FALSE,
  layer = "counts",
  assay = NULL,
  exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"),
  exp_legend_title = NULL,
  limits = NULL,
  lib_normalize = identical(layer, "counts"),
  libsize = NULL,
  feature_split = NULL,
  feature_split_by = NULL,
  n_split = NULL,
  split_order = NULL,
  split_method = c("kmeans", "hclust", "mfuzz"),
  decreasing = FALSE,
  fuzzification = NULL,
  cluster_features_by = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
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
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  group_palette = "Paired",
  group_palcolor = NULL,
  cell_split_palette = "simspec",
  cell_split_palcolor = NULL,
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

- features:

  The features to include in the heatmap. Default is `NULL`.

- cells:

  A character vector specifying the cells to include in the heatmap.
  Default is `NULL`.

- group.by:

  A character vector specifying the groups to group by. Default is
  `NULL`.

- split.by:

  A character vector specifying the variable to split the heatmap by.
  Default is `NULL`.

- within_groups:

  Whether to create separate heatmap scales for each group or within
  each group. Default is `FALSE`.

- max_cells:

  An integer, maximum number of cells to sample per group. Default is
  `100`.

- cell_order:

  A vector of cell names defining the order of cells. Default is `NULL`.

- border:

  Whether to add a border to the heatmap. Default is `TRUE`.

- flip:

  Whether to flip the heatmap. Default is `FALSE`.

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

- cluster_features_by:

  A character vector specifying which group.by to use when clustering
  features. Default is `NULL`. By default, this parameter is set to
  NULL, which means that all groups will be used.

- cluster_rows:

  Whether to cluster rows in the heatmap. Default is `FALSE`.

- cluster_columns:

  Whether to cluster columns in the heatmap. Default is `FALSE`.

- cluster_row_slices:

  Whether to cluster row slices in the heatmap. Default is `FALSE`.

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

- heatmap_palette:

  A character vector specifying the palette to use for the heatmap.
  Default is `"RdBu"`.

- heatmap_palcolor:

  A character vector specifying the heatmap color to use. Default is
  `NULL`.

- group_palette:

  A character vector specifying the palette to use for groups. Default
  is `"Paired"`.

- group_palcolor:

  A character vector specifying the group color to use. Default is
  `NULL`.

- cell_split_palette:

  A character vector specifying the palette to use for cell splits.
  Default is `"simspec"`.

- cell_split_palcolor:

  A character vector specifying the cell split color to use. Default is
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
  Default is [`list()`](https://rdrr.io/r/base/list.html).

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
  function. Default is [`list()`](https://rdrr.io/r/base/list.html).

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
#> Warning: The `slot` argument of `Assays()` is deprecated as of SeuratObject 5.0.0.
#> ℹ Please use `LayerData()` instead.
#> ℹ The deprecated feature was likely used in the scop package.
#>   Please report the issue at <https://github.com/mengxu98/scop/issues>.
#> ⠙ [2025-12-05 08:30:54] Running [1/5] Processing: Ductal  ETA:  1s
#> ✔ [2025-12-05 08:30:54] Completed 5 tasks in 871ms
#> 
de_filter <- dplyr::filter(
  pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
ht1 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  group.by = "CellType",
  split.by = "Phase",
  cell_split_palette = "Dark2"
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
ht1$plot


thisplot::panel_fix(
  ht1$plot,
  height = 4,
  width = 6,
  raster = TRUE,
  dpi = 50
)


ht2 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  group.by = c("CellType", "SubCellType"),
  n_split = 4,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  ht_params = list(row_gap = grid::unit(0, "mm")),
  use_raster = FALSE
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.

ht2$plot


ht3 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  feature_split = de_filter$group1,
  group.by = "CellType",
  species = "Mus_musculus",
  db = "GO_BP",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#>  
#> → Will install 51 packages.
#> → All 51 packages (0 B) are cached.
#> + AnnotationDbi       1.72.0  [bld]
#> + AnnotationHub       4.0.0   [bld][cmp]
#> + BiocBaseUtils       1.12.0  [bld]
#> + BiocFileCache       3.0.0   [bld]
#> + BiocVersion         3.22.0  [bld]
#> + Biostrings          2.78.0  [bld][cmp]
#> + DBI                 1.2.3   
#> + DOSE                4.4.0   [bld]
#> + GO.db               3.22.0  [bld]
#> + GOSemSim            2.36.0  [bld][cmp]
#> + KEGGREST            1.50.0  [bld]
#> + R.methodsS3         1.8.2   
#> + R.oo                1.27.1  
#> + R.utils             2.13.0  
#> + RSQLite             2.4.5   
#> + ape                 5.8-1   
#> + bit                 4.6.0   
#> + bit64               4.6.0-1 
#> + blob                1.2.4   
#> + clipr               0.8.0    + ✔ libx11-dev
#> + clusterProfiler     4.18.2  [bld]
#> + dbplyr              2.5.1   
#> + enrichplot          1.30.4  [bld]
#> + fastmatch           1.1-6   
#> + fgsea               1.36.0  [bld][cmp]
#> + filelock            1.0.3   
#> + fontBitstreamVera   0.1.1   
#> + fontLiberation      0.1.0   
#> + fontquiver          0.2.1   
#> + gdtools             0.4.4    + ✔ libcairo2-dev, ✔ libfontconfig1-dev, ✔ libfreetype6-dev
#> + ggforce             0.5.0   
#> + ggiraph             0.9.2    + ✔ libpng-dev
#> + ggnewscale          0.5.2   
#> + ggtangle            0.0.9   
#> + ggtree              4.0.1   [bld]
#> + gson                0.1.0   
#> + hms                 1.1.4   
#> + httr2               1.2.1   
#> + org.Hs.eg.db        3.22.0  [bld]
#> + quarto              1.5.1   
#> + qvalue              2.42.0  [bld]
#> + readr               2.1.6   
#> + rstudioapi          0.17.1  
#> + scatterpie          0.2.6   
#> + systemfonts         1.3.1    + ✔ libfontconfig1-dev, ✔ libfreetype6-dev
#> + tidydr              0.0.6   
#> + tidytree            0.4.6   
#> + treeio              1.34.0  [bld]
#> + tweenr              2.0.3   
#> + tzdb                0.5.0   
#> + vroom               1.6.7   
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 51 pkgs are cached
#> ✔ Got BiocVersion 3.22.0 (source) (1.11 kB)
#> ✔ Got blob 1.2.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (47.47 kB)
#> ✔ Got BiocBaseUtils 1.12.0 (source) (232.56 kB)
#> ✔ Got clipr 0.8.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (51.24 kB)
#> ✔ Got filelock 1.0.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (24.70 kB)
#> ✔ Got bit64 4.6.0-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (492.54 kB)
#> ✔ Got BiocFileCache 3.0.0 (source) (744.19 kB)
#> ✔ Got bit 4.6.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (628.10 kB)
#> ✔ Got ggtangle 0.0.9 (x86_64-pc-linux-gnu-ubuntu-24.04) (257.09 kB)
#> ✔ Got clusterProfiler 4.18.2 (source) (632.38 kB)
#> ✔ Got AnnotationHub 4.0.0 (source) (1.00 MB)
#> ✔ Got scatterpie 0.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (150.90 kB)
#> ✔ Got dbplyr 2.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.24 MB)
#> ✔ Got RSQLite 2.4.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.35 MB)
#> ✔ Got tzdb 0.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (666.60 kB)
#> ✔ Got GOSemSim 2.36.0 (source) (610.99 kB)
#> ✔ Got fontBitstreamVera 0.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (699.47 kB)
#> ✔ Got hms 1.1.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (103.38 kB)
#> ✔ Got fastmatch 1.1-6 (x86_64-pc-linux-gnu-ubuntu-24.04) (35.95 kB)
#> ✔ Got gdtools 0.4.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (203.79 kB)
#> ✔ Got rstudioapi 0.17.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (317.69 kB)
#> ✔ Got ggtree 4.0.1 (source) (370.24 kB)
#> ✔ Got gson 0.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (232.54 kB)
#> ✔ Got KEGGREST 1.50.0 (source) (239.73 kB)
#> ✔ Got qvalue 2.42.0 (source) (2.77 MB)
#> ✔ Got ape 5.8-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.97 MB)
#> ✔ Got vroom 1.6.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (966.05 kB)
#> ✔ Got fontquiver 0.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.28 MB)
#> ✔ Got R.methodsS3 1.8.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (82.67 kB)
#> ✔ Got tweenr 2.0.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (461.68 kB)
#> ✔ Got treeio 1.34.0 (source) (701.64 kB)
#> ✔ Got readr 2.1.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (868.61 kB)
#> ✔ Got ggiraph 0.9.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.67 MB)
#> ✔ Got ggforce 0.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.95 MB)
#> ✔ Got tidytree 0.4.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (343.02 kB)
#> ✔ Got R.oo 1.27.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (996.09 kB)
#> ✔ Got quarto 1.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (544.37 kB)
#> ✔ Got ggnewscale 0.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (351.15 kB)
#> ✔ Got DBI 1.2.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (916.93 kB)
#> ✔ Got systemfonts 1.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (816.45 kB)
#> ✔ Got AnnotationDbi 1.72.0 (source) (4.38 MB)
#> ✔ Got httr2 1.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (782.87 kB)
#> ✔ Got DOSE 4.4.0 (source) (5.75 MB)
#> ✔ Got R.utils 2.13.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.45 MB)
#> ✔ Got fgsea 1.36.0 (source) (6.17 MB)
#> ✔ Got fontLiberation 0.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.54 MB)
#> ✔ Got Biostrings 2.78.0 (source) (12.82 MB)
#> ✔ Got GO.db 3.22.0 (source) (25.25 MB)
#> ✔ Got org.Hs.eg.db 3.22.0 (source) (104.97 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install libx11-dev libcairo2-dev libfontconfig1-dev libfreetype6-dev libpng-dev libcurl4-openssl-dev libssl-dev make libglpk-dev libxml2-dev pandoc libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libx11-dev is already the newest version (2:1.8.7-1build1).
#> libx11-dev set to manually installed.
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
#> 0 upgraded, 0 newly installed, 0 to remove and 49 not upgraded.
#> ℹ Building BiocBaseUtils 1.12.0
#> ℹ Building BiocVersion 3.22.0
#> ℹ Building Biostrings 2.78.0
#> ℹ Building qvalue 2.42.0
#> ✔ Built BiocVersion 3.22.0 (1.5s)
#> ✔ Installed ape 5.8-1  (103ms)
#> ✔ Installed bit 4.6.0  (42ms)
#> ✔ Installed bit64 4.6.0-1  (1.1s)
#> ✔ Built BiocBaseUtils 1.12.0 (3.2s)
#> ✔ Installed blob 1.2.4  (117ms)
#> ✔ Installed clipr 0.8.0  (163ms)
#> ✔ Installed DBI 1.2.3  (135ms)
#> ✔ Installed dbplyr 2.5.1  (134ms)
#> ✔ Installed fastmatch 1.1-6  (113ms)
#> ℹ Building fgsea 1.36.0
#> ✔ Installed filelock 1.0.3  (203ms)
#> ✔ Installed fontBitstreamVera 0.1.1  (53ms)
#> ✔ Installed fontLiberation 0.1.0  (94ms)
#> ✔ Installed fontquiver 0.2.1  (1.1s)
#> ✔ Installed gdtools 0.4.4  (121ms)
#> ✔ Built qvalue 2.42.0 (5.4s)
#> ✔ Installed ggforce 0.5.0  (110ms)
#> ✔ Installed ggiraph 0.9.2  (131ms)
#> ✔ Installed ggnewscale 0.5.2  (118ms)
#> ✔ Installed ggtangle 0.0.9  (114ms)
#> ✔ Installed gson 0.1.0  (126ms)
#> ✔ Installed hms 1.1.4  (146ms)
#> ✔ Installed httr2 1.2.1  (232ms)
#> ✔ Installed quarto 1.5.1  (149ms)
#> ✔ Installed R.methodsS3 1.8.2  (138ms)
#> ✔ Installed R.oo 1.27.1  (137ms)
#> ✔ Installed R.utils 2.13.0  (130ms)
#> ✔ Installed readr 2.1.6  (139ms)
#> ✔ Installed RSQLite 2.4.5  (161ms)
#> ℹ Building BiocFileCache 3.0.0
#> ✔ Installed rstudioapi 0.17.1  (231ms)
#> ✔ Installed scatterpie 0.2.6  (59ms)
#> ✔ Installed systemfonts 1.3.1  (55ms)
#> ✔ Installed tidydr 0.0.6  (62ms)
#> ✔ Installed tidytree 0.4.6  (36ms)
#> ℹ Building treeio 1.34.0
#> ✔ Built BiocFileCache 3.0.0 (8.4s)
#> ✔ Installed tweenr 2.0.3  (75ms)
#> ✔ Installed tzdb 0.5.0  (102ms)
#> ✔ Installed vroom 1.6.7  (120ms)
#> ✔ Installed BiocBaseUtils 1.12.0  (91ms)
#> ✔ Installed BiocFileCache 3.0.0  (58ms)
#> ✔ Installed BiocVersion 3.22.0  (51ms)
#> ✔ Installed qvalue 2.42.0  (125ms)
#> ✔ Built treeio 1.34.0 (9.4s)
#> ✔ Installed treeio 1.34.0  (63ms)
#> ℹ Building ggtree 4.0.1
#> ✔ Built ggtree 4.0.1 (12s)
#> ✔ Installed ggtree 4.0.1  (111ms)
#> ✔ Built Biostrings 2.78.0 (32.4s)
#> ✔ Installed Biostrings 2.78.0  (170ms)
#> ℹ Building KEGGREST 1.50.0
#> ✔ Built fgsea 1.36.0 (33.1s)
#> ✔ Installed fgsea 1.36.0  (132ms)
#> ✔ Built KEGGREST 1.50.0 (5s)
#> ✔ Installed KEGGREST 1.50.0  (26ms)
#> ℹ Building AnnotationDbi 1.72.0
#> ✔ Built AnnotationDbi 1.72.0 (11.6s)
#> ✔ Installed AnnotationDbi 1.72.0  (70ms)
#> ℹ Building AnnotationHub 4.0.0
#> ℹ Building GO.db 3.22.0
#> ℹ Building org.Hs.eg.db 3.22.0
#> ✔ Built AnnotationHub 4.0.0 (12.3s)
#> ✔ Installed AnnotationHub 4.0.0  (58ms)
#> ✔ Built GO.db 3.22.0 (37.4s)
#> ✔ Installed GO.db 3.22.0  (490ms)
#> ℹ Building GOSemSim 2.36.0
#> ✔ Built GOSemSim 2.36.0 (13.2s)
#> ✔ Installed GOSemSim 2.36.0  (41ms)
#> ℹ Building DOSE 4.4.0
#> ✔ Built DOSE 4.4.0 (11.6s)
#> ✔ Installed DOSE 4.4.0  (55ms)
#> ℹ Building enrichplot 1.30.4
#> ✔ Built enrichplot 1.30.4 (11.5s)
#> ✔ Installed enrichplot 1.30.4  (23ms)
#> ℹ Building clusterProfiler 4.18.2
#> ✔ Built clusterProfiler 4.18.2 (12.3s)
#> ✔ Installed clusterProfiler 4.18.2  (31ms)
#> ✔ Built org.Hs.eg.db 3.22.0 (3m 56.2s)
#> ✔ Installed org.Hs.eg.db 3.22.0  (2.1s)
#> ✔ 1 pkg + 155 deps: kept 105, added 51, dld 49 (195.15 MB) [5m 13.9s]
#> Error in loadNamespace(x): there is no package called ‘R.cache’
ht3$plot
#> Error: object 'ht3' not found

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Error in loadNamespace(x): there is no package called ‘slingshot’
ht4 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  nlabel = 10,
  cell_order = names(sort(pancreas_sub$Lineage1)),
  cell_annotation = c("SubCellType", "Lineage1"),
  cell_annotation_palette = c("Paired", "cividis")
)
#> Error in FeatureHeatmap(pancreas_sub, features = de_filter$gene, nlabel = 10,     cell_order = names(sort(pancreas_sub$Lineage1)), cell_annotation = c("SubCellType",         "Lineage1"), cell_annotation_palette = c("Paired", "cividis")): Cell_annotation: Lineage1 is not in the Seurat object.
ht4$plot
#> Error: object 'ht4' not found

if (FALSE) { # \dontrun{
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = c("CSPA", "TF")
)

ht5 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  n_split = 4,
  group.by = "CellType",
  heatmap_palette = "viridis",
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"), c("forestgreen")
  ),
  cell_annotation = c("Phase", "G2M_score"),
  cell_annotation_palette = c("Dark2", "Purples")
)
ht5$plot

ht6 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  n_split = 4,
  group.by = "CellType",
  heatmap_palette = "viridis",
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"), c("forestgreen")
  ),
  cell_annotation = c("Phase", "G2M_score"),
  cell_annotation_palette = c("Dark2", "Purples"),
  flip = TRUE,
  column_title_rot = 45
)
ht6$plot
} # }
```
