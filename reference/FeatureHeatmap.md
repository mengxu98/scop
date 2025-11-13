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
  Default is an empty list.

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

[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 11:55:40] Start standard scop workflow...
#> ℹ [2025-11-13 11:55:40] Checking a list of <Seurat> object...
#> ! [2025-11-13 11:55:41] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 11:55:41] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 11:55:42] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 11:55:43] Use the separate HVF from srt_list
#> ℹ [2025-11-13 11:55:43] Number of available HVF: 2000
#> ℹ [2025-11-13 11:55:43] Finished check
#> ℹ [2025-11-13 11:55:43] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 11:55:44] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 11:55:45] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 11:55:45] Reorder clusters...
#> ℹ [2025-11-13 11:55:45] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 11:55:45] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 11:55:45] UMAP will return its model
#> ℹ [2025-11-13 11:55:48] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 11:55:48] UMAP will return its model
#> ✔ [2025-11-13 11:55:52] Run scop standard workflow done
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group_by = "CellType"
)
#> ✔ [2025-11-13 11:55:52] immunogenomics/presto installed successfully
#> ℹ [2025-11-13 11:55:53] Data type is log-normalized
#> ℹ [2025-11-13 11:55:53] Start differential expression test
#> ℹ [2025-11-13 11:55:53] Find all markers(wilcox) among 5 groups...
#> ℹ [2025-11-13 11:55:53] Using 1 core
#> ⠙ [2025-11-13 11:55:53] Running [1/5] ETA:  1s
#> ✔ [2025-11-13 11:55:53] Completed 5 tasks in 962ms
#> 
#> ℹ [2025-11-13 11:55:53] Building results
#> ✔ [2025-11-13 11:55:54] Differential expression test completed
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
#> ℹ [2025-11-13 11:56:07] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2025-11-13 11:56:07] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2025-11-13 11:56:07] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

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
#> ℹ [2025-11-13 11:56:35] Start Enrichment analysis
#> ✔ [2025-11-13 11:56:35] clusterProfiler installed successfully
#> ℹ [2025-11-13 11:56:35] Species: "Mus_musculus"
#> ℹ [2025-11-13 11:56:35] Loading cached: GO_BP version: 3.22.0 nterm:15169 created: 2025-11-13 11:51:23
#> ℹ [2025-11-13 11:56:37] Permform enrichment...
#> ℹ [2025-11-13 11:56:37] Using 1 core
#> ⠙ [2025-11-13 11:56:37] Running [1/5] ETA:  1m
#> ⠹ [2025-11-13 11:56:37] Running [2/5] ETA: 50s
#> ⠸ [2025-11-13 11:56:37] Running [3/5] ETA: 33s
#> ⠼ [2025-11-13 11:56:37] Running [4/5] ETA: 16s
#> ✔ [2025-11-13 11:56:37] Completed 5 tasks in 1m 14.5s
#> 
#> ℹ [2025-11-13 11:56:37] Building results
#> ✔ [2025-11-13 11:57:52] Enrichment analysis done
#> ✔ [2025-11-13 11:57:52] simplifyEnrichment installed successfully
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> ℹ [2025-11-13 11:58:30] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2025-11-13 11:58:30] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2025-11-13 11:58:30] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht3$plot


pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

ht4 <- FeatureHeatmap(
  pancreas_sub,
  features = de_filter$gene,
  nlabel = 10,
  cell_order = names(sort(pancreas_sub$Lineage1)),
  cell_annotation = c("SubCellType", "Lineage1"),
  cell_annotation_palette = c("Paired", "cividis")
)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
ht4$plot


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
