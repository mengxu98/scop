# The Cell Correlation Heatmap

This function generates a heatmap to visualize the similarity between
different cell types or conditions. It takes in Seurat objects or
expression matrices as input and calculates pairwise similarities or
distance.

## Usage

``` r
CellCorHeatmap(
  srt_query,
  srt_ref = NULL,
  bulk_ref = NULL,
  query_group = NULL,
  ref_group = NULL,
  query_assay = NULL,
  ref_assay = NULL,
  query_reduction = NULL,
  ref_reduction = NULL,
  query_dims = 1:30,
  ref_dims = 1:30,
  query_collapsing = !is.null(query_group),
  ref_collapsing = TRUE,
  features = NULL,
  features_type = c("HVF", "DE"),
  feature_source = "both",
  nfeatures = 2000,
  DEtest_param = list(max.cells.per.ident = 200, test.use = "wilcox"),
  DE_threshold = "p_val_adj < 0.05",
  distance_metric = "cosine",
  k = 30,
  filter_lowfreq = 0,
  prefix = "KNNPredict",
  exp_legend_title = NULL,
  border = TRUE,
  flip = FALSE,
  limits = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  row_names_rot = 0,
  column_names_rot = 90,
  row_title = NULL,
  column_title = NULL,
  row_title_side = "left",
  column_title_side = "top",
  row_title_rot = 90,
  column_title_rot = 0,
  nlabel = 0,
  label_cutoff = 0,
  label_by = "row",
  label_size = 10,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  query_group_palette = "Paired",
  query_group_palcolor = NULL,
  ref_group_palette = "simspec",
  ref_group_palcolor = NULL,
  query_annotation = NULL,
  query_annotation_palette = "Paired",
  query_annotation_palcolor = NULL,
  query_cell_annotation_params = if (flip) {
     list(height = grid::unit(10, "mm"))
 }
    else {
     list(width = grid::unit(10, "mm"))
 },
  ref_annotation = NULL,
  ref_annotation_palette = "Paired",
  ref_annotation_palcolor = NULL,
  ref_cell_annotation_params = if (flip) {
     list(width = grid::unit(10, "mm"))
 }
    else {
     list(height = grid::unit(10, "mm"))
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

- srt_query:

  A Seurat object or count matrix representing the query dataset. This
  dataset will be used to calculate the similarities between cells.

- srt_ref:

  A Seurat object or count matrix representing the reference dataset. If
  provided, the similarities will be calculated between cells from the
  query and reference datasets. If not provided, the similarities will
  be calculated within the query dataset.

- bulk_ref:

  A count matrix representing bulk data. If provided, the similarities
  will be calculated between cells from the query dataset and bulk data.

- query_group:

  The grouping variable in the query dataset. This variable will be used
  to group cells in the heatmap rows. If not provided, all cells will be
  treated as one group.

- ref_group:

  The grouping variable in the reference dataset. This variable will be
  used to group cells in the heatmap columns. If not provided, all cells
  will be treated as one group.

- query_assay:

  The assay to use for the query dataset. If not provided, the default
  assay of the query dataset will be used.

- ref_assay:

  The assay to use for the reference dataset. If not provided, the
  default assay of the reference dataset will be used.

- query_reduction:

  The dimensionality reduction method to use for the query dataset. If
  not provided, no dimensionality reduction will be applied to the query
  dataset.

- ref_reduction:

  The dimensionality reduction method to use for the reference dataset.
  If not provided, no dimensionality reduction will be applied to the
  reference dataset.

- query_dims:

  The dimensions to use for the query dataset. If not provided, the
  first 30 dimensions will be used.

- ref_dims:

  The dimensions to use for the reference dataset. If not provided, the
  first 30 dimensions will be used.

- query_collapsing:

  Whether to collapse cells within each query group before calculating
  similarities. If set to TRUE, the similarities will be calculated
  between query groups rather than individual cells.

- ref_collapsing:

  Detail description of the `query_collapsing` argument.

- features:

  A vector of feature names to include in the heatmap. If not provided,
  a default set of highly variable features (HVF) will be used.

- features_type:

  The type of features to use. Options are `"HVF"` for highly variable
  features, `"DE"` for differentially expressed features between query
  and reference groups.

- feature_source:

  The source of features to use. Options are `"query"` to use only
  features from the query dataset, `"ref"` to use only features from the
  reference dataset, or `"both"` to use features from both datasets. If
  not provided or set to "both", features will be selected from both
  datasets.

- nfeatures:

  The maximum number of features to include in the heatmap. Default is
  2000.

- DEtest_param:

  The parameters to use for differential expression testing. This should
  be a list with two elements: `"max.cells.per.ident"` specifying the
  maximum number of cells per group for differential expression testing,
  and `"test.use"` specifying the statistical test to use for
  differential expression testing. Default parameters will be used.

- DE_threshold:

  The threshold for differential expression. Only features with adjusted
  p-values below this threshold will be considered differentially
  expressed.

- distance_metric:

  The distance metric to use for calculating similarities between cells.
  This can be any of the following: `"cosine"`, `"pearson"`,
  `"spearman"`, `"correlation"`, `"jaccard"`, `"ejaccard"`, `"dice"`,
  `"edice"`, `"hamman"`, `"simple matching"`, or `"faith"`. Dhe default
  is `"cosine"`.

- k:

  The number of nearest neighbors to use for calculating similarities.
  Default is 30.

- filter_lowfreq:

  The minimum frequency threshold for selecting query dataset features.
  Features with a frequency below this threshold will be excluded from
  the heatmap. Default is 0.

- prefix:

  The prefix to use for the KNNPredict tool layer in the query dataset.
  This can be used to avoid conflicts with other tools in the Seurat
  object. The default is `"KNNPredict"`.

- exp_legend_title:

  The title for the color legend in the heatmap. If not provided, a
  default title based on the similarity metric will be used.

- border:

  Whether to add a border around each heatmap cell. The default is
  `TRUE`.

- flip:

  Whether to flip the orientation of the heatmap. If set to TRUE, the
  rows and columns of the heatmap will be swapped. This can be useful
  for visualizing large datasets in a more compact form. The default is
  `FALSE`.

- limits:

  The limits for the color scale in the heatmap. If not provided, the
  default is to use the range of similarity values.

- cluster_rows:

  Whether to cluster the rows of the heatmap. If set to TRUE, the rows
  will be rearranged based on hierarchical clustering. The default is
  `FALSE`.

- cluster_columns:

  Whether to cluster the columns of the heatmap. If set to TRUE, the
  columns will be rearranged based on hierarchical clustering. The
  default is `FALSE`.

- show_row_names:

  Whether to show the row names in the heatmap. The default is `FALSE`.

- show_column_names:

  Whether to show the column names in the heatmap. The default is
  `FALSE`.

- row_names_side:

  The side of the heatmap to show the row names. Options are `"left"` or
  `"right"`. Default is `"left"`.

- column_names_side:

  The side of the heatmap to show the column names. Options are `"top"`
  or `"bottom"`. If not provided, Default is `"top"`.

- row_names_rot:

  The rotation angle of the row names. If not provided, Default is `0`
  degrees.

- column_names_rot:

  The rotation angle of the column names. If not provided, the default
  is 90 degrees.

- row_title:

  The title for the row names in the heatmap. If not provided, the
  default is to use the query grouping variable.

- column_title:

  The title for the column names in the heatmap. Default is to use the
  reference grouping variable.

- row_title_side:

  The side of the heatmap to show the row title. Options are `"top"` or
  `"bottom"`. Default is `"left"`.

- column_title_side:

  The side of the heatmap to show the column title. Options are `"left"`
  or `"right"`. Default is `"top"`.

- row_title_rot:

  The rotation angle of the row title. Default is `90` degrees.

- column_title_rot:

  The rotation angle of the column title. Default is `0` degrees.

- nlabel:

  The maximum number of labels to show on each side of the heatmap. If
  set to 0, no labels will be shown. This can be useful for reducing
  clutter in large heatmaps. Default is 0.

- label_cutoff:

  The similarity cutoff for showing labels. Only cells with similarity
  values above this cutoff will have labels. Default is 0.

- label_by:

  The dimension to use for labeling cells. Options are `"row"` to label
  cells by row, `"column"` to label cells by column, or `"both"` to
  label cells by both row and column. Default is `"row"`.

- label_size:

  The size of the labels. Default is `10`.

- heatmap_palette:

  The color palette to use for the heatmap. This can be any of the
  palettes available in the circlize package. Default is `"RdBu"`.

- heatmap_palcolor:

  The specific colors to use for the heatmap palette. This should be a
  vector of color names or RGB values. Default is `NULL`.

- query_group_palette:

  The color palette to use for the query group legend. This can be any
  of the palettes available in the circlize package. Default is
  `"Paired"`.

- query_group_palcolor:

  The specific colors to use for the query group palette. This should be
  a vector of color names or RGB values. Default is `NULL`.

- ref_group_palette:

  The color palette to use for the reference group legend. This can be
  any of the palettes available in the circlize package. Default is
  `"simspec"`.

- ref_group_palcolor:

  The specific colors to use for the reference group palette. This
  should be a vector of color names or RGB values. Default is `NULL`.

- query_annotation:

  A vector of cell metadata column names or assay feature names to use
  for highlighting specific cells in the heatmap. Each element of the
  vector will create a separate cell annotation track in the heatmap. If
  not provided, no cell annotations will be shown.

- query_annotation_palette:

  The color palette to use for the query cell annotation tracks. This
  can be any of the palettes available in the circlize package. If a
  single color palette is provided, it will be used for all cell
  annotation tracks. If multiple color palettes are provided, each track
  will be assigned a separate palette. Default is `"Paired"`.

- query_annotation_palcolor:

  The specific colors to use for the query cell annotation palettes.
  This should be a list of vectors, where each vector contains the
  colors for a specific cell annotation track. If a single color vector
  is provided, it will be used for all cell annotation tracks. Default
  is `NULL`.

- query_cell_annotation_params:

  Additional parameters to customize the appearance of the query cell
  annotation tracks. This should be a list with named elements, where
  the names correspond to parameter names in the
  [ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  function. Any conflicting parameters will override the defaults set by
  this function.

- ref_annotation:

  A vector of cell metadata column names or assay feature names to use
  for highlighting specific cells in the heatmap. Each element of the
  vector will create a separate cell annotation track in the heatmap. If
  not provided, no cell annotations will be shown.

- ref_annotation_palette:

  The color palette to use for the reference cell annotation tracks.
  This can be any of the palettes available in the circlize package. If
  a single color palette is provided, it will be used for all cell
  annotation tracks. If multiple color palettes are provided, each track
  will be assigned a separate palette. Default is `"Paired"`.

- ref_annotation_palcolor:

  The specific colors to use for the reference cell annotation palettes.
  This should be a list of vectors, where each vector contains the
  colors for a specific cell annotation track. If a single color vector
  is provided, it will be used for all cell annotation tracks. If
  multiple color vectors are provided, each track will be assigned a
  separate color vector. Default is `NULL`.

- ref_cell_annotation_params:

  Detail description of the `query_cell_annotation_params` argument.

- use_raster:

  Whether to use raster images for rendering the heatmap. If set to
  `TRUE`, the heatmap will be rendered as a raster image using the
  raster_device argument. Default is determined based on the number of
  rows and columns in the heatmap.

- raster_device:

  The raster device to use for rendering the heatmap. This should be a
  character string specifying the device name, such as `"png"`,
  `"jpeg"`, or `"pdf"`. Default is `"png"`.

- raster_by_magick:

  Whether to use the magick package for rendering rasters. If set to
  `TRUE`, the magick package will be used instead of the raster package.
  This can be useful for rendering large heatmaps more efficiently. If
  the magick package is not installed, this argument will be ignored.

- height:

  The height of the heatmap in the specified units. If not provided, the
  height will be automatically determined based on the number of rows in
  the heatmap and the default unit.

- width:

  The width of the heatmap in the specified units. If not provided, the
  width will be automatically determined based on the number of columns
  in the heatmap and the default unit.

- units:

  The units to use for the width and height of the heatmap. Default is
  `"inch"`, Options are `"mm"`, `"cm"`, or `"inch"`.

- seed:

  The random seed to use for reproducible results. Default is `11`.

- ht_params:

  Additional parameters to customize the appearance of the heatmap. This
  should be a list with named elements, where the names correspond to
  parameter names in the
  [ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  function. Any conflicting parameters will override the defaults set by
  this function.

## Value

A list with the following elements:

- `plot:` The heatmap plot as a ggplot object.

- `features:` The features used in the heatmap.

- `simil_matrix:` The similarity matrix used to generate the heatmap.

- `simil_name:` The name of the similarity metric used to generate the
  heatmap.

- `cell_metadata:` The cell metadata used to generate the heatmap.

## See also

[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md),
[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md)

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
ht1 <- CellCorHeatmap(
  srt_query = pancreas_sub,
  query_group = "SubCellType"
)
#> As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
#> This message is displayed once per session.
ht1$plot


data(panc8_sub)
# Simply convert genes from human to mouse and preprocess the data
genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames
)
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]

ht2 <- CellCorHeatmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  nlabel = 3,
  label_cutoff = 0.6,
  query_group = "SubCellType",
  ref_group = "celltype",
  query_annotation = "Phase",
  query_annotation_palette = "Set2",
  ref_annotation = "tech",
  ref_annotation_palette = "Set3"
)
ht2$plot


ht3 <- CellCorHeatmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "SubCellType",
  query_collapsing = FALSE,
  cluster_rows = TRUE,
  ref_group = "celltype",
  ref_collapsing = FALSE,
  cluster_columns = TRUE
)
ht3$plot


ht4 <- CellCorHeatmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  show_row_names = TRUE,
  show_column_names = TRUE,
  query_group = "SubCellType",
  ref_group = "celltype",
  query_annotation = c(
    "Sox9", "Rbp4", "Gcg", "Nap1l2", "Xist"
  ),
  ref_annotation = c(
    "Sox9", "Rbp4", "Gcg", "Nap1l2", "Xist"
  )
)
ht4$plot
```
