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
  highly variable features (HVF) will be used.

- features_type:

  A character vector specifying the type of features to be used for the
  KNN prediction. Must be one of "HVF" (highly variable features) or
  "DE" (differentially expressed features). Default is `"HVF"`.

- feature_source:

  The source of the features to be used. Must be one of "both", "query",
  or "ref". Default is `"both"`.

- nfeatures:

  The maximum number of features to include in the heatmap. Default is
  `2000`.

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

  The distance metric to use for calculating the pairwise distances
  between cells. Options include: "pearson", "spearman", "cosine",
  "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman",
  "simple matching", and "faith". Additional distance metrics can also
  be used, such as "euclidean", "manhattan", "hamming", etc.

- k:

  The number of nearest neighbors to use for calculating similarities.
  Default is `30`.

- filter_lowfreq:

  The minimum frequency threshold for selecting query dataset features.
  Features with a frequency below this threshold will be excluded from
  the heatmap. Default is `0`.

- prefix:

  The prefix to use for the KNNPredict tool layer in the query dataset.
  This can be used to avoid conflicts with other tools in the Seurat
  object. Default is `"KNNPredict"`.

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
  clutter in large heatmaps. Default is `0`.

- label_cutoff:

  The similarity cutoff for showing labels. Only cells with similarity
  values above this cutoff will have labels. Default is `0`.

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

  Whether to use the `magick` package for rendering rasters. If set to
  `TRUE`, the `magick` package will be used instead of the raster
  package. This can be useful for rendering large heatmaps more
  efficiently. The `magick` package will automatically be installed if
  it is not installed.

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

  Random seed for reproducibility. Default is `11`.

- ht_params:

  Additional parameters to customize the appearance of the heatmap. This
  should be a list with named elements, where the names correspond to
  parameter names in the
  [ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  function. Any conflicting parameters will override the defaults set by
  this function. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A list with the following elements:

- `plot`: The heatmap plot as a ggplot object.

- `features`: The features used in the heatmap.

- `simil_matrix`: The similarity matrix used to generate the heatmap.

- `simil_name`: The name of the similarity metric used to generate the
  heatmap.

- `cell_metadata`: The cell metadata used to generate the heatmap.

## See also

[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md),
[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Computation failed in `stat_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Computation failed in `stat_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Computation failed in `stat_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Computation failed in `stat_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: Groups with fewer than two datapoints have been dropped.
#> ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
#> Warning: Computation failed in `stat_ydensity()`.
#> Caused by error in `$<-.data.frame`:
#> ! replacement has 1 row, data has 0
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
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
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.
ht4$plot
```
