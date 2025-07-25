% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DynamicHeatmap.R
\name{DynamicHeatmap}
\alias{DynamicHeatmap}
\title{Heatmap plot for dynamic features along lineages}
\usage{
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
  Ensembl_version = 103,
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
}
\arguments{
\item{srt}{A Seurat object.}

\item{lineages}{A character vector specifying the lineages to plot.}

\item{features}{A character vector specifying the features to plot.
By default, this parameter is set to NULL, and the dynamic features will be determined by the parameters
\code{min_expcells}, \code{r.sq}, \code{dev.expl}, \code{padjust} and \code{num_intersections}.}

\item{use_fitted}{A logical indicating whether to use fitted values. Default is FALSE.}

\item{border}{A logical indicating whether to add a border to the heatmap. Default is TRUE.}

\item{flip}{A logical indicating whether to flip the heatmap. Default is FALSE.}

\item{min_expcells}{A numeric value specifying the minimum number of expected cells. Default is 20.}

\item{r.sq}{A numeric value specifying the R-squared threshold. Default is 0.2.}

\item{dev.expl}{A numeric value specifying the deviance explained threshold. Default is 0.2.}

\item{padjust}{A numeric value specifying the p-value adjustment threshold. Default is 0.05.}

\item{num_intersections}{This parameter is a numeric vector used to determine the number of intersections among lineages. It helps in selecting which dynamic features will be used. By default, when this parameter is set to NULL, all dynamic features that pass the specified threshold will be used for each lineage.}

\item{cell_density}{A numeric value is used to define the cell density within each cell bin. By default, this parameter is set to 1, which means that all cells will be included within each cell bin.}

\item{cell_bins}{A numeric value specifying the number of cell bins. Default is 100.}

\item{order_by}{A character vector specifying the order of the heatmap. Default is "peaktime".}

\item{layer}{A character vector specifying the layer in the Seurat object to use.
Default is "counts".}

\item{assay}{A character vector specifying the assay in the Seurat object to use.
Default is NULL.}

\item{exp_method}{A character vector specifying the method for calculating expression values.
Default is "zscore" with options "zscore", "raw", "fc", "log2fc", "log1p".}

\item{exp_legend_title}{A character vector specifying the title for the legend of expression value.
Default is NULL.}

\item{limits}{A two-length numeric vector specifying the limits for the color scale.
Default is NULL.}

\item{lib_normalize}{A logical value indicating whether to normalize the data by library size.}

\item{libsize}{A numeric vector specifying the library size for each cell.
Default is NULL.}

\item{family}{A character specifying the model used to calculate the dynamic features if needed. By default, this parameter is set to NULL, and the appropriate family will be automatically determined.}

\item{cluster_features_by}{A character vector specifying which lineage to use when clustering features. By default, this parameter is set to NULL, which means that all lineages will be used.}

\item{cluster_rows}{A logical value indicating whether to cluster rows in the heatmap.
Default is FALSE.}

\item{cluster_row_slices}{A logical value indicating whether to cluster row slices in the heatmap.
Default is FALSE.}

\item{cluster_columns}{A logical value indicating whether to cluster columns in the heatmap.
Default is FALSE.}

\item{cluster_column_slices}{A logical value indicating whether to cluster column slices in the heatmap.
Default is FALSE.}

\item{show_row_names}{A logical value indicating whether to show row names in the heatmap.
Default is FALSE.}

\item{show_column_names}{A logical value indicating whether to show column names in the heatmap.
Default is FALSE.}

\item{row_names_side}{A character vector specifying the side to place row names.}

\item{column_names_side}{A character vector specifying the side to place column names.}

\item{row_names_rot}{A numeric value specifying the rotation angle for row names.
Default is 0.}

\item{column_names_rot}{A numeric value specifying the rotation angle for column names.
Default is 90.}

\item{row_title}{A character vector specifying the title for rows.
Default is NULL.}

\item{column_title}{A character vector specifying the title for columns.
Default is NULL.}

\item{row_title_side}{A character vector specifying the side to place row title.
Default is "left".}

\item{column_title_side}{A character vector specifying the side to place column title.
Default is "top".}

\item{row_title_rot}{A numeric value specifying the rotation angle for row title.
Default is 0.}

\item{column_title_rot}{A numeric value specifying the rotation angle for column title.}

\item{feature_split}{A factor specifying how to split the features.
Default is NULL.}

\item{feature_split_by}{A character vector specifying which group.by to use when splitting features (into n_split feature clusters).
Default is NULL.}

\item{n_split}{An integer specifying the number of feature splits (feature clusters) to create.
Default is NULL.}

\item{split_order}{A numeric vector specifying the order of splits. Default is NULL.}

\item{split_method}{A character vector specifying the method for splitting features.
Default is "kmeans" with options "kmeans", "hclust", "mfuzz").}

\item{decreasing}{A logical value indicating whether to sort feature splits in decreasing order.
Default is FALSE.}

\item{fuzzification}{A numeric value specifying the fuzzification coefficient.
Default is NULL.}

\item{anno_terms}{A logical value indicating whether to include term annotations.
Default is FALSE.}

\item{anno_keys}{A logical value indicating whether to include key annotations.
Default is FALSE.}

\item{anno_features}{A logical value indicating whether to include feature annotations.
Default is FALSE.}

\item{terms_width}{A unit specifying the width of term annotations.
Default is unit(4, "in").}

\item{terms_fontsize}{A numeric vector specifying the font size(s) for term annotations.
Default is 8.}

\item{keys_width}{A unit specifying the width of key annotations.
Default is unit(2, "in").}

\item{keys_fontsize}{A two-length numeric vector specifying the minimum and maximum font size(s) for key annotations.
Default is c(6, 10).}

\item{features_width}{A unit specifying the width of feature annotations.
Default is unit(2, "in").}

\item{features_fontsize}{A two-length numeric vector specifying the minimum and maximum font size(s) for feature annotations.
Default is c(6, 10).}

\item{IDtype}{A character vector specifying the type of IDs for features.
Default is "symbol".}

\item{species}{A character vector specifying the species for features.
Default is "Homo_sapiens".}

\item{db_update}{A logical value indicating whether to update the database.
Default is FALSE.}

\item{db_version}{A character vector specifying the version of the database.
Default is "latest".}

\item{db_combine}{A logical value indicating whether to use a combined database.
Default is FALSE.}

\item{convert_species}{A logical value indicating whether to use a species-converted database if annotation is missing for \code{species}.
Default is FALSE.}

\item{Ensembl_version}{An integer specifying the Ensembl version.
Default is 103.}

\item{mirror}{A character vector specifying the mirror for the Ensembl database.
Default is NULL.}

\item{db}{A character vector specifying the database to use.
Default is "GO_BP".}

\item{TERM2GENE}{A data.frame specifying the TERM2GENE mapping for the database.
Default is NULL.}

\item{TERM2NAME}{A data.frame specifying the TERM2NAME mapping for the database.
Default is NULL.}

\item{minGSSize}{An integer specifying the minimum gene set size for the database.
Default is 10.}

\item{maxGSSize}{An integer specifying the maximum gene set size for the database.
Default is 500.}

\item{GO_simplify}{A logical value indicating whether to simplify gene ontology terms.
Default is FALSE.}

\item{GO_simplify_cutoff}{A character vector specifying the cutoff for GO simplification.
Default is "p.adjust < 0.05".}

\item{simplify_method}{A character vector specifying the method for GO simplification.
Default is "Wang".}

\item{simplify_similarityCutoff}{A numeric value specifying the similarity cutoff for GO simplification.
Default is 0.7.}

\item{pvalueCutoff}{A numeric vector specifying the p-value cutoff(s) for significance.
Default is NULL.}

\item{padjustCutoff}{A numeric value specifying the adjusted p-value cutoff for significance.
Default is 0.05.}

\item{topTerm}{An integer specifying the number of top terms to include.
Default is 5.}

\item{show_termid}{A logical value indicating whether to show term IDs.
Default is FALSE.}

\item{topWord}{An integer specifying the number of top words to include.
Default is 20.}

\item{words_excluded}{A character vector specifying the words to exclude.
Default is NULL.}

\item{nlabel}{An integer specifying the number of labels to include.
Default is 0.}

\item{features_label}{A character vector specifying the features to label.
Default is NULL.}

\item{label_size}{A numeric value specifying the size of labels.
Default is 10.}

\item{label_color}{A character vector specifying the color of labels.
Default is "black".}

\item{pseudotime_label}{A numeric vector specifying the pseudotime label. Default is NULL.}

\item{pseudotime_label_color}{A character string specifying the pseudotime label color. Default is "black".}

\item{pseudotime_label_linetype}{A numeric value specifying the pseudotime label line type. Default is 2.}

\item{pseudotime_label_linewidth}{A numeric value specifying the pseudotime label line width. Default is 3.}

\item{heatmap_palette}{A character vector specifying the palette to use for the heatmap.
Default is "RdBu".}

\item{heatmap_palcolor}{A character vector specifying the heatmap color to use.
Default is NULL.}

\item{pseudotime_palette}{A character vector specifying the color palette to use for pseudotime.}

\item{pseudotime_palcolor}{A list specifying the colors to use for the pseudotime in the heatmap.}

\item{feature_split_palette}{A character vector specifying the palette to use for feature splits.
Default is "simspec".}

\item{feature_split_palcolor}{A character vector specifying the feature split color to use.
Default is NULL.}

\item{cell_annotation}{A character vector specifying the cell annotation(s) to include.
Default is NULL.}

\item{cell_annotation_palette}{A character vector specifying the palette to use for cell annotations.
The length of the vector should match the number of cell_annotation. Default is "Paired".}

\item{cell_annotation_palcolor}{A list of character vector specifying the cell annotation color(s) to use.
The length of the list should match the number of cell_annotation. Default is NULL.}

\item{cell_annotation_params}{A list specifying additional parameters for cell annotations.
Default is a list with width = unit(1, "cm") if flip is TRUE, else a list with height = unit(1, "cm").}

\item{feature_annotation}{A character vector specifying the feature annotation(s) to include.
Default is NULL.}

\item{feature_annotation_palette}{A character vector specifying the palette to use for feature annotations.
The length of the vector should match the number of feature_annotation.
Default is "Dark2".}

\item{feature_annotation_palcolor}{A list of character vector specifying the feature annotation color to use.
The length of the list should match the number of feature_annotation.
Default is NULL.}

\item{feature_annotation_params}{A list specifying additional parameters for feature annotations.
Default is an empty list.}

\item{separate_annotation}{A character vector of names of annotations to be displayed in separate annotation blocks. Each name should match a column name in the metadata of the Seurat object.}

\item{separate_annotation_palette}{A character vector specifying the color palette to use for separate annotations.}

\item{separate_annotation_palcolor}{A list specifying the colors to use for each level of the separate annotations.}

\item{separate_annotation_params}{A list of other parameters to be passed to the HeatmapAnnotation function when creating the separate annotation blocks.}

\item{reverse_ht}{A logical indicating whether to reverse the heatmap. Default is NULL.}

\item{use_raster}{A logical value indicating whether to use a raster device for plotting.
Default is NULL.}

\item{raster_device}{A character vector specifying the raster device to use.
Default is "png".}

\item{raster_by_magick}{A logical value indicating whether to use the 'magick' package for raster.
Default is FALSE.}

\item{height}{A numeric vector specifying the height(s) of the heatmap body.
Default is NULL.}

\item{width}{A numeric vector specifying the width(s) of the heatmap body.
Default is NULL.}

\item{units}{A character vector specifying the units for the height and width.
Default is "inch".}

\item{seed}{An integer specifying the random seed. Default is 11.}

\item{ht_params}{A list specifying additional parameters passed to the \link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap} function.
Default is an empty list.}
}
\description{
Heatmap plot for dynamic features along lineages
}
\examples{
data(pancreas_sub)
# pancreas_sub <- AnnotateFeatures(
#   srt = pancreas_sub,
#   species = "Mus_musculus",
#   db = c("TF", "CSPA")
# )

pancreas_sub <- RunSlingshot(
  srt = pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
pancreas_sub <- RunDynamicFeatures(
  srt = pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200
)

ht1 <- DynamicHeatmap(
  srt = pancreas_sub,
  lineages = "Lineage1",
  n_split = 5,
  split_method = "kmeans-peaktime",
  cell_annotation = "SubCellType"
)
ht1$plot

panel_fix(ht1$plot, raster = TRUE, dpi = 50)

ht2 <- DynamicHeatmap(
  srt = pancreas_sub,
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
ht2$plot

panel_fix(
  ht2$plot,
  height = 5,
  width = 5,
  raster = TRUE,
  dpi = 50
)

ht3 <- DynamicHeatmap(
  srt = pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_split = 5,
  split_method = "kmeans",
  cluster_rows = TRUE,
  cell_annotation = "SubCellType"
)
ht3$plot

ht4 <- DynamicHeatmap(
  srt = pancreas_sub,
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
ht4$plot

ht5 <- DynamicHeatmap(
  srt = pancreas_sub,
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
ht5$plot

ht6 <- DynamicHeatmap(
  srt = pancreas_sub,
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
ht6$plot
}
\seealso{
\code{\link{RunDynamicFeatures}} \code{\link{RunDynamicEnrichment}}
}
