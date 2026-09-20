# Heatmap plot for dynamic features along lineages

Heatmap plot for dynamic features along lineages

## Usage

``` r
DynamicHeatmap(
  object,
  lineages,
  features = NULL,
  use_fitted = FALSE,
  border = TRUE,
  heatmap_border = NULL,
  cell_annotation_border = NULL,
  feature_annotation_border = NULL,
  heatmap_border_palcolor = "black",
  cell_annotation_border_palcolor = "black",
  feature_annotation_border_palcolor = "black",
  heatmap_border_size = 1,
  cell_annotation_border_size = 1,
  feature_annotation_border_size = 1,
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
  suffix = lineages,
  n_candidates = 1000,
  minfreq = 5,
  fit_method = c("gam", "pretsa"),
  knot = 0,
  max_knot_allowed = 10,
  padjust_method = "fdr",
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
  terms_stat_width = grid::unit(1.35, "in"),
  terms_fontsize = 8,
  terms_stat = "none",
  terms_stat_digits = 2,
  terms_stat_label = "value",
  terms_stat_axis = FALSE,
  terms_stat_background_palcolor = NULL,
  terms_stat_border = NULL,
  terms_stat_border_palcolor = NULL,
  terms_stat_border_size = NULL,
  terms_stat_label_palcolor = NULL,
  terms_group_background = FALSE,
  terms_background_palcolor = "grey98",
  terms_background_alpha = 1,
  terms_border = TRUE,
  terms_border_palcolor = "black",
  terms_border_size = 0.8,
  terms_text_palcolor = NULL,
  terms_bar_palcolor = NULL,
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
  cell_annotation_palette = "Chinese",
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
  separate_annotation_palette = "Chinese",
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
  raster_by_magick = TRUE,
  height = NULL,
  width = NULL,
  units = "inch",
  cores = 1,
  verbose = TRUE,
  seed = 11,
  legend.position = "right",
  ht_params = list(),
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- lineages:

  Lineage names for which dynamic features should be calculated.

- features:

  Features to plot. By default, this parameter is set to NULL, and the
  dynamic features will be determined by the parameters `min_expcells`,
  `r.sq`, `dev.expl`, `padjust` and `num_intersections`.

- use_fitted:

  Whether to use fitted values.

- border:

  Draw borders. Kept for compatibility; more specific `*_border`
  arguments inherit this when `NULL`.

- heatmap_border, cell_annotation_border, feature_annotation_border:

  Borders for the heatmap body, annotations, and their matching legends.
  `NULL` inherits `border`.

- heatmap_border_palcolor, cell_annotation_border_palcolor,
  feature_annotation_border_palcolor:

  Border colors for the matching body, annotation, and legend.

- heatmap_border_size, cell_annotation_border_size,
  feature_annotation_border_size:

  Border line widths for the matching body, annotation, and legend.

- flip:

  Flip rows and columns.

- min_expcells:

  The minimum number of expected cells.

- r.sq:

  The R-squared threshold.

- dev.expl:

  The deviance explained threshold.

- padjust:

  The p-value adjustment threshold.

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

  The number of cell bins.

- order_by:

  The order of the heatmap.

- layer:

  Assay layer to use.

- assay:

  Assay to use. `NULL` uses the default assay.

- exp_method:

  Expression transform: `"zscore"`, `"raw"`, `"fc"`, `"log2fc"`, or
  `"log1p"`.

- exp_legend_title:

  Legend title for expression.

- limits:

  Color-scale limits (length 2).

- lib_normalize, libsize:

  Library-size normalization and per-cell library sizes.

- family:

  A character or character vector specifying the family of distributions
  to use for the GAM. If family is set to NULL, the appropriate family
  will be automatically determined based on the data. If length(family)
  is 1, the same family will be used for all features. Otherwise, family
  must have the same length as features.

- suffix:

  Suffix to append to the output layer names for each lineage. Default
  is the lineage names.

- n_candidates:

  A number of candidate features to select when features is `NULL`.

- minfreq:

  Minimum frequency threshold for candidate features. Features with a
  frequency less than minfreq will be excluded.

- fit_method:

  The method used for fitting features. Either `"gam"` (generalized
  additive models) or `"pretsa"` (Pattern recognition in Temporal and
  Spatial Analyses).

- knot:

  For `fit_method = "pretsa"`: B-spline knots. `0` or `"auto"`.

- max_knot_allowed:

  For `fit_method = "pretsa"` when `knot = "auto"`: max knots.

- padjust_method:

  The method used for p-value adjustment.

- cluster_features_by:

  Which lineage to use when clustering features. By default, this
  parameter is set to `NULL`, which means that all lineages will be
  used.

- cluster_rows, cluster_columns, cluster_row_slices,
  cluster_column_slices:

  Heatmap clustering.

- show_row_names:

  Whether to draw row/column names for the heatmap body.

- show_column_names:

  Whether to draw row/column names for the heatmap body.

- row_names_side, column_names_side, row_names_rot, column_names_rot:

  Name placement.

- row_title, column_title, row_title_side, column_title_side,
  row_title_rot, column_title_rot:

  Slice titles.

- feature_split, feature_split_by, n_split, split_order, split_method,
  decreasing:

  Feature splitting. `split_method` is `"kmeans"`, `"hclust"`, or
  `"mfuzz"`.

- fuzzification:

  Mfuzz fuzzification coefficient.

- anno_terms, anno_keys, anno_features:

  Enrichment annotations.

- terms_width, terms_stat_width, terms_fontsize:

  Term annotation size.

- terms_stat:

  Enrichment statistic for term bars: `"none"`, `"score"` (`-log10` of
  the active p-value), or an enrichment column such as `"p.adjust"`.

- terms_stat_digits, terms_stat_label, terms_stat_axis:

  Statistic labels (`"none"`, `"value"`, `"significance"`, `"both"`) and
  shared axis.

- terms_stat_background_palcolor, terms_stat_border,
  terms_stat_border_palcolor, terms_stat_border_size,
  terms_stat_label_palcolor:

  Statistic-panel appearance. `NULL` inherits the matching `terms_*`
  setting.

- terms_group_background, terms_background_palcolor,
  terms_background_alpha, terms_border, terms_border_palcolor,
  terms_border_size, terms_text_palcolor, terms_bar_palcolor:

  Term-block appearance. `terms_text_palcolor = NULL` maps text to
  enrichment significance; `terms_bar_palcolor = NULL` matches bar color
  to term text.

- keys_width, keys_fontsize, features_width, features_fontsize:

  Key and feature annotations.

- IDtype, species, db_combine, mirror, db, TERM2GENE, TERM2NAME,
  minGSSize, maxGSSize:

  Gene-set database (see
  [PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md)).

- db_update:

  Force a refresh. `FALSE` loads the cache when available.

- db_version:

  Database version to retrieve.

- convert_species:

  Use a species-converted database when the annotation is missing for
  `species`.

- Ensembl_version:

  Ensembl version. `NULL` uses the latest.

- GO_simplify, GO_simplify_cutoff, simplify_method,
  simplify_similarityCutoff:

  GO simplification.

- pvalueCutoff, padjustCutoff, topTerm, show_termid, topWord,
  words_excluded:

  Enrichment filters.

- nlabel, features_label, label_size, label_color:

  Feature labels.

- pseudotime_label:

  The pseudotime label.

- pseudotime_label_color:

  The pseudotime label color.

- pseudotime_label_linetype:

  The pseudotime label line type.

- pseudotime_label_linewidth:

  The pseudotime label line width.

- heatmap_palette:

  Palette used for CNV heatmap values.

- heatmap_palcolor:

  Palette used for CNV heatmap values.

- pseudotime_palette:

  The color palette to use for pseudotime.

- pseudotime_palcolor:

  The colors to use for the pseudotime in the heatmap.

- feature_split_palette:

  Split colors.

- feature_split_palcolor:

  Custom colors for feature-split annotation labels.

- cell_annotation, cell_annotation_palette, cell_annotation_palcolor,
  cell_annotation_params:

  Cell annotations. Palette length should match `cell_annotation`.

- feature_annotation, feature_annotation_palette,
  feature_annotation_palcolor, feature_annotation_params:

  Feature annotations. Palette length should match `feature_annotation`.

- separate_annotation:

  Names of the annotations to be displayed in separate annotation
  blocks. Each name should match a column name in the metadata of the
  `Seurat` object.

- separate_annotation_palette:

  The color palette to use for separate annotations.

- separate_annotation_palcolor:

  The colors to use for each level of the separate annotations.

- separate_annotation_params:

  Other parameters to
  [ComplexHeatmap::HeatmapAnnotation](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)
  when creating a separate annotation blocks.

- reverse_ht:

  Whether to reverse the heatmap.

- use_raster, raster_device, raster_by_magick:

  Raster device (`NULL` chooses automatically).

- width, height, units:

  Heatmap size. `NULL` sizes from matrix dimensions.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Optional integer seed. When supplied, every input receives a
  deterministic independent L'Ecuyer-CMRG random-number stream, making
  results reproducible across worker counts and scheduling order. The
  caller's random number state is restored when the call finishes.

- legend.position:

  Legend side (`"right"`, `"left"`, `"top"`, `"bottom"`). Gap to the
  heatmap grows automatically when long row names are on the right.

- ht_params:

  Extra arguments passed to
  [ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
  overriding defaults.

- ...:

  Additional arguments passed to helper functions.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 21:33:31] Start standard processing workflow...
#> ℹ [2026-09-20 21:33:31] Checking a list of <Seurat>...
#> ! [2026-09-20 21:33:32] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 21:33:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:33:32] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 21:33:32] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 21:33:32] Number of available HVF: 2000
#> ℹ [2026-09-20 21:33:32] Finished check
#> ℹ [2026-09-20 21:33:32] Perform `ScaleData()`
#> ℹ [2026-09-20 21:33:32] Perform pca linear dimension reduction
#> ℹ [2026-09-20 21:33:32] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 21:33:32] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 21:33:33] Reorder clusters...
#> ℹ [2026-09-20 21:33:33] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 21:33:33] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 21:33:38] Standard processing workflow completed

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"), ,
  fit_method = "pretsa",
  n_candidates = 200
)
#> ℹ [2026-09-20 21:33:40] Start find dynamic features
#> ℹ [2026-09-20 21:33:41] Data type is raw counts
#> ℹ [2026-09-20 21:33:41] Number of candidate features (union): 242
#> ℹ [2026-09-20 21:33:42] Data type is raw counts
#> ℹ [2026-09-20 21:33:42] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 21:33:42] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-09-20 21:33:42] Find dynamic features done

ht1 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = "Lineage1",
  n_split = 5,
  split_method = "kmeans-peaktime",
  cell_annotation = "SubCellType",
  width = 2,
  height = 3
)
#> ℹ [2026-09-20 21:33:42] [1] 154 features from Lineage1 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
ht1$plot


thisplot::panel_fix(ht1$plot, raster = TRUE, dpi = 50)


ht2 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = "Lineage1",
  features = c(
    "Sox9",
    "Neurod2",
    "Isl1",
    "Rbp4",
    "Pyy",
    "S_score",
    "G2M_score"
  ),
  cell_annotation = "SubCellType"
)
#> ℹ [2026-09-20 21:33:45] Start find dynamic features
#> ℹ [2026-09-20 21:33:45] Data type is raw counts
#> ℹ [2026-09-20 21:33:45] Number of candidate features (union): 2
#> ℹ [2026-09-20 21:33:46] Data type is raw counts
#> ! [2026-09-20 21:33:46] Negative values detected
#> ! [2026-09-20 21:33:46] Negative values detected
#> ℹ [2026-09-20 21:33:46] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 21:33:46] Using 1 core
#> ⠙ [2026-09-20 21:33:46] Running for S_score [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 21:33:46] Completed 2 tasks in 82ms
#> 
#> ℹ [2026-09-20 21:33:46] Building results
#> ✔ [2026-09-20 21:33:46] Find dynamic features done
#> ℹ [2026-09-20 21:33:46] Some features were missing in at least one lineage: 
#> ℹ                       Isl1,Neurod2,Pyy,Rbp4,Sox9...
ht2$plot


ht3 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = c("Lineage1", "Lineage2"),
  n_split = 5,
  nlabel = 10,
  split_method = "kmeans",
  cluster_rows = TRUE,
  cell_annotation = "SubCellType",
  width = 1,
  height = 2
)
#> ℹ [2026-09-20 21:33:46] [1] 176 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
ht3$plot


ht4 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  cell_annotation = "SubCellType",
  n_split = 3,
  nlabel = 10,
  split_method = "mfuzz",
  species = "Mus_musculus",
  db = "GO_BP",
  anno_terms = TRUE,
  width = 1,
  height = 2
)
#> ℹ [2026-09-20 21:33:51] [1] 176 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
#> ℹ [2026-09-20 21:33:51] Start Enrichment analysis
#> ℹ [2026-09-20 21:33:51] Species: "Mus_musculus"
#> ℹ [2026-09-20 21:33:51] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-09-20 21:23:40
#> ℹ [2026-09-20 21:33:53] Permform enrichment...
#> ✔ [2026-09-20 21:33:55] Enrichment analysis done

ht5 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = "Lineage1",
  cell_annotation = "SubCellType",
  n_split = 2,
  split_method = "mfuzz",
  species = "Mus_musculus",
  db = "GO_BP",
  cores = 2,
  nlabel = 10,
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE,
  width = 1,
  height = 2,
  terms_width = grid::unit(1, "in"),
  terms_fontsize = 6,
  keys_width = grid::unit(0.5, "in"),
  keys_fontsize = c(3, 6),
  features_width = grid::unit(0.5, "in"),
  features_fontsize = c(3, 6)
)
#> ℹ [2026-09-20 21:33:58] [1] 154 features from Lineage1 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
#> ℹ [2026-09-20 21:33:59] Start Enrichment analysis
#> ℹ [2026-09-20 21:33:59] Species: "Mus_musculus"
#> ℹ [2026-09-20 21:33:59] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-09-20 21:23:40
#> ℹ [2026-09-20 21:34:00] Permform enrichment...
#> ✔ [2026-09-20 21:34:02] Enrichment analysis done

pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = c("CSPA", "TF")
)
#> ℹ [2026-09-20 21:34:20] Species: "Mus_musculus"
#> ℹ [2026-09-20 21:34:20] Loading cached: TF version: AnimalTFDB4 nterm:2 created: 2026-09-20 21:04:17
#> ℹ [2026-09-20 21:34:21] Preparing database: CSPA
ht6 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  use_fitted = TRUE,
  n_split = 3,
  nlabel = 10,
  split_method = "mfuzz",
  heatmap_palette = "viridis",
  cell_annotation = c(
    "SubCellType", "Phase", "G2M_score"
  ),
  cell_annotation_palette = c(
    "Chinese", "simspec", "Purples"
  ),
  separate_annotation = list(
    "SubCellType", c("Arxes1", "Ncoa2")
  ),
  separate_annotation_palette = c(
    "Chinese", "Set1"
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
  pseudotime_label_color = "red",
  width = 1,
  height = 2
)
#> ℹ [2026-09-20 21:34:23] [1] 176 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
#> ℹ [2026-09-20 21:34:24] Start find dynamic features
#> ℹ [2026-09-20 21:34:24] Data type is raw counts
#> ℹ [2026-09-20 21:34:25] Number of candidate features (union): 2
#> ℹ [2026-09-20 21:34:25] Data type is raw counts
#> ℹ [2026-09-20 21:34:25] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 21:34:25] Using 1 core
#> ⠙ [2026-09-20 21:34:25] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 21:34:25] Completed 2 tasks in 263ms
#> 
#> ℹ [2026-09-20 21:34:25] Building results
#> ✔ [2026-09-20 21:34:25] Find dynamic features done
#> ℹ [2026-09-20 21:34:26] Start find dynamic features
#> ℹ [2026-09-20 21:34:26] Data type is raw counts
#> ℹ [2026-09-20 21:34:26] Number of candidate features (union): 2
#> ℹ [2026-09-20 21:34:27] Data type is raw counts
#> ℹ [2026-09-20 21:34:27] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-09-20 21:34:27] Using 1 core
#> ⠙ [2026-09-20 21:34:27] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 21:34:27] Completed 2 tasks in 228ms
#> 
#> ℹ [2026-09-20 21:34:27] Building results
#> ✔ [2026-09-20 21:34:27] Find dynamic features done
#> Picking joint bandwidth of 19.9
#> Picking joint bandwidth of 17.5
#> Picking joint bandwidth of 19.9
#> Picking joint bandwidth of 17.5

ht7 <- DynamicHeatmap(
  pancreas_sub,
  exp_legend_title = "Z-score",
  lineages = c("Lineage1", "Lineage2"),
  reverse_ht = "Lineage1",
  use_fitted = TRUE,
  n_split = 3,
  nlabel = 10,
  split_method = "mfuzz",
  heatmap_palette = "viridis",
  cell_annotation = c(
    "SubCellType", "Phase", "G2M_score"
  ),
  cell_annotation_palette = c(
    "Chinese", "simspec", "Purples"
  ),
  separate_annotation = list(
    "SubCellType", c("Arxes1", "Ncoa2")
  ),
  separate_annotation_palette = c("Chinese", "Set1"),
  separate_annotation_params = list(width = grid::unit(10, "mm")),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"),
    c("forestgreen")
  ),
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  flip = TRUE,
  column_title_rot = 90,
  width = 2,
  height = 1
)
#> ℹ [2026-09-20 21:34:32] [1] 176 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
#> ℹ [2026-09-20 21:34:32] Start find dynamic features
#> ℹ [2026-09-20 21:34:33] Data type is raw counts
#> ℹ [2026-09-20 21:34:33] Number of candidate features (union): 2
#> ℹ [2026-09-20 21:34:33] Data type is raw counts
#> ℹ [2026-09-20 21:34:33] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 21:34:33] Using 1 core
#> ⠙ [2026-09-20 21:34:33] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 21:34:33] Completed 2 tasks in 225ms
#> 
#> ℹ [2026-09-20 21:34:33] Building results
#> ✔ [2026-09-20 21:34:34] Find dynamic features done
#> ℹ [2026-09-20 21:34:34] Start find dynamic features
#> ℹ [2026-09-20 21:34:34] Data type is raw counts
#> ℹ [2026-09-20 21:34:36] Number of candidate features (union): 2
#> ℹ [2026-09-20 21:34:36] Data type is raw counts
#> ℹ [2026-09-20 21:34:36] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-09-20 21:34:36] Using 1 core
#> ⠙ [2026-09-20 21:34:36] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 21:34:36] Completed 2 tasks in 144ms
#> 
#> ℹ [2026-09-20 21:34:36] Building results
#> ✔ [2026-09-20 21:34:36] Find dynamic features done
#> Picking joint bandwidth of 19.9
#> Picking joint bandwidth of 17.5
#> Picking joint bandwidth of 19.9
#> Picking joint bandwidth of 17.5
```
