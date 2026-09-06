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
  row_names_wrap = NULL,
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
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  cell_split_palette = "simspec",
  cell_split_palcolor = NULL,
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
  use_raster = NULL,
  raster_device = "png",
  raster_by_magick = FALSE,
  height = NULL,
  width = NULL,
  units = "inch",
  cores = 1,
  seed = 11,
  legend.position = "right",
  ht_params = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- features:

  Features to plot.

- cells:

  Cell names to include.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column to facet by.

- within_groups:

  Separate color scales per group.

- max_cells:

  An integer, maximum number of cells to sample per group.

- cell_order:

  A vector of cell names defining the order of cells.

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

- feature_split, feature_split_by, n_split, split_order, split_method,
  decreasing:

  Feature splitting. `split_method` is `"kmeans"`, `"hclust"`, or
  `"mfuzz"`.

- fuzzification:

  Mfuzz fuzzification coefficient.

- cluster_features_by:

  Grouping used when clustering features. `NULL` uses all groups.

- cluster_rows, cluster_columns, cluster_row_slices,
  cluster_column_slices:

  Heatmap clustering.

- show_row_names, show_column_names, row_names_wrap:

  Show names. `row_names_wrap` wraps displayed names (underscores as
  spaces) without changing feature IDs.

- row_names_side, column_names_side, row_names_rot, column_names_rot:

  Name placement.

- row_title, column_title, row_title_side, column_title_side,
  row_title_rot, column_title_rot:

  Slice titles.

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

- heatmap_palette, heatmap_palcolor, group_palette, group_palcolor:

  Heatmap and group colors.

- cell_split_palette, cell_split_palcolor, feature_split_palette,
  feature_split_palcolor:

  Split colors.

- cell_annotation, cell_annotation_palette, cell_annotation_palcolor,
  cell_annotation_params:

  Cell annotations. Palette length should match `cell_annotation`.

- feature_annotation, feature_annotation_palette,
  feature_annotation_palcolor, feature_annotation_params:

  Feature annotations. Palette length should match `feature_annotation`.

- use_raster, raster_device, raster_by_magick:

  Raster device (`NULL` chooses automatically).

- width, height, units:

  Heatmap size. `NULL` sizes from matrix dimensions.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

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

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments passed to helper functions.

## See also

[RunDEtest](https://mengxu98.github.io/scop/reference/RunDEtest.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:29:31] Start standard processing workflow...
#> ℹ [2026-09-06 21:29:32] Checking a list of <Seurat>...
#> ! [2026-09-06 21:29:32] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:29:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:29:32] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:29:32] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:29:32] Number of available HVF: 2000
#> ℹ [2026-09-06 21:29:32] Finished check
#> ℹ [2026-09-06 21:29:32] Perform `ScaleData()`
#> ℹ [2026-09-06 21:29:32] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:29:32] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:29:33] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:29:33] Reorder clusters...
#> ℹ [2026-09-06 21:29:33] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:29:33] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:29:39] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-09-06 21:29:39] Data type is log-normalized
#> ℹ [2026-09-06 21:29:39] Start differential expression test
#> ℹ [2026-09-06 21:29:39] Find all markers(wilcox) among [1] 5 groups...
#> ℹ [2026-09-06 21:29:39] Using 1 core
#> ⠙ [2026-09-06 21:29:39] Running for Ductal [1/5] ■■          20% | ETA: 15s
#> ⠹ [2026-09-06 21:29:39] Running for Ngn3-high-EP [2/5] ■■■■        40% | ETA:  …
#> ⠸ [2026-09-06 21:29:39] Running for Endocrine [3/5] ■■■■■■      60% | ETA:  6s
#> ⠼ [2026-09-06 21:29:39] Running for Ngn3-low-EP [4/5] ■■■■■■■■    80% | ETA:  3s
#> ✔ [2026-09-06 21:29:39] Completed 5 tasks in 14.8s
#> 
#> ℹ [2026-09-06 21:29:39] Building results
#> ✔ [2026-09-06 21:29:54] Differential expression test completed
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
#> ℹ [2026-09-06 21:30:06] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:30:06] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:30:06] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
#> ℹ [2026-09-06 21:30:34] Start Enrichment analysis
#> ℹ [2026-09-06 21:30:34] Species: "Mus_musculus"
#> ℹ [2026-09-06 21:30:34] Loading cached: GO_BP version: 3.23.0 nterm:14957 created: 2026-09-06 21:27:14
#> ℹ [2026-09-06 21:30:36] Permform enrichment...
#> ℹ [2026-09-06 21:30:37] Using 1 core
#> ⠙ [2026-09-06 21:30:37] Running for 1 [1/5] ■■          20% | ETA:  3s
#> ✔ [2026-09-06 21:30:37] Completed 5 tasks in 2.9s
#> 
#> ℹ [2026-09-06 21:30:37] Building results
#> ✔ [2026-09-06 21:30:41] Enrichment analysis done
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> ℹ [2026-09-06 21:31:24] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:31:24] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:31:24] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
  cell_annotation_palette = c("Chinese", "cividis")
)
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
ht4$plot


pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = c("CSPA", "TF")
)
#> ℹ [2026-09-06 21:31:38] Species: "Mus_musculus"
#> ℹ [2026-09-06 21:31:38] Loading cached: CSPA version: CSPA nterm:1 created: 2026-09-06 21:25:21
#> ℹ [2026-09-06 21:31:39] Loading cached: TF version: AnimalTFDB4 nterm:2 created: 2026-09-06 20:52:04

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
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> ℹ [2026-09-06 21:31:43] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:31:43] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:31:43] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 columns You can control `use_raster` argument by explicitly
#> setting TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> ℹ [2026-09-06 21:31:54] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:31:54] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:31:55] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht6$plot
```
