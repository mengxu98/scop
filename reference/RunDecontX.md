# Run ambient RNA decontamination with decontX

Run ambient RNA decontamination with decontX

## Usage

``` r
RunDecontX(
  srt,
  assay = "RNA",
  group.by = NULL,
  batch = NULL,
  background = NULL,
  background_assay = NULL,
  bg_batch = NULL,
  assay_name = "decontXcounts",
  store_assay = TRUE,
  round_counts = FALSE,
  seed = 11,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for decontamination. Default is
  `"RNA"`.

- group.by:

  Cell cluster labels passed to `decontX::decontX()`. Can be `NULL`, a
  meta.data column name, or a vector aligned to cells. Default is
  `NULL`.

- batch:

  Batch labels passed to `decontX::decontX()`. Can be `NULL`, a
  meta.data column name, or a vector aligned to cells. Default is
  `NULL`.

- background:

  Optional background / empty-droplet input passed to
  `decontX::decontX()`. Can be a `Seurat` object,
  `SingleCellExperiment`, or count matrix. Default is `NULL`.

- background_assay:

  Assay name used when `background` is a `Seurat` object or
  `SingleCellExperiment`. Default is `NULL`, which falls back to `assay`
  for `Seurat` background and `"counts"` for `SingleCellExperiment`
  background.

- bg_batch:

  Batch labels for `background` passed to `decontX::decontX()`. Can be
  `NULL`, a metadata column name, or a vector aligned to the background
  droplets. Default is `NULL`.

- assay_name:

  Name of the assay used to store decontaminated counts. Default is
  `"decontXcounts"`.

- store_assay:

  Whether to store decontaminated counts as a new assay. Default is
  `TRUE`.

- round_counts:

  Whether to round decontaminated counts before creating the assay.
  Default is `FALSE`.

- seed:

  Random seed for reproducibility. Default is `11`.

- ...:

  Additional arguments passed to `decontX::decontX()`.

## Value

Returns a Seurat object with decontX contamination estimates stored in
the meta.data, and optional decontaminated counts stored in a new assay.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 16:38:43] Start standard processing workflow...
#> ℹ [2026-04-02 16:38:44] Checking a list of <Seurat>...
#> ! [2026-04-02 16:38:44] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:38:44] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:38:45] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:38:46] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:38:46] Number of available HVF: 2000
#> ℹ [2026-04-02 16:38:46] Finished check
#> ℹ [2026-04-02 16:38:46] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:38:47] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:38:50] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:38:51] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:38:51] Reorder clusters...
#> ℹ [2026-04-02 16:38:51] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:38:51] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:38:51] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:38:51] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:38:51] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:38:54] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:38:57] Standard processing workflow completed
pancreas_sub <- RunDecontX(
  pancreas_sub,
  group.by = "CellType"
)
#> ℹ [2026-04-02 16:38:58] Data type is raw counts
#> ℹ [2026-04-02 16:39:05] Running decontX
#> Error in loadNamespace(x): there is no package called ‘decontX’

FeatureStatPlot(
  pancreas_sub,
  stat.by = "decontX_contamination"
)
#> ! [2026-04-02 16:39:05] "decontX_contamination" are not found
#> Error in ExpressionStatPlot(exp.data = exp.data, meta.data = meta.data,     stat.by = stat.by, group.by = group.by, split.by = split.by,     bg.by = bg.by, plot.by = "group", fill.by = fill.by, cells = cells,     keep_empty = keep_empty, individual = individual, plot_type = plot_type,     palette = palette, palcolor = palcolor, alpha = alpha, bg_palette = bg_palette,     bg_palcolor = bg_palcolor, bg_alpha = bg_alpha, add_box = add_box,     box_color = box_color, box_width = box_width, box_ptsize = box_ptsize,     add_point = add_point, pt.color = pt.color, pt.size = pt.size,     pt.alpha = pt.alpha, jitter.width = jitter.width, jitter.height = jitter.height,     add_trend = add_trend, trend_color = trend_color, trend_linewidth = trend_linewidth,     trend_ptsize = trend_ptsize, add_stat = add_stat, stat_color = stat_color,     stat_size = stat_size, stat_stroke = stat_stroke, stat_shape = stat_shape,     add_line = add_line, line_color = line_color, line_size = line_size,     line_type = line_type, cells.highlight = cells.highlight,     cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,     alpha.highlight = alpha.highlight, calculate_coexp = calculate_coexp,     same.y.lims = same.y.lims, y.min = y.min, y.max = y.max,     y.trans = y.trans, y.nbreaks = y.nbreaks, sort = sort, stack = stack,     flip = flip, comparisons = comparisons, ref_group = ref_group,     pairwise_method = pairwise_method, multiplegroup_comparisons = multiplegroup_comparisons,     multiple_method = multiple_method, sig_label = sig_label,     sig_labelsize = sig_labelsize, aspect.ratio = aspect.ratio,     title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,     legend.position = legend.position, legend.direction = legend.direction,     legend.title = legend.title, theme_use = theme_use, theme_args = theme_args,     force = force, seed = seed): `stat.by` must be type of numeric variable

FeatureDimPlot(
  pancreas_sub,
  features = "decontX_contamination"
)
#> ! [2026-04-02 16:39:05] "decontX_contamination" are not in the features of <Seurat>
#> Error in FeatureDimPlot(pancreas_sub, features = "decontX_contamination"): There are no valid features present.
```
