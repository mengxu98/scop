# Differential gene test

This function utilizes the Seurat package to perform a differential
expression (DE) test on gene expression data. Users have the flexibility
to specify custom cell groups, marker types, and various options for DE
analysis.

## Usage

``` r
RunDEtest(
  srt,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  features = NULL,
  feature_type = c("gene", "peak", "cCRE"),
  analysis_level = c("cell", "pseudobulk"),
  markers_type = c("all", "paired", "conserved", "disturbed"),
  grouping.var = NULL,
  meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
  test.use = "wilcox",
  only.pos = TRUE,
  fc.threshold = 1.5,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  norm.method = "LogNormalize",
  sample_col = NULL,
  condition_col = NULL,
  p.adjust.method = "bonferroni",
  layer = "data",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  A grouping variable in the dataset to define the groups or conditions
  for the differential test. If not provided, the function uses the
  "active.ident" variable in the Seurat object.

- group1:

  A vector of cell IDs or a character vector specifying the cells that
  belong to the first group. If both group.by and group1 are provided,
  group1 takes precedence. For pseudobulk analysis, this parameter is
  interpreted as the first condition label.

- group2:

  A vector of cell IDs or a character vector specifying the cells that
  belong to the second group. This parameter is only used when group.by
  or group1 is provided. For pseudobulk analysis, this parameter is
  interpreted as the second condition label.

- cells1:

  A vector of cell IDs specifying the cells that belong to group1. If
  provided, group1 is ignored.

- cells2:

  A vector of cell IDs specifying the cells that belong to group2. This
  parameter is only used when cells1 is provided.

- features:

  A vector of feature names specifying the features to consider for the
  differential test. If not provided, all features in the dataset are
  considered.

- feature_type:

  Feature type used for differential testing. Default is `"gene"`.

- analysis_level:

  Analysis level used for differential testing. Default is `"cell"`.

- markers_type:

  A character value specifying the type of markers to find. Possible
  values are "all", "paired", "conserved", and "disturbed". Pseudobulk
  analysis currently supports only `"all"`.

- grouping.var:

  A character value specifying the grouping variable for finding
  conserved or disturbed markers. This parameter is only used when
  markers_type is "conserved" or "disturbed".

- meta.method:

  A character value specifying the method to use for combining p-values
  in the conserved markers test. Possible values are "maximump",
  "minimump", "wilkinsonp", "meanp", "sump", and "votep".

- test.use:

  Differential testing method. For pseudobulk analysis, only
  `"limma_voom"` and `"edgeR"` are currently supported.

- only.pos:

  Only return positive markers (FALSE by default)

- fc.threshold:

  A numeric value used to filter genes for testing based on their
  average fold change between/among the two groups. Default is `1.5`.

- base:

  The base with respect to which logarithms are computed.

- pseudocount.use:

  Pseudocount to add to averaged expression values when calculating
  logFC. 1 by default.

- mean.fxn:

  Function to use for fold change or average difference calculation. The
  default depends on the the value of `fc.slot`:

  - "counts" : difference in the log of the mean counts, with
    pseudocount.

  - "data" : difference in the log of the average exponentiated data,
    with pseudocount. This adjusts for differences in sequencing depth
    between cells, and assumes that "data" has been log-normalized.

  - "scale.data" : difference in the means of scale.data.

- min.pct:

  only test genes that are detected in a minimum fraction of min.pct
  cells in either of the two populations. Meant to speed up the function
  by not testing genes that are very infrequently expressed. Default is
  0.01

- min.diff.pct:

  only test genes that show a minimum difference in the fraction of
  detection between the two groups. Set to -Inf by default

- max.cells.per.ident:

  Down sample each identity class to a max number. Default is no
  downsampling. Not activated by default (set to Inf)

- latent.vars:

  Variables to test, used only when `test.use` is one of 'LR',
  'negbinom', 'poisson', or 'MAST'

- min.cells.feature:

  Minimum number of cells expressing the feature in at least one of the
  two groups, currently only used for poisson and negative binomial
  tests

- min.cells.group:

  Minimum number of cells in one of the groups

- norm.method:

  Normalization method for fold change calculation when layer is 'data'.
  Default is `"LogNormalize"`.

- sample_col:

  Metadata column storing biological sample IDs for pseudobulk analysis.
  Required when `analysis_level = "pseudobulk"`.

- condition_col:

  Metadata column storing condition labels for pseudobulk analysis.
  Required when `analysis_level = "pseudobulk"`.

- p.adjust.method:

  A character value specifying the method to use for adjusting p-values.
  Default is `"bonferroni"`.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Assay to use in differential expression testing

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- ...:

  Additional arguments to pass to the
  [Seurat::FindMarkers](https://satijalab.org/seurat/reference/FindMarkers.html)
  function.

## See also

[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-27 16:30:56] Start standard processing workflow...
#> ℹ [2026-04-27 16:30:57] Checking a list of <Seurat>...
#> ! [2026-04-27 16:30:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:30:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-27 16:30:59] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-27 16:31:00] Use the separate HVF from `srt_list`
#> ℹ [2026-04-27 16:31:00] Number of available HVF: 2000
#> ℹ [2026-04-27 16:31:00] Finished check
#> ℹ [2026-04-27 16:31:00] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-27 16:31:01] Perform pca linear dimension reduction
#> ℹ [2026-04-27 16:31:01] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-27 16:31:02] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-27 16:31:02] Reorder clusters...
#> ℹ [2026-04-27 16:31:02] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-27 16:31:02] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-27 16:31:02] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-27 16:31:06] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-27 16:31:11] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType"
)
#> ℹ [2026-04-27 16:31:11] Data type is log-normalized
#> ℹ [2026-04-27 16:31:11] Start differential expression test
#> ℹ [2026-04-27 16:31:11] Find all markers(wilcox) among [1] 8 groups...
#> ℹ [2026-04-27 16:31:11] Using 1 core
#> ⠙ [2026-04-27 16:31:11] Running for Ductal [1/8] ■           12% | ETA:  1s
#> ⠹ [2026-04-27 16:31:11] Running for Pre-endocrine [5/8] ■■■■■■      62% | ETA: …
#> ✔ [2026-04-27 16:31:11] Completed 8 tasks in 1.3s
#> 
#> ℹ [2026-04-27 16:31:11] Building results
#> ✔ [2026-04-27 16:31:12] Differential expression test completed
AllMarkers <- dplyr::filter(
  pancreas_sub@tools$DEtest_SubCellType$AllMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
ht1 <- GroupHeatmap(
  pancreas_sub,
  features = AllMarkers$gene,
  feature_split = AllMarkers$group1,
  group.by = "SubCellType"
)
#> ℹ [2026-04-27 16:31:15] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:31:15] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:31:15] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht1$plot


TopMarkers <- AllMarkers |>
  dplyr::group_by(gene) |>
  dplyr::top_n(1, avg_log2FC) |>
  dplyr::group_by(group1) |>
  dplyr::top_n(3, avg_log2FC)
ht2 <- GroupHeatmap(
  pancreas_sub,
  features = TopMarkers$gene,
  feature_split = TopMarkers$group1,
  group.by = "SubCellType",
  show_row_names = TRUE
)
#> ℹ [2026-04-27 16:31:18] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:31:18] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:31:18] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht2$plot


pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType",
  markers_type = "paired",
  cores = 2
)
#> ℹ [2026-04-27 16:31:19] Data type is log-normalized
#> ℹ [2026-04-27 16:31:19] Start differential expression test
#> ℹ [2026-04-27 16:31:19] Find paired markers(wilcox) among [1] 8 groups...
#> ℹ [2026-04-27 16:31:19] Using 2 cores
#> ⠙ [2026-04-27 16:31:19] Running for 1 [1/56]              2% | ETA: 14s
#> ⠹ [2026-04-27 16:31:19] Running for 13 [13/56] ■■          23% | ETA:  6s
#> ⠸ [2026-04-27 16:31:19] Running for 38 [38/56] ■■■■■■      68% | ETA:  2s
#> ✔ [2026-04-27 16:31:19] Completed 56 tasks in 7s
#> 
#> ℹ [2026-04-27 16:31:19] Building results
#> ✔ [2026-04-27 16:31:26] Differential expression test completed
PairedMarkers <- dplyr::filter(
  pancreas_sub@tools$DEtest_SubCellType$PairedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
ht3 <- GroupHeatmap(
  pancreas_sub,
  features = PairedMarkers$gene,
  feature_split = PairedMarkers$group1,
  group.by = "SubCellType"
)
#> ℹ [2026-04-27 16:32:17] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:32:17] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:32:17] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht3$plot


data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-04-27 16:32:20] Run integration workflow...
#> ℹ [2026-04-27 16:32:21] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-27 16:32:21] Checking a list of <Seurat>...
#> ! [2026-04-27 16:32:22] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:32:22] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-04-27 16:32:24] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-04-27 16:32:24] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:32:24] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-04-27 16:32:26] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-04-27 16:32:26] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:32:26] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-04-27 16:32:28] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-04-27 16:32:28] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:32:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-04-27 16:32:30] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-04-27 16:32:31] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:32:31] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-04-27 16:32:33] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-04-27 16:32:33] Use the separate HVF from `srt_list`
#> ℹ [2026-04-27 16:32:34] Number of available HVF: 2000
#> ℹ [2026-04-27 16:32:34] Finished check
#> ℹ [2026-04-27 16:32:37] Perform Uncorrected integration
#> ℹ [2026-04-27 16:32:38] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-27 16:32:38] Perform "pca" linear dimension reduction
#> ℹ [2026-04-27 16:32:41] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-04-27 16:32:42] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-27 16:32:42] Reorder clusters...
#> ℹ [2026-04-27 16:32:42] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-27 16:32:42] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-04-27 16:32:47] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-04-27 16:32:52] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ✔ [2026-04-27 16:32:59] Uncorrected integration completed
CellDimPlot(
  panc8_sub,
  group.by = c("celltype", "tech")
)


panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-04-27 16:33:00] Data type is log-normalized
#> ℹ [2026-04-27 16:33:00] Start differential expression test
#> ℹ [2026-04-27 16:33:00] Find conserved markers(wilcox) among [1] 13 groups...
#> ℹ [2026-04-27 16:33:00] Using 2 cores
#> ⠙ [2026-04-27 16:33:00] Running for delta [1/13]              8% | ETA: 17s
#> ⠹ [2026-04-27 16:33:00] Running for acinar [3/13] ■■          23% | ETA: 10s
#> ⠸ [2026-04-27 16:33:00] Running for ductal [6/13] ■■■■        46% | ETA:  6s
#> ✔ [2026-04-27 16:33:00] Completed 13 tasks in 8.8s
#> 
#> ℹ [2026-04-27 16:33:00] Building results
#> ✔ [2026-04-27 16:33:09] Differential expression test completed
ConservedMarkers1 <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
ht4 <- GroupHeatmap(
  panc8_sub,
  layer = "data",
  features = ConservedMarkers1$gene,
  feature_split = ConservedMarkers1$group1,
  group.by = "tech",
  split.by = "celltype",
  within_groups = TRUE
)
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> ℹ [2026-04-27 16:33:20] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:33:20] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:33:20] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht4$plot


panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "tech",
  grouping.var = "celltype",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-04-27 16:33:30] Data type is log-normalized
#> ℹ [2026-04-27 16:33:30] Start differential expression test
#> ℹ [2026-04-27 16:33:30] Find conserved markers(wilcox) among [1] 5 groups...
#> ℹ [2026-04-27 16:33:30] Using 2 cores
#> ⠙ [2026-04-27 16:33:30] Running for celseq [1/5] ■■          20% | ETA: 10s
#> ⠹ [2026-04-27 16:33:30] Running for smartseq2 [3/5] ■■■■■■      60% | ETA:  4s
#> ✔ [2026-04-27 16:33:30] Completed 5 tasks in 8.7s
#> 
#> ℹ [2026-04-27 16:33:30] Building results
#> ✔ [2026-04-27 16:33:39] Differential expression test completed
ConservedMarkers2 <- dplyr::filter(
  panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
ht4 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = ConservedMarkers2$gene,
  feature_split = ConservedMarkers2$group1,
  group.by = "tech",
  split.by = "celltype"
)
#> ℹ [2026-04-27 16:33:42] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:33:42] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:33:42] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht4$plot


panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "disturbed",
  cores = 2
)
#> ℹ [2026-04-27 16:33:48] Data type is log-normalized
#> ℹ [2026-04-27 16:33:48] Start differential expression test
#> ℹ [2026-04-27 16:33:48] Find disturbed markers(wilcox) among [1] 13 groups...
#> ℹ [2026-04-27 16:33:48] Using 2 cores
#> ⠙ [2026-04-27 16:33:48] Running for gamma [1/13]              8% | ETA: 39s
#> ⠹ [2026-04-27 16:33:48] Running for acinar [3/13] ■■          23% | ETA: 22s
#> ⠸ [2026-04-27 16:33:48] Running for alpha [5/13] ■■■         38% | ETA: 17s
#> ⠼ [2026-04-27 16:33:48] Running for mast [7/13] ■■■■■       54% | ETA: 11s
#> ⠴ [2026-04-27 16:33:48] Running for quiescent-stellate [10/13] ■■■■■■■     77% …
#> ✔ [2026-04-27 16:33:48] Completed 13 tasks in 17.3s
#> 
#> ℹ [2026-04-27 16:33:48] Building results
#> ✔ [2026-04-27 16:34:05] Differential expression test completed
DisturbedMarkers <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
)
ht5 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = DisturbedMarkers$gene,
  feature_split = DisturbedMarkers$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> `use_raster` is automatically set to TRUE for a matrix with more than
#> 2000 rows. You can control `use_raster` argument by explicitly setting
#> TRUE/FALSE to it.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
#> ℹ [2026-04-27 16:34:17] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:34:17] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:34:17] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht5$plot


gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
DisturbedMarkers_specific <- DisturbedMarkers[
  DisturbedMarkers$gene %in% gene_specific,
]
ht6 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> ℹ [2026-04-27 16:34:35] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:34:35] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:34:35] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht6$plot


ht7 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  aggregate_fun = function(x) mean(expm1(x)) + 1,
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  grouping.var = "tech",
  numerator = "smartseq2"
)
#> ! [2026-04-27 16:34:47] When 'grouping.var' is specified, 'exp_method' can only be 'log2fc'
#> ℹ [2026-04-27 16:34:49] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-04-27 16:34:49] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-04-27 16:34:49] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht7$plot


pbmc_small <- UpdateSeuratObject(pbmc_small)
#> Error in UpdateSeuratObject(pbmc_small): could not find function "UpdateSeuratObject"
pbmc_small[["sample"]] <- rep(c("S1", "S2", "S3", "S4"), length.out = ncol(pbmc_small))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'pbmc_small' not found
pbmc_small[["condition"]] <- rep(c("ctrl", "ctrl", "case", "case"), length.out = ncol(pbmc_small))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'pbmc_small' not found
pbmc_small <- RunDEtest(
  pbmc_small,
  analysis_level = "pseudobulk",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "limma_voom",
  layer = "counts",
  cores = 1
)
#> Error: object 'pbmc_small' not found
pbmc_small <- RunDEtest(
  pbmc_small,
  analysis_level = "pseudobulk",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "edgeR",
  layer = "counts",
  cores = 1
)
#> Error: object 'pbmc_small' not found
edgeR_markers <- pbmc_small@tools$DEtest_pseudobulk$AllMarkers_edgeR
#> Error: object 'pbmc_small' not found

# \donttest{
data(pbmcmultiome_sub)
pbmcmultiome_sub[["sample"]] <- rep(
  c("S1", "S2", "S3", "S4"),
  length.out = ncol(pbmcmultiome_sub)
)
pbmcmultiome_sub[["condition"]] <- rep(
  c("ctrl", "ctrl", "case", "case"),
  length.out = ncol(pbmcmultiome_sub)
)
pbmcmultiome_sub <- RunDEtest(
  pbmcmultiome_sub,
  assay = "peaks",
  layer = "counts",
  feature_type = "peak",
  analysis_level = "pseudobulk",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "edgeR",
  cores = 1
)
#> ℹ [2026-04-27 16:34:52] Start pseudobulk differential testing
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-04-27 16:34:53] Pseudobulk differential testing completed
peak_markers <- pbmcmultiome_sub@tools$DEtest_pseudobulk$AllMarkers_edgeR
# }
```
