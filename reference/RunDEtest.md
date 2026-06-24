# Differential gene test

Perform differential expression testing on a `Seurat` object or a
`SummarizedExperiment` bulk object. Users have the flexibility to
specify custom cell groups, marker types, and various options for DE
analysis.

## Usage

``` r
RunDEtest(object = NULL, ..., srt = NULL)

# S3 method for class 'Seurat'
RunDEtest(
  object,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  ident.1 = NULL,
  ident.2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  feature_type = c("gene", "peak", "cCRE"),
  markers_type = c("all", "paired", "conserved", "disturbed"),
  grouping.var = NULL,
  meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
  test.use = "wilcox",
  only.pos = TRUE,
  fc.threshold = 1.5,
  logfc.threshold = NULL,
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
  bulk_assay = "counts",
  p.adjust.method = "bonferroni",
  layer = "data",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
)

# S3 method for class 'SummarizedExperiment'
RunDEtest(
  object,
  group.by = NULL,
  group1 = NULL,
  group2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  features = NULL,
  feature_type = c("gene", "peak", "cCRE"),
  markers_type = c("all", "paired", "conserved", "disturbed"),
  grouping.var = NULL,
  meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
  test.use = "edgeR",
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
  bulk_assay = "counts",
  p.adjust.method = "bonferroni",
  layer = "counts",
  assay = NULL,
  seed = 11,
  verbose = TRUE,
  cores = 1,
  ...
)
```

## Arguments

- object:

  A `Seurat` object or a `SummarizedExperiment` object.

- ...:

  Additional arguments to pass to the
  [Seurat::FindMarkers](https://satijalab.org/seurat/reference/FindMarkers.html)
  function.

- srt:

  Compatibility alias for `object`.

- group.by:

  A grouping variable in the dataset to define the groups or conditions
  for the differential test. If not provided, the function uses the
  "active.ident" variable in the Seurat object.

- group1:

  A vector of cell IDs or a character vector specifying the cells that
  belong to the first group. If both group.by and group1 are provided,
  group1 takes precedence. For sample-level methods (`"edgeR"`,
  `"limma"`, `"DESeq2"`, and `"dream"`), this parameter is interpreted
  as the first condition label.

- group2:

  A vector of cell IDs or a character vector specifying the cells that
  belong to the second group. This parameter is only used when group.by
  or group1 is provided. For sample-level methods (`"edgeR"`, `"limma"`,
  `"DESeq2"`, and `"dream"`), this parameter is interpreted as the
  second condition label.

- ident.1, ident.2:

  Seurat-style aliases for `group1` and `group2`.

- cells1:

  A vector of cell IDs specifying the cells that belong to group1. If
  provided, group1 is ignored.

- cells2:

  A vector of cell IDs specifying the cells that belong to group2. This
  parameter is only used when cells1 is provided.

- cells.1, cells.2:

  Seurat-style aliases for `cells1` and `cells2`.

- features:

  A vector of feature names specifying the features to consider for the
  differential test. If not provided, all features in the dataset are
  considered.

- feature_type:

  Feature type used for differential testing. Default is `"gene"`.

- markers_type:

  A character value specifying the type of markers to find. Possible
  values are "all", "paired", "conserved", and "disturbed". Sample-level
  methods (`"edgeR"`, `"limma"`, `"DESeq2"`, and `"dream"`) currently
  support only `"all"`.

- grouping.var:

  A character value specifying the grouping variable for finding
  conserved or disturbed markers. This parameter is only used when
  markers_type is "conserved" or "disturbed".

- meta.method:

  A character value specifying the method to use for combining p-values
  in the conserved markers test. Possible values are "maximump",
  "minimump", "wilkinsonp", "meanp", "sump", and "votep".

- test.use:

  Differential testing method. `"edgeR"`, `"limma"`, `"DESeq2"`, and
  `"dream"` run sample-level pseudobulk differential testing on `Seurat`
  input and bulk DE on `SummarizedExperiment` input.

- only.pos:

  Only return positive markers (FALSE by default)

- fc.threshold:

  A numeric value used to filter genes for testing based on their
  average fold change between/among the two groups. Default is `1.5`.

- logfc.threshold:

  Seurat-style log fold-change threshold. When provided, it is converted
  to `fc.threshold = base^logfc.threshold`.

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

  Metadata column storing biological sample IDs. Required when
  `test.use` is a sample-level pseudobulk method on `Seurat`.

- condition_col:

  Metadata column storing condition labels. Required when `test.use` is
  a sample-level pseudobulk method on `Seurat`, and required for
  `SummarizedExperiment` input.

- bulk_assay:

  Assay name used as the bulk counts matrix for `SummarizedExperiment`
  input.

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

## See also

[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-24 04:05:05] Start standard processing workflow...
#> ℹ [2026-06-24 04:05:06] Checking a list of <Seurat>...
#> ! [2026-06-24 04:05:06] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:05:06] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:05:06] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:05:06] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 04:05:06] Number of available HVF: 2000
#> ℹ [2026-06-24 04:05:06] Finished check
#> ℹ [2026-06-24 04:05:06] Perform `ScaleData()`
#> ℹ [2026-06-24 04:05:06] Perform pca linear dimension reduction
#> ℹ [2026-06-24 04:05:07] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-24 04:05:07] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-24 04:05:07] Reorder clusters...
#> ℹ [2026-06-24 04:05:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-24 04:05:08] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-24 04:05:15] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType",
  only.pos = FALSE
)
#> ℹ [2026-06-24 04:05:15] Data type is log-normalized
#> ℹ [2026-06-24 04:05:15] Start differential expression test
#> ℹ [2026-06-24 04:05:15] Find all markers(wilcox) among [1] 8 groups...
#> ℹ [2026-06-24 04:05:15] Using 1 core
#> ⠙ [2026-06-24 04:05:15] Running for Ductal [1/8] ■           12% | ETA:  2s
#> ✔ [2026-06-24 04:05:15] Completed 8 tasks in 2.3s
#> 
#> ℹ [2026-06-24 04:05:15] Building results
#> ✔ [2026-06-24 04:05:17] Differential expression test completed
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
#> ℹ [2026-06-24 04:05:20] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-06-24 04:05:20] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-06-24 04:05:20] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

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
#> ℹ [2026-06-24 04:05:23] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-06-24 04:05:23] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-06-24 04:05:23] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht2$plot


pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType",
  markers_type = "paired",
  cores = 2
)
#> ℹ [2026-06-24 04:05:24] Data type is log-normalized
#> ℹ [2026-06-24 04:05:24] Start differential expression test
#> ℹ [2026-06-24 04:05:24] Find paired markers(wilcox) among [1] 8 groups...
#> ℹ [2026-06-24 04:05:24] Using 2 cores
#> ⠙ [2026-06-24 04:05:24] Running for 1 [1/56]              2% | ETA: 25s
#> ⠹ [2026-06-24 04:05:24] Running for 10 [9/56] ■           16% | ETA: 14s
#> ⠸ [2026-06-24 04:05:24] Running for 21 [21/56] ■■■         38% | ETA:  9s
#> ⠼ [2026-06-24 04:05:24] Running for 33 [33/56] ■■■■■       59% | ETA:  6s
#> ⠴ [2026-06-24 04:05:24] Running for 46 [45/56] ■■■■■■■■    80% | ETA:  3s
#> ✔ [2026-06-24 04:05:24] Completed 56 tasks in 14.1s
#> 
#> ℹ [2026-06-24 04:05:24] Building results
#> ✔ [2026-06-24 04:05:38] Differential expression test completed
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
#> ℹ [2026-06-24 04:06:24] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-06-24 04:06:24] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-06-24 04:06:24] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht3$plot


data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-06-24 04:06:27] Run integration workflow...
#> ℹ [2026-06-24 04:06:28] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-06-24 04:06:29] Checking a list of <Seurat>...
#> ! [2026-06-24 04:06:29] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:06:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-06-24 04:06:29] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-06-24 04:06:29] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:06:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-06-24 04:06:29] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-06-24 04:06:29] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:06:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-06-24 04:06:29] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-06-24 04:06:29] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:06:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-06-24 04:06:29] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-06-24 04:06:30] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:06:30] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-06-24 04:06:30] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-06-24 04:06:30] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 04:06:30] Number of available HVF: 2000
#> ℹ [2026-06-24 04:06:31] Finished check
#> ℹ [2026-06-24 04:06:33] Perform Uncorrected integration
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-06-24 04:06:33] Perform `Seurat::ScaleData()`
#> Error: ScaleData.Seurat requires an Assay5 object with a data layer.
CellDimPlot(
  panc8_sub,
  group.by = c("celltype", "tech")
)
#> Error in DefaultReduction(srt): Unable to find any reductions

panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "conserved",
  cores = 2
)
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-06-24 04:07:06] Infinite values detected
#> ! [2026-06-24 04:07:06] Data in the 'data' layer is unknown. Please check the data type
#> ℹ [2026-06-24 04:07:06] Start differential expression test
#> ℹ [2026-06-24 04:07:06] Find conserved markers(wilcox) among [1] 13 groups...
#> ℹ [2026-06-24 04:07:06] Using 2 cores
#> ⠙ [2026-06-24 04:07:06] Running for delta [1/13]              8% | ETA:  3s
#> ⠹ [2026-06-24 04:07:06] Running for quiescent-stellate [10/13] ■■■■■■■     77% …
#> ✔ [2026-06-24 04:07:06] Completed 13 tasks in 2.2s
#> 
#> ℹ [2026-06-24 04:07:06] Building results
#> ! [2026-06-24 04:07:06] Found 11 failed results
#> ℹ [2026-06-24 04:07:09] ✖ Error details:
#> ℹ                       ✖ error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds (11): "delta", "gamma", "acinar" and 8 more
#> Error in x[, c("avg_log2FC", "pct.1", "pct.2", "max_pval", "p_val", "p_val_adj",     "gene", "group1", "group2")]: incorrect number of dimensions
ConservedMarkers1 <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht4 <- GroupHeatmap(
  panc8_sub,
  layer = "data",
  features = ConservedMarkers1$gene,
  feature_split = ConservedMarkers1$group1,
  group.by = "tech",
  split.by = "celltype",
  within_groups = TRUE
)
#> Error: object 'ConservedMarkers1' not found
ht4$plot
#> Error: object 'ht4' not found

panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "tech",
  grouping.var = "celltype",
  markers_type = "conserved",
  cores = 2
)
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-06-24 04:07:09] Infinite values detected
#> ! [2026-06-24 04:07:09] Data in the 'data' layer is unknown. Please check the data type
#> ℹ [2026-06-24 04:07:09] Start differential expression test
#> ℹ [2026-06-24 04:07:09] Find conserved markers(wilcox) among [1] 5 groups...
#> ℹ [2026-06-24 04:07:09] Using 2 cores
#> ℹ [2026-06-24 04:07:09] Building results
#> ! [2026-06-24 04:07:09] Found 5 failed results
#> ℹ [2026-06-24 04:07:10] ✖ Error details:
#> ℹ                       ✖ error in evaluating the argument 'x' in selecting a method for function 'rowSums': subscript out of bounds (5): "celseq", "celseq2", "smartseq2" and 2 more
#> Error in x[, c("avg_log2FC", "pct.1", "pct.2", "max_pval", "p_val", "p_val_adj",     "gene", "group1", "group2")]: incorrect number of dimensions
ConservedMarkers2 <- dplyr::filter(
  panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht4 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = ConservedMarkers2$gene,
  feature_split = ConservedMarkers2$group1,
  group.by = "tech",
  split.by = "celltype"
)
#> Error: object 'ConservedMarkers2' not found
ht4$plot
#> Error: object 'ht4' not found

panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "disturbed",
  cores = 2
)
#> Warning: Layer ‘data’ is empty
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> ! [2026-06-24 04:07:11] Infinite values detected
#> ! [2026-06-24 04:07:11] Data in the 'data' layer is unknown. Please check the data type
#> ℹ [2026-06-24 04:07:11] Start differential expression test
#> ℹ [2026-06-24 04:07:11] Find disturbed markers(wilcox) among [1] 13 groups...
#> ℹ [2026-06-24 04:07:11] Using 2 cores
#> ⠙ [2026-06-24 04:07:11] Running for delta [1/13]              8% | ETA: 15s
#> ⠹ [2026-06-24 04:07:11] Running for alpha [5/13] ■■■         38% | ETA:  7s
#> ✔ [2026-06-24 04:07:11] Completed 13 tasks in 6.6s
#> 
#> ℹ [2026-06-24 04:07:11] Building results
#> ! [2026-06-24 04:07:11] Found 9 failed results
#> ℹ [2026-06-24 04:07:17] ✖ Error details:
#> ℹ                       ✖ undefined columns selected (9): "delta", "gamma", "acinar" and 6 more
#> Error in `[.data.frame`(DisturbedMarkers, , "group1"): undefined columns selected
DisturbedMarkers <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
)
#> Error in UseMethod("filter"): no applicable method for 'filter' applied to an object of class "NULL"
ht5 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = DisturbedMarkers$gene,
  feature_split = DisturbedMarkers$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> Error: object 'DisturbedMarkers' not found
ht5$plot
#> Error: object 'ht5' not found

gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': object 'DisturbedMarkers' not found
DisturbedMarkers_specific <- DisturbedMarkers[
  DisturbedMarkers$gene %in% gene_specific,
]
#> Error: object 'DisturbedMarkers' not found
ht6 <- GroupHeatmap(
  srt = panc8_sub,
  layer = "data",
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> Error: object 'DisturbedMarkers_specific' not found
ht6$plot
#> Error: object 'ht6' not found

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
#> Error: object 'DisturbedMarkers_specific' not found
ht7$plot
#> Error: object 'ht7' not found

cell_index <- ave(
  seq_along(pancreas_sub$CellType),
  pancreas_sub$CellType,
  FUN = seq_along
)
pancreas_sub[["sample"]] <- paste0(
  "S",
  (cell_index - 1) %% 4 + 1
)
pancreas_sub[["condition"]] <- ifelse(
  pancreas_sub$sample %in% c("S1", "S2"),
  "ctrl",
  "case"
)
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "limma",
  fc.threshold = 1,
  layer = "counts",
  only.pos = FALSE
)
#> ℹ [2026-06-24 04:07:18] Start sample-level differential testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-06-24 04:07:20] Sample-level differential testing completed
DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  test.use = "limma",
  group_use = "Ductal",
  plot_type = "volcano",
  x_metric = "avg_log2FC",
  y_metric = "p_val"
)


pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "edgeR",
  layer = "counts",
  fc.threshold = 1,
  only.pos = FALSE
)
#> ℹ [2026-06-24 04:07:21] Start sample-level differential testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-06-24 04:07:24] Sample-level differential testing completed
DEtestPlot(
  pancreas_sub,
  group.by = "CellType",
  test.use = "edgeR",
  group_use = "Ductal",
  plot_type = "volcano",
  x_metric = "avg_log2FC",
  y_metric = "p_val",
  DE_threshold = "abs(avg_log2FC) > log2(1.5) & p_val < 0.05"
)


data(islet_bulk)
bulk_out <- RunDEtest(
  islet_bulk,
  condition_col = "condition",
  group1 = "control",
  group2 = "bfa",
  test.use = "edgeR",
  only.pos = FALSE,
  fc.threshold = 1
)
#> calcNormFactors has been renamed to normLibSizes
DEtestPlot(
  bulk_out,
  test.use = "edgeR",
  plot_type = "volcano",
  x_metric = "avg_log2FC",
  y_metric = "p_val"
)
```
