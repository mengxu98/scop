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
  min.cells.sample = 10,
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
  min.cells.sample = 10,
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

  Additional arguments passed to the marker backend. Arguments not
  supported by scop's sparse Wilcoxon path are delegated to
  [Seurat::FindMarkers](https://satijalab.org/seurat/reference/FindMarkers.html).

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

  Feature type used for differential testing.

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
  average fold change between/among the two groups.

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
  For a raw `ChromatinAssay`, `RunDEtest()` first applies
  [`Signac::RunTFIDF()`](https://stuartlab.org/signac/reference/RunTFIDF.html);
  other assays default to `"LogNormalize"`.

- sample_col:

  Metadata column storing biological sample IDs. Required when
  `test.use` is a sample-level pseudobulk method on `Seurat`.

- condition_col:

  Metadata column storing condition labels. Required when `test.use` is
  a sample-level pseudobulk method on `Seurat`, and required for
  `SummarizedExperiment` input.

- min.cells.sample:

  Minimum number of cells required for a sample to be kept in
  sample-level pseudobulk testing. Applied independently within each
  `group.by` level (for example each cell type). Samples below the
  threshold are dropped before aggregation. Default `10`. Set to `1` to
  keep every sample that has at least one cell. Ignored for cell-level
  tests and for `SummarizedExperiment` input.

- bulk_assay:

  Assay name used as the bulk counts matrix for `SummarizedExperiment`
  input.

- p.adjust.method:

  A character value specifying the method to use for adjusting p-values.

- layer:

  Assay layer to use.

- assay:

  Assay to use in differential expression testing

- seed:

  Optional integer seed. When supplied, every input receives a
  deterministic independent L'Ecuyer-CMRG random-number stream, making
  results reproducible across worker counts and scheduling order. The
  caller's random number state is restored when the call finishes.

- verbose:

  Whether to print the message. Default is `TRUE`.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

## See also

[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:54:24] Start standard processing workflow...
#> ℹ [2026-09-06 21:54:25] Checking a list of <Seurat>...
#> ! [2026-09-06 21:54:25] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:54:25] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:54:25] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:54:25] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:54:25] Number of available HVF: 2000
#> ℹ [2026-09-06 21:54:25] Finished check
#> ℹ [2026-09-06 21:54:25] Perform `ScaleData()`
#> ℹ [2026-09-06 21:54:25] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:54:25] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:54:26] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:54:26] Reorder clusters...
#> ℹ [2026-09-06 21:54:26] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:54:26] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:54:32] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType",
  only.pos = FALSE
)
#> ℹ [2026-09-06 21:54:32] Data type is log-normalized
#> ℹ [2026-09-06 21:54:33] Start differential expression test
#> ℹ [2026-09-06 21:54:33] Find all markers(wilcox) among [1] 8 groups...
#> ℹ [2026-09-06 21:54:33] Using 1 core
#> ⠙ [2026-09-06 21:54:33] Running for Ductal [1/8] ■           12% | ETA:  1m
#> ⠹ [2026-09-06 21:54:33] Running for Ngn3-high-EP [2/8] ■■          25% | ETA:  …
#> ⠸ [2026-09-06 21:54:33] Running for Beta [3/8] ■■■         38% | ETA: 46s
#> ⠼ [2026-09-06 21:54:33] Running for Ngn3-low-EP [4/8] ■■■■■       50% | ETA: 36s
#> ⠴ [2026-09-06 21:54:33] Running for Pre-endocrine [5/8] ■■■■■■      62% | ETA: …
#> ⠦ [2026-09-06 21:54:33] Running for Alpha [6/8] ■■■■■■■     75% | ETA: 17s
#> ⠧ [2026-09-06 21:54:33] Running for Epsilon [7/8] ■■■■■■■■    88% | ETA:  8s
#> ✔ [2026-09-06 21:54:33] Completed 8 tasks in 1m 12.6s
#> 
#> ℹ [2026-09-06 21:54:33] Building results
#> ✔ [2026-09-06 21:55:45] Differential expression test completed
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
#> ℹ [2026-09-06 21:55:46] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:55:46] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:55:46] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
#> ℹ [2026-09-06 21:55:49] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:55:49] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:55:49] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht2$plot


pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "SubCellType",
  markers_type = "paired",
  cores = 2
)
#> ℹ [2026-09-06 21:55:50] Data type is log-normalized
#> ℹ [2026-09-06 21:55:50] Start differential expression test
#> ℹ [2026-09-06 21:55:50] Find paired markers(wilcox) among [1] 8 groups...
#> ℹ [2026-09-06 21:55:50] Using 2 cores
#> ⠙ [2026-09-06 21:55:50] Running for 5 [1/28]              4% | ETA: 23m
#> ⠹ [2026-09-06 21:55:50] Running for 1 [4/28] ■           14% | ETA:  5m
#> ⠸ [2026-09-06 21:55:50] Running for 12 [8/28] ■■          29% | ETA:  3m
#> ⠼ [2026-09-06 21:55:50] Running for 15 [15/28] ■■■■■       54% | ETA:  1m
#> ⠴ [2026-09-06 21:55:50] Running for 22 [18/28] ■■■■■■      64% | ETA: 46s
#> ⠦ [2026-09-06 21:55:50] Running for 18 [21/28] ■■■■■■■     75% | ETA: 28s
#> ⠧ [2026-09-06 21:55:50] Running for 25 [25/28] ■■■■■■■■    89% | ETA: 11s
#> ✔ [2026-09-06 21:55:50] Completed 28 tasks in 1m 32.4s
#> 
#> ℹ [2026-09-06 21:55:50] Building results
#> ✔ [2026-09-06 21:57:22] Differential expression test completed
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
#> ℹ [2026-09-06 21:57:28] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:57:28] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:57:28] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht3$plot


data(panc8_sub)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-09-06 21:57:32] Run integration workflow...
#> ℹ [2026-09-06 21:57:32] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-06 21:57:33] Checking a list of <Seurat>...
#> ! [2026-09-06 21:57:33] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:57:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-09-06 21:57:33] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-09-06 21:57:33] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:57:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-09-06 21:57:33] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-09-06 21:57:33] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:57:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-09-06 21:57:33] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-09-06 21:57:33] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:57:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-09-06 21:57:34] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-09-06 21:57:34] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:57:34] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 21:57:34] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-06 21:57:34] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:57:34] Number of available HVF: 2000
#> ℹ [2026-09-06 21:57:34] Finished check
#> ℹ [2026-09-06 21:57:36] Perform Uncorrected integration
#> ℹ [2026-09-06 21:57:36] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-06 21:57:37] Perform "pca" linear dimension reduction
#> ℹ [2026-09-06 21:57:37] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-06 21:57:38] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-06 21:57:38] Reorder clusters...
#> ℹ [2026-09-06 21:57:38] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:57:38] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-09-06 21:57:43] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-09-06 21:57:48] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ✔ [2026-09-06 21:57:54] Uncorrected integration completed
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
#> ℹ [2026-09-06 21:57:55] Data type is log-normalized
#> ℹ [2026-09-06 21:57:55] Start differential expression test
#> ℹ [2026-09-06 21:57:55] Find conserved markers(wilcox) among [1] 13 groups...
#> ℹ [2026-09-06 21:57:55] Using 2 cores
#> ⠙ [2026-09-06 21:57:55] Running for ductal [1/13]              8% | ETA:  9m
#> ⠹ [2026-09-06 21:57:55] Running for acinar [2/13] ■           15% | ETA:  5m
#> ⠸ [2026-09-06 21:57:55] Running for beta [4/13] ■■■         31% | ETA:  3m
#> ⠼ [2026-09-06 21:57:55] Running for macrophage [6/13] ■■■■        46% | ETA:  1m
#> ⠴ [2026-09-06 21:57:55] Running for activated-stellate [7/13] ■■■■■       54% |…
#> ⠦ [2026-09-06 21:57:55] Running for endothelial [10/13] ■■■■■■■     77% | ETA: …
#> ⠧ [2026-09-06 21:57:55] Running for mast [12/13] ■■■■■■■■■   92% | ETA: 10s
#> ✔ [2026-09-06 21:57:55] Completed 13 tasks in 1m 60s
#> 
#> ℹ [2026-09-06 21:57:55] Building results
#> ✔ [2026-09-06 21:59:55] Differential expression test completed
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
#> ℹ [2026-09-06 21:59:59] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 21:59:59] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 21:59:59] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht4$plot


panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "tech",
  grouping.var = "celltype",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-09-06 22:00:02] Data type is log-normalized
#> ℹ [2026-09-06 22:00:02] Start differential expression test
#> ℹ [2026-09-06 22:00:02] Find conserved markers(wilcox) among [1] 5 groups...
#> ℹ [2026-09-06 22:00:02] Using 2 cores
#> ⠙ [2026-09-06 22:00:02] Running for celseq [1/5] ■■          20% | ETA:  3m
#> ⠹ [2026-09-06 22:00:02] Running for celseq2 [2/5] ■■■■        40% | ETA:  1m
#> ⠸ [2026-09-06 22:00:02] Running for smartseq2 [3/5] ■■■■■■      60% | ETA: 49s
#> ⠼ [2026-09-06 22:00:02] Running for indrop [4/5] ■■■■■■■■    80% | ETA: 22s
#> ✔ [2026-09-06 22:00:02] Completed 5 tasks in 1m 29.1s
#> 
#> ℹ [2026-09-06 22:00:02] Building results
#> ✔ [2026-09-06 22:01:31] Differential expression test completed
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
#> ℹ [2026-09-06 22:01:32] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 22:01:32] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 22:01:32] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht4$plot


panc8_sub <- RunDEtest(
  srt = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "disturbed",
  cores = 2
)
#> ℹ [2026-09-06 22:01:37] Data type is log-normalized
#> ℹ [2026-09-06 22:01:37] Start differential expression test
#> ℹ [2026-09-06 22:01:37] Find disturbed markers(wilcox) among [1] 13 groups...
#> ℹ [2026-09-06 22:01:37] Using 2 cores
#> ⠙ [2026-09-06 22:01:37] Running for ductal [1/13]              8% | ETA:  8m
#> ⠹ [2026-09-06 22:01:37] Running for acinar [2/13] ■           15% | ETA:  5m
#> ⠸ [2026-09-06 22:01:37] Running for activated-stellate [4/13] ■■■         31% |…
#> ⠼ [2026-09-06 22:01:37] Running for endothelial [7/13] ■■■■■       54% | ETA:  …
#> ⠴ [2026-09-06 22:01:37] Running for mast [10/13] ■■■■■■■     77% | ETA: 26s
#> ⠦ [2026-09-06 22:01:37] Running for beta [12/13] ■■■■■■■■■   92% | ETA:  8s
#> ✔ [2026-09-06 22:01:37] Completed 13 tasks in 1m 36.5s
#> 
#> ℹ [2026-09-06 22:01:37] Building results
#> ✔ [2026-09-06 22:03:14] Differential expression test completed
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
#> ℹ [2026-09-06 22:03:16] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 22:03:16] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 22:03:16] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
#> ℹ [2026-09-06 22:03:31] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 22:03:31] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 22:03:31] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
#> ! [2026-09-06 22:03:42] When 'grouping.var' is specified, 'exp_method' can only be 'log2fc'
#> ℹ [2026-09-06 22:03:44] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-09-06 22:03:44] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-09-06 22:03:44] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht7$plot


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
#> ℹ [2026-09-06 22:03:45] Start sample-level differential testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-09-06 22:03:47] Sample-level differential testing completed
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
  min.cells.sample = 10,
  only.pos = FALSE
)
#> ℹ [2026-09-06 22:03:47] Start sample-level differential testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-09-06 22:03:50] Sample-level differential testing completed
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
