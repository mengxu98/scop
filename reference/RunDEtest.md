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
  ...,
  design = NULL,
  contrast = NULL
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
  ...,
  design = NULL,
  contrast = NULL
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

  Deprecated alias for `object`; it will be removed in scop 1.0.0.

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

- design:

  Optional one-sided formula for sample-level `test.use = "edgeR"`,
  evaluated using cell metadata aggregated by `sample_col`, or bulk
  `colData`. Supply together with `contrast`, using exact argument
  names. Available for RNA/gene and ATAC peak/cCRE inputs; omitted
  arguments preserve legacy behavior. Every formula variable must be
  constant within an aggregate sample. For paired cell types, use
  distinct donor-by-cell-type sample IDs and a donor block such as
  `~ donor + cell_type`. Missing combinations are never zero-filled. The
  model must have full rank and positive residual degrees of freedom,
  with at least two biological samples per selected condition. Formula
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms are not
  supported and are rejected explicitly.

- contrast:

  Optional finite named numeric coefficient contrast for `design`. Names
  must match every fitted model-matrix column exactly, including
  zero-weight columns; order is aligned by name. Positive log-fold
  changes follow this vector, not necessarily `group2` versus `group1`.
  Only one condition pair is fitted per `group.by` level. Invalid
  explicit models fail rather than being skipped.

## Details

Without `design` and `contrast`, existing defaults and comparison
behavior are unchanged. Explicit edgeR models use native
`filterByExpr(design = ...)`, TMM normalization, dispersion estimation,
`glmQLFit(robust = TRUE)`, and `glmQLFTest(contrast = ...)`. Both
`p_val_adj` (the requested adjustment) and `p_val_BH` use the full
tested feature universe before fold-change display filtering. Set
`p.adjust.method = "BH"` to also use BH for `p_val_adj`; its default
remains `"bonferroni"`. A fold-change display filter is not a
minimum-effect-size hypothesis test.

Seurat results retain one receipt per group under
`object@tools[["DEtest_pseudobulk"]]$edgeR_models`, or
`object@tools[[paste0("DEtest_", group.by)]]$edgeR_models`. Each receipt
contains `full_results`, input/tested/filtered feature IDs, raw
pseudobulk counts, sample metadata and cell membership, cell counts,
library sizes, normalization factors, design matrix, contrast,
direction, residual degrees of freedom and backend versions. Bulk
receipts are under
`S4Vectors::metadata(object)$DEtest$methods$de_edgeR_qlf$details$groups`,
in each comparison's `details$edgeR_model` (without per-cell
membership). Explicit `peak`/`cCRE` models require a `ChromatinAssay` or
`RangedSummarizedExperiment`; full and displayed results include 1-based
closed coordinates and supplied genome labels. No coordinates or count
units are guessed from feature names. Cells and technical libraries must
not be relabeled as independent biological replicates; metadata cannot
establish independence.

## See also

[VolcanoPlot](https://mengxu98.github.io/scop/reference/VolcanoPlot.md),
[RunEnrichment](https://mengxu98.github.io/scop/reference/RunEnrichment.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GroupHeatmap](https://mengxu98.github.io/scop/reference/GroupHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-10-25 21:05:36] Start standard processing workflow...
#> ℹ [2026-10-25 21:05:36] Checking a list of <Seurat>...
#> ! [2026-10-25 21:05:36] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:05:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-25 21:05:36] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-25 21:05:36] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:05:36] Number of available HVF: 2000
#> ℹ [2026-10-25 21:05:36] Finished check
#> ℹ [2026-10-25 21:05:36] Perform `ScaleData()`
#> ℹ [2026-10-25 21:05:36] Perform pca linear dimension reduction
#> ℹ [2026-10-25 21:05:37] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-25 21:05:37] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-25 21:05:37] Reorder clusters...
#> ℹ [2026-10-25 21:05:37] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:05:37] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-25 21:05:45] Standard processing workflow completed
pancreas_sub <- RunDEtest(
  object = pancreas_sub,
  group.by = "SubCellType",
  only.pos = FALSE
)
#> ℹ [2026-10-25 21:05:45] Data type is log-normalized
#> ℹ [2026-10-25 21:05:45] Start differential expression test
#> ℹ [2026-10-25 21:05:45] Find all markers(wilcox) among [1] 8 groups...
#> ✔ [2026-10-25 21:05:45] Differential expression test completed
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
#> ℹ [2026-10-25 21:05:47] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:05:47] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:05:47] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
#> ℹ [2026-10-25 21:05:48] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:05:48] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:05:48] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht2$plot


pancreas_sub <- RunDEtest(
  object = pancreas_sub,
  group.by = "SubCellType",
  markers_type = "paired",
  cores = 2
)
#> ℹ [2026-10-25 21:05:49] Data type is log-normalized
#> ℹ [2026-10-25 21:05:49] Start differential expression test
#> ℹ [2026-10-25 21:05:49] Find paired markers(wilcox) among [1] 8 groups...
#> ℹ [2026-10-25 21:05:49] Using 2 cores
#> ⠙ [2026-10-25 21:05:49] Running for 5 [1/28]              4% | ETA:  1m
#> ✔ [2026-10-25 21:05:49] Completed 28 tasks in 3.2s
#> 
#> ℹ [2026-10-25 21:05:49] Building results
#> ✔ [2026-10-25 21:05:52] Differential expression test completed
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
#> ℹ [2026-10-25 21:05:59] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:05:59] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:05:59] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht3$plot


data(panc8_sub)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-10-25 21:06:01] Run integration workflow...
#> ℹ [2026-10-25 21:06:01] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-25 21:06:01] Checking a list of <Seurat>...
#> ! [2026-10-25 21:06:02] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:06:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-10-25 21:06:02] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-10-25 21:06:02] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:06:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-10-25 21:06:02] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-10-25 21:06:02] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:06:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-10-25 21:06:02] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-10-25 21:06:02] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:06:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-10-25 21:06:02] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-10-25 21:06:02] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-25 21:06:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:06:03] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-25 21:06:03] Use the separate HVF from `srt_list`
#> ℹ [2026-10-25 21:06:03] Number of available HVF: 2000
#> ℹ [2026-10-25 21:06:03] Finished check
#> ℹ [2026-10-25 21:06:04] Perform Uncorrected integration
#> ℹ [2026-10-25 21:06:05] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-25 21:06:06] Perform "pca" linear dimension reduction
#> ℹ [2026-10-25 21:06:06] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-25 21:06:07] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-25 21:06:07] Reorder clusters...
#> ℹ [2026-10-25 21:06:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-25 21:06:07] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-10-25 21:06:13] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-10-25 21:06:19] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ✔ [2026-10-25 21:06:25] Uncorrected integration completed
CellDimPlot(
  panc8_sub,
  group.by = c("celltype", "tech")
)


panc8_sub <- RunDEtest(
  object = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-10-25 21:06:26] Data type is log-normalized
#> ℹ [2026-10-25 21:06:26] Start differential expression test
#> ℹ [2026-10-25 21:06:26] Find conserved markers(wilcox) among [1] 13 groups...
#> ℹ [2026-10-25 21:06:26] Using 2 cores
#> ⠙ [2026-10-25 21:06:26] Running for ductal [1/13]              8% | ETA:  2m
#> ⠹ [2026-10-25 21:06:26] Running for mast [10/13] ■■■■■■■     77% | ETA:  3s
#> ✔ [2026-10-25 21:06:26] Completed 13 tasks in 11.8s
#> 
#> ℹ [2026-10-25 21:06:26] Building results
#> ✔ [2026-10-25 21:06:38] Differential expression test completed
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
#> ℹ [2026-10-25 21:06:42] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:06:42] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:06:42] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht4$plot


panc8_sub <- RunDEtest(
  object = panc8_sub,
  group.by = "tech",
  grouping.var = "celltype",
  markers_type = "conserved",
  cores = 2
)
#> ℹ [2026-10-25 21:06:45] Data type is log-normalized
#> ℹ [2026-10-25 21:06:45] Start differential expression test
#> ℹ [2026-10-25 21:06:45] Find conserved markers(wilcox) among [1] 5 groups...
#> ℹ [2026-10-25 21:06:45] Using 2 cores
#> ⠙ [2026-10-25 21:06:45] Running for celseq [1/5] ■■          20% | ETA: 35s
#> ✔ [2026-10-25 21:06:45] Completed 5 tasks in 11s
#> 
#> ℹ [2026-10-25 21:06:45] Building results
#> ✔ [2026-10-25 21:06:56] Differential expression test completed
ConservedMarkers2 <- dplyr::filter(
  panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)
ht4 <- GroupHeatmap(
  object = panc8_sub,
  layer = "data",
  features = ConservedMarkers2$gene,
  feature_split = ConservedMarkers2$group1,
  group.by = "tech",
  split.by = "celltype"
)
#> ℹ [2026-10-25 21:06:57] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:06:57] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:06:57] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht4$plot


panc8_sub <- RunDEtest(
  object = panc8_sub,
  group.by = "celltype",
  grouping.var = "tech",
  markers_type = "disturbed",
  cores = 2
)
#> ℹ [2026-10-25 21:06:59] Data type is log-normalized
#> ℹ [2026-10-25 21:06:59] Start differential expression test
#> ℹ [2026-10-25 21:06:59] Find disturbed markers(wilcox) among [1] 13 groups...
#> ℹ [2026-10-25 21:06:59] Using 2 cores
#> ⠙ [2026-10-25 21:06:59] Running for ductal [1/13]              8% | ETA:  3m
#> ⠹ [2026-10-25 21:06:59] Running for beta [4/13] ■■■         31% | ETA: 35s
#> ✔ [2026-10-25 21:06:59] Completed 13 tasks in 17s
#> 
#> ℹ [2026-10-25 21:06:59] Building results
#> ✔ [2026-10-25 21:07:16] Differential expression test completed
DisturbedMarkers <- dplyr::filter(
  panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2"
)
ht5 <- GroupHeatmap(
  object = panc8_sub,
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
#> ℹ [2026-10-25 21:07:18] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:07:18] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:07:18] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht5$plot


gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
DisturbedMarkers_specific <- DisturbedMarkers[
  DisturbedMarkers$gene %in% gene_specific,
]
ht6 <- GroupHeatmap(
  object = panc8_sub,
  layer = "data",
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  split.by = "tech"
)
#> ℹ [2026-10-25 21:07:22] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:07:22] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:07:22] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
ht6$plot


ht7 <- GroupHeatmap(
  object = panc8_sub,
  layer = "data",
  aggregate_fun = function(x) mean(expm1(x)) + 1,
  features = DisturbedMarkers_specific$gene,
  feature_split = DisturbedMarkers_specific$group1,
  group.by = "celltype",
  grouping.var = "tech",
  numerator = "smartseq2"
)
#> ! [2026-10-25 21:07:24] When 'grouping.var' is specified, 'exp_method' can only be 'log2fc'
#> ℹ [2026-10-25 21:07:25] The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ [2026-10-25 21:07:26] The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ [2026-10-25 21:07:26] If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.
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
  object = pancreas_sub,
  group.by = "CellType",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "limma",
  fc.threshold = 1,
  layer = "counts",
  only.pos = FALSE
)
#> ℹ [2026-10-25 21:07:26] Start sample-level differential testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-10-25 21:07:28] Sample-level differential testing completed
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
  object = pancreas_sub,
  group.by = "CellType",
  sample_col = "sample",
  condition_col = "condition",
  test.use = "edgeR",
  layer = "counts",
  fc.threshold = 1,
  min.cells.sample = 10,
  only.pos = FALSE
)
#> ℹ [2026-10-25 21:07:28] Start sample-level differential testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> ✔ [2026-10-25 21:07:31] Sample-level differential testing completed
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
  object = islet_bulk,
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


# Explicit paired sample-level model (synthetic biological-sample counts).
if (requireNamespace("edgeR", quietly = TRUE) &&
  requireNamespace("statmod", quietly = TRUE)) {
  local({
    withr::local_preserve_seed()
    set.seed(11)
    samples <- data.frame(
      donor = factor(rep(paste0("D", 1:4), each = 2)),
      condition = factor(rep(c("control", "treated"), 4)),
      row.names = paste0("sample", 1:8)
    )
    counts <- matrix(rnbinom(200 * 8, mu = 60, size = 12),
      nrow = 200,
      dimnames = list(paste0("gene", 1:200), rownames(samples))
    )
    bulk <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts), colData = samples
    )
    columns <- colnames(model.matrix(~ donor + condition, samples))
    contrast <- setNames(rep(0, length(columns)), columns)
    contrast["conditiontreated"] <- 1
    fitted <- RunDEtest(bulk,
      condition_col = "condition",
      group1 = "control", group2 = "treated", test.use = "edgeR",
      design = ~ donor + condition, contrast = contrast,
      only.pos = FALSE, fc.threshold = 1, p.adjust.method = "BH", verbose = FALSE
    )
    results <- S4Vectors::metadata(fitted)$DEtest$results
    stopifnot(all(results$p_val_BH >= 0 & results$p_val_BH <= 1))
  })
}
#> calcNormFactors has been renamed to normLibSizes
# End explicit paired sample-level model.
```
