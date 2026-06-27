# Annotate single cells using SingleR

Annotate single cells using SingleR

## Usage

``` r
RunSingleR(
  srt_query,
  srt_ref,
  query_group = NULL,
  ref_group = NULL,
  query_assay = "RNA",
  ref_assay = "RNA",
  genes = "de",
  de.method = "wilcox",
  sd.thresh = 1,
  de.n = NULL,
  aggr.ref = FALSE,
  aggr.args = list(),
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  prune = TRUE,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- query_group:

  A character vector specifying the column name in the `srt_query`
  metadata that represents the cell grouping.

- ref_group:

  A character vector specifying the column name in the `srt_ref`
  metadata that represents the cell grouping.

- query_assay:

  A character vector specifying the assay to be used for the query data.
  Default is the default assay of the `srt_query` object.

- ref_assay:

  A character vector specifying the assay to be used for the reference
  data. Default is the default assay of the `srt_ref` object.

- genes:

  A string containing `"de"`, indicating that markers should be
  calculated from `ref`. For back compatibility, other string values are
  allowed but will be ignored with a deprecation warning.

  Alternatively, if `ref` is *not* a list, `genes` can be either:

  - A list of lists of character vectors containing DE genes between
    pairs of labels.

  - A list of character vectors containing marker genes for each label.

  If `ref` *is* a list, `genes` can be a list of length equal to `ref`.
  Each element of the list should be one of the two above choices
  described for non-list `ref`, containing markers for labels in the
  corresponding entry of `ref`.

- de.method:

  String specifying how DE genes should be detected between pairs of
  labels. Defaults to `"classic"`, which sorts genes by the log-fold
  changes and takes the top `de.n`. Other options are `"wilcox"` and
  `"t"`, see Details. Ignored if `genes` is a list of markers/DE genes.

- sd.thresh:

  Deprecated and ignored.

- de.n:

  An integer scalar specifying the number of DE genes to use when
  `genes="de"`. If `de.method="classic"`, defaults to
  `500 * (2/3) ^ log2(N)` where `N` is the number of unique labels.
  Otherwise, defaults to 10. Ignored if `genes` is a list of markers/DE
  genes.

- aggr.ref, aggr.args:

  Arguments controlling the aggregation of the references prior to
  annotation, see
  [`trainSingleR`](https://rdrr.io/pkg/SingleR/man/trainSingleR.html).

- quantile:

  "quantile" parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- fine.tune:

  `"fine.tune"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- tune.thresh:

  `"tune.thresh"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- prune:

  `"prune"` parameter in
  [SingleR::SingleR](https://rdrr.io/pkg/SingleR/man/SingleR.html)
  function.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

An annotate `Seurat` object. The annotation results are stored in the
`singler_annotation` column of the meta data, and the corresponding
scores are stored in the `singler_score` column.

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunKNNMap](https://mengxu98.github.io/scop/reference/RunKNNMap.md)

## Examples

``` r
data(panc8_sub)
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
#> ℹ [2026-06-27 20:26:37] Rename features for the assay: RNA
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-06-27 20:26:38] Start standard processing workflow...
#> ℹ [2026-06-27 20:26:38] Checking a list of <Seurat>...
#> ! [2026-06-27 20:26:38] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-27 20:26:38] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 20:26:38] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 20:26:38] Use the separate HVF from `srt_list`
#> ℹ [2026-06-27 20:26:38] Number of available HVF: 2000
#> ℹ [2026-06-27 20:26:38] Finished check
#> ℹ [2026-06-27 20:26:38] Perform `ScaleData()`
#> ℹ [2026-06-27 20:26:38] Perform pca linear dimension reduction
#> ℹ [2026-06-27 20:26:39] Use stored estimated dimensions 1:27 for Standardpca
#> ℹ [2026-06-27 20:26:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-27 20:26:40] Reorder clusters...
#> ℹ [2026-06-27 20:26:40] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-27 20:26:40] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-27 20:26:47] Standard processing workflow completed

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-06-27 20:26:47] Start standard processing workflow...
#> ℹ [2026-06-27 20:26:48] Checking a list of <Seurat>...
#> ! [2026-06-27 20:26:48] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-27 20:26:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 20:26:48] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-27 20:26:48] Use the separate HVF from `srt_list`
#> ℹ [2026-06-27 20:26:48] Number of available HVF: 2000
#> ℹ [2026-06-27 20:26:48] Finished check
#> ℹ [2026-06-27 20:26:48] Perform `ScaleData()`
#> ℹ [2026-06-27 20:26:48] Perform pca linear dimension reduction
#> ℹ [2026-06-27 20:26:49] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-27 20:26:49] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-27 20:26:49] Reorder clusters...
#> ℹ [2026-06-27 20:26:50] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-27 20:26:50] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-27 20:26:56] Standard processing workflow completed
pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "Standardpca_SNN_res.0.6",
  ref_group = "celltype"
)
#> ℹ [2026-06-27 20:26:56] Start SingleR annotation
#> ℹ [2026-06-27 20:27:44] Data type is log-normalized
#> ℹ [2026-06-27 20:27:44] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-06-27 20:27:44] Data type is log-normalized
#> ℹ [2026-06-27 20:27:44] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-06-27 20:27:45] Perform "SingleRCluster"
#> Detected a large SingleCellExperiment as the reference dataset, consider
#> setting 'aggr.ref = TRUE' for speed in trainSingleR(). If you know better, this
#> hint can be disabled with 'hint.sce=FALSE'.
#> ✔ [2026-06-27 20:27:46] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "SubCellType")
)


pancreas_sub <- RunSingleR(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = NULL,
  ref_group = "celltype"
)
#> ℹ [2026-06-27 20:27:46] Start SingleR annotation
#> ℹ [2026-06-27 20:27:47] Data type is log-normalized
#> ℹ [2026-06-27 20:27:47] Detected `srt_query` data type: "log_normalized_counts"
#> ℹ [2026-06-27 20:27:48] Data type is log-normalized
#> ℹ [2026-06-27 20:27:48] Detected `srt_ref` data type: "log_normalized_counts"
#> ℹ [2026-06-27 20:27:48] Perform "SingleRCell"
#> Detected a large SingleCellExperiment as the reference dataset, consider
#> setting 'aggr.ref = TRUE' for speed in trainSingleR(). If you know better, this
#> hint can be disabled with 'hint.sce=FALSE'.
#> ✔ [2026-06-27 20:27:51] SingleR annotation completed
CellDimPlot(
  pancreas_sub,
  group.by = c("singler_annotation", "SubCellType"),
  label = TRUE
)


FeatureStatPlot(
  pancreas_sub,
  stat.by = "singler_score",
  group.by = "singler_annotation"
)
#> ! [2026-06-27 20:27:52] Removed 2 groups with < 2 observations for violin plot: "sp--gp-activated-stellate" and "sp--gp-macrophage"
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's colour values.


ht1 <- CellCorHeatmap(
  srt_query = pancreas_sub,
  srt_ref = pancreas_sub,
  query_group = "SubCellType",
  cluster_rows = TRUE,
  ref_group = "singler_annotation",
  cluster_columns = TRUE,
  width = 2,
  height = 2
)
#> ℹ [2026-06-27 20:27:52] Drop [1] 19 cells with NA in the ref_group
#> ℹ [2026-06-27 20:27:53] Use the HVF to calculate distance metric
#> ℹ [2026-06-27 20:27:53] Use [1] 2000 features to calculate distance.
#> ℹ [2026-06-27 20:27:53] Detected query data type: "log_normalized_counts"
#> ℹ [2026-06-27 20:27:53] Detected reference data type: "log_normalized_counts"
#> ℹ [2026-06-27 20:27:53] Calculate similarity...
#> ℹ [2026-06-27 20:27:53] Use raw method to find neighbors
#> ℹ [2026-06-27 20:27:53] Predict cell type...
ht1$plot
```
