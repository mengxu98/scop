# Calculates dynamic features for lineages

Calculates dynamic features for lineages

## Usage

``` r
RunDynamicFeatures(
  srt,
  lineages,
  features = NULL,
  suffix = lineages,
  n_candidates = 1000,
  minfreq = 5,
  family = NULL,
  layer = "counts",
  assay = NULL,
  libsize = NULL,
  fit_method = c("gam", "pretsa"),
  knot = 0,
  max_knot_allowed = 10,
  padjust_method = "fdr",
  cores = 1,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- lineages:

  A character vector specifying the lineage names for which dynamic
  features should be calculated.

- features:

  A character vector of features to use. If `NULL`, n_candidates must be
  provided.

- suffix:

  A character vector specifying the suffix to append to the output layer
  names for each lineage. Default is the lineage names.

- n_candidates:

  A number of candidate features to select when features is `NULL`.
  Default is `1000`.

- minfreq:

  An integer specifying the minimum frequency threshold for candidate
  features. Features with a frequency less than minfreq will be
  excluded. Default is `5`.

- family:

  A character or character vector specifying the family of distributions
  to use for the GAM. If family is set to NULL, the appropriate family
  will be automatically determined based on the data. If length(family)
  is 1, the same family will be used for all features. Otherwise, family
  must have the same length as features.

- layer:

  Which layer to use. Default is `"counts"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- libsize:

  A numeric or numeric vector specifying the library size correction
  factors for each cell. If NULL, the library size correction factors
  will be calculated based on the expression matrix. If length(libsize)
  is 1, the same value will be used for all cells. Otherwise, libsize
  must have the same length as the number of cells in srt. Default is
  `NULL`.

- fit_method:

  The method used for fitting features. Either `"gam"` (generalized
  additive models) or `"pretsa"` (Pattern recognition in Temporal and
  Spatial Analyses). Default is `"gam"`.

- knot:

  For `fit_method = "pretsa"`: B-spline knots. `0` or `"auto"`. Default
  is `0`.

- max_knot_allowed:

  For `fit_method = "pretsa"` when `knot = "auto"`: max knots. Default
  is `10`.

- padjust_method:

  The method used for p-value adjustment. Default is `"fdr"`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

Returns the modified Seurat object with the calculated dynamic features
stored in the tools slot.

## References

Zhuang, H., Ji, Z. PreTSA: computationally efficient modeling of
temporal and spatial gene expression patterns. Genome Biol (2026).
https://doi.org/10.1186/s13059-026-03994-3

## See also

[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md),
[DynamicPlot](https://mengxu98.github.io/scop/reference/DynamicPlot.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-04-02 16:39:47] Start standard processing workflow...
#> ℹ [2026-04-02 16:39:48] Checking a list of <Seurat>...
#> ! [2026-04-02 16:39:48] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:39:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:39:50] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:39:50] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:39:50] Number of available HVF: 2000
#> ℹ [2026-04-02 16:39:50] Finished check
#> ℹ [2026-04-02 16:39:51] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:39:51] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:39:55] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:39:55] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:39:55] Reorder clusters...
#> ℹ [2026-04-02 16:39:55] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:39:55] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:39:55] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:39:55] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:39:55] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:39:59] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:40:02] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Error in loadNamespace(x): there is no package called ‘slingshot’

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "gam"
)
#> ℹ [2026-04-02 16:40:06] Start find dynamic features
#> ℹ [2026-04-02 16:40:08] Data type is raw counts
#> Error in subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])]): No cells found

names(
  pancreas_sub@tools$DynamicFeatures_Lineage1
)
#> NULL
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#> NULL
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> Error in DynamicHeatmap(pancreas_sub, lineages = c("Lineage1", "Lineage2"),     cell_annotation = "SubCellType", n_split = 3, reverse_ht = "Lineage1"): Lineages: Lineage1 is not in the meta data of <Seurat>
ht$plot
#> Error: object 'ht' not found

DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ℹ [2026-04-02 16:40:08] Start find dynamic features
#> ℹ [2026-04-02 16:40:09] Data type is raw counts
#> Error in subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])]): No cells found

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "pretsa"
)
#> ℹ [2026-04-02 16:40:09] Start find dynamic features
#> ℹ [2026-04-02 16:40:10] Data type is raw counts
#> Error in subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])]): No cells found
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#> NULL
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> Error in DynamicHeatmap(pancreas_sub, lineages = c("Lineage1", "Lineage2"),     cell_annotation = "SubCellType", n_split = 3, reverse_ht = "Lineage1"): Lineages: Lineage1 is not in the meta data of <Seurat>
ht$plot
#> Error: object 'ht' not found

DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ℹ [2026-04-02 16:40:10] Start find dynamic features
#> ℹ [2026-04-02 16:40:11] Data type is raw counts
#> Error in subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])]): No cells found
```
