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
#> ℹ [2026-04-21 07:34:50] Start standard processing workflow...
#> ℹ [2026-04-21 07:34:51] Checking a list of <Seurat>...
#> ! [2026-04-21 07:34:51] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-21 07:34:51] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:34:53] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 07:34:54] Use the separate HVF from `srt_list`
#> ℹ [2026-04-21 07:34:54] Number of available HVF: 2000
#> ℹ [2026-04-21 07:34:54] Finished check
#> ℹ [2026-04-21 07:34:54] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-21 07:34:55] Perform pca linear dimension reduction
#> ℹ [2026-04-21 07:34:55] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-21 07:34:56] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-21 07:34:56] Reorder clusters...
#> ℹ [2026-04-21 07:34:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-21 07:34:56] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-21 07:34:56] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-21 07:35:00] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-21 07:35:05] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_path()`).


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "gam"
)
#> ℹ [2026-04-21 07:35:07] Start find dynamic features
#> ℹ [2026-04-21 07:35:08] Data type is raw counts
#> ℹ [2026-04-21 07:35:10] Number of candidate features (union): 236
#> ℹ [2026-04-21 07:35:11] Data type is raw counts
#> ℹ [2026-04-21 07:35:11] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-21 07:35:11] Using 1 core
#> ⠙ [2026-04-21 07:35:11] Running for Gcg [1/236]              0% | ETA:  9s
#> ⠹ [2026-04-21 07:35:11] Running for Ttr [23/236]             10% | ETA:  6s
#> ⠸ [2026-04-21 07:35:11] Running for Mest [66/236] ■■          28% | ETA:  9s
#> ⠼ [2026-04-21 07:35:11] Running for Rap1b [152/236] ■■■■■■      64% | ETA:  4s
#> ⠴ [2026-04-21 07:35:11] Running for Abcc8 [230/236] ■■■■■■■■■   97% | ETA:  0s
#> ✔ [2026-04-21 07:35:11] Completed 236 tasks in 9.9s
#> 
#> ℹ [2026-04-21 07:35:11] Building results
#> ℹ [2026-04-21 07:35:21] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-21 07:35:21] Using 1 core
#> ⠙ [2026-04-21 07:35:21] Running for Irs4 [61/236] ■■          26% | ETA: 10s
#> ⠹ [2026-04-21 07:35:21] Running for Spc25 [117/236] ■■■■        50% | ETA:  6s
#> ⠸ [2026-04-21 07:35:21] Running for Ctxn2 [190/236] ■■■■■■■■    81% | ETA:  2s
#> ✔ [2026-04-21 07:35:21] Completed 236 tasks in 10.4s
#> 
#> ℹ [2026-04-21 07:35:21] Building results
#> ✔ [2026-04-21 07:35:31] Find dynamic features done

names(
  pancreas_sub@tools$DynamicFeatures_Lineage1
)
#> [1] "DynamicFeatures" "raw_matrix"      "fitted_matrix"   "upr_matrix"     
#> [5] "lwr_matrix"      "libsize"         "lineages"        "family"         
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime  valleytime pvalue
#> Gcg       Gcg        170 0.6058287 0.7774265 22.79123 14.27994337      0
#> Ghrl     Ghrl        160 0.3614421 0.6994200 19.39321 13.28320185      0
#> Iapp     Iapp        261 0.3110362 0.7734046 22.79123  0.07106129      0
#> Pyy       Pyy        394 0.4063756 0.7946711 20.69916  0.07106129      0
#> Rbp4     Rbp4        352 0.4421445 0.7584099 19.87018 12.35433740      0
#> Chgb     Chgb        184 0.2311685 0.6920917 21.26291  0.07106129      0
#>      padjust
#> Gcg        0
#> Ghrl       0
#> Iapp       0
#> Pyy        0
#> Rbp4       0
#> Chgb       0
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-04-21 07:35:31] [1] 183 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Chgb,Lrpprc,Slc38a5,Cck,Cdkn1a...
#> ℹ [2026-04-21 07:35:32] 
#> ℹ                       The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ                       The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ                       If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht$plot


DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ℹ [2026-04-21 07:35:35] Start find dynamic features
#> ℹ [2026-04-21 07:35:36] Data type is raw counts
#> ℹ [2026-04-21 07:35:37] Number of candidate features (union): 2
#> ℹ [2026-04-21 07:35:37] Data type is raw counts
#> ℹ [2026-04-21 07:35:37] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-21 07:35:37] Using 1 core
#> ⠙ [2026-04-21 07:35:37] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-21 07:35:37] Completed 2 tasks in 166ms
#> 
#> ℹ [2026-04-21 07:35:37] Building results
#> ✔ [2026-04-21 07:35:37] Find dynamic features done
#> ℹ [2026-04-21 07:35:37] Start find dynamic features
#> ℹ [2026-04-21 07:35:39] Data type is raw counts
#> ℹ [2026-04-21 07:35:39] Number of candidate features (union): 2
#> ℹ [2026-04-21 07:35:40] Data type is raw counts
#> ℹ [2026-04-21 07:35:40] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-21 07:35:40] Using 1 core
#> ⠙ [2026-04-21 07:35:40] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-21 07:35:40] Completed 2 tasks in 190ms
#> 
#> ℹ [2026-04-21 07:35:40] Building results
#> ✔ [2026-04-21 07:35:40] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "pretsa"
)
#> ℹ [2026-04-21 07:35:41] Start find dynamic features
#> ℹ [2026-04-21 07:35:42] Data type is raw counts
#> ℹ [2026-04-21 07:35:44] Number of candidate features (union): 236
#> ℹ [2026-04-21 07:35:44] Data type is raw counts
#> ℹ [2026-04-21 07:35:44] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-21 07:35:44] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-04-21 07:35:44] Find dynamic features done
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime  valleytime        pvalue
#> Gcg       Gcg        170 0.5593818 0.5593818 22.79123 13.44902051 1.820103e-116
#> Ghrl     Ghrl        160 0.1503809 0.1503809 18.70697  5.16726925  4.504737e-23
#> Iapp     Iapp        261 0.6564627 0.6564627 22.79123  0.07106129 6.131719e-152
#> Pyy       Pyy        394 0.7021626 0.7021626 22.79123  8.15806399 2.734634e-172
#> Rbp4     Rbp4        352 0.6342685 0.6342685 22.79123  7.97512852 5.147270e-143
#> Chgb     Chgb        184 0.5134108 0.5134108 22.79123  9.01545802 2.511481e-102
#>            padjust
#> Gcg  6.818162e-116
#> Ghrl  6.365976e-23
#> Iapp 4.989951e-151
#> Pyy  3.396703e-171
#> Rbp4 3.470730e-142
#> Chgb 7.056064e-102
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-04-21 07:35:45] [1] 172 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Chgb,Gast,Lrpprc,Slc38a5,Cck,Cdkn1a...
#> ℹ [2026-04-21 07:35:46] 
#> ℹ                       The size of the heatmap is fixed because certain elements are not scalable.
#> ℹ                       The width and height of the heatmap are determined by the size of the current viewport.
#> ℹ                       If you want to have more control over the size, you can manually set the parameters 'width' and 'height'.

ht$plot


DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ℹ [2026-04-21 07:35:48] Start find dynamic features
#> ℹ [2026-04-21 07:35:49] Data type is raw counts
#> ℹ [2026-04-21 07:35:50] Number of candidate features (union): 2
#> ℹ [2026-04-21 07:35:50] Data type is raw counts
#> ℹ [2026-04-21 07:35:51] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-21 07:35:51] Using 1 core
#> ⠙ [2026-04-21 07:35:51] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-21 07:35:51] Completed 2 tasks in 162ms
#> 
#> ℹ [2026-04-21 07:35:51] Building results
#> ✔ [2026-04-21 07:35:51] Find dynamic features done
#> ℹ [2026-04-21 07:35:51] Start find dynamic features
#> ℹ [2026-04-21 07:35:52] Data type is raw counts
#> ℹ [2026-04-21 07:35:53] Number of candidate features (union): 2
#> ℹ [2026-04-21 07:35:53] Data type is raw counts
#> ℹ [2026-04-21 07:35:53] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-21 07:35:53] Using 1 core
#> ℹ [2026-04-21 07:35:53] Building results
#> ✔ [2026-04-21 07:35:53] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
