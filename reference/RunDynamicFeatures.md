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
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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
#> ℹ [2026-04-27 16:43:20] Start standard processing workflow...
#> ℹ [2026-04-27 16:43:21] Checking a list of <Seurat>...
#> ! [2026-04-27 16:43:21] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-27 16:43:21] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-27 16:43:24] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-27 16:43:24] Use the separate HVF from `srt_list`
#> ℹ [2026-04-27 16:43:24] Number of available HVF: 2000
#> ℹ [2026-04-27 16:43:25] Finished check
#> ℹ [2026-04-27 16:43:25] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-27 16:43:25] Perform pca linear dimension reduction
#> ℹ [2026-04-27 16:43:26] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-27 16:43:26] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-27 16:43:26] Reorder clusters...
#> ℹ [2026-04-27 16:43:26] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-27 16:43:26] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-27 16:43:26] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-27 16:43:31] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-27 16:43:36] Standard processing workflow completed
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
#> ℹ [2026-04-27 16:43:39] Start find dynamic features
#> ℹ [2026-04-27 16:43:40] Data type is raw counts
#> ℹ [2026-04-27 16:43:42] Number of candidate features (union): 236
#> ℹ [2026-04-27 16:43:42] Data type is raw counts
#> ℹ [2026-04-27 16:43:42] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-27 16:43:42] Using 1 core
#> ⠙ [2026-04-27 16:43:42] Running for Gcg [1/236]              0% | ETA:  9s
#> ⠹ [2026-04-27 16:43:42] Running for Cdca8 [76/236] ■■■         32% | ETA:  5s
#> ⠸ [2026-04-27 16:43:42] Running for Lig1 [162/236] ■■■■■■      69% | ETA:  3s
#> ✔ [2026-04-27 16:43:42] Completed 236 tasks in 8.1s
#> 
#> ℹ [2026-04-27 16:43:42] Building results
#> ℹ [2026-04-27 16:43:50] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-27 16:43:50] Using 1 core
#> ⠙ [2026-04-27 16:43:50] Running for Slc38a5 [10/236]              4% | ETA:  7s
#> ⠹ [2026-04-27 16:43:50] Running for Irs4 [61/236] ■■          26% | ETA: 11s
#> ⠸ [2026-04-27 16:43:50] Running for Mki67 [106/236] ■■■■        45% | ETA:  8s
#> ⠼ [2026-04-27 16:43:50] Running for Pou3f1 [157/236] ■■■■■■      67% | ETA:  5s
#> ⠴ [2026-04-27 16:43:50] Running for G6pc2 [213/236] ■■■■■■■■■   90% | ETA:  1s
#> ✔ [2026-04-27 16:43:50] Completed 236 tasks in 13.7s
#> 
#> ℹ [2026-04-27 16:43:50] Building results
#> ✔ [2026-04-27 16:44:04] Find dynamic features done

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
#> ℹ [2026-04-27 16:44:04] [1] 183 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Chgb,Lrpprc,Slc38a5,Cck,Cdkn1a...
#> ℹ [2026-04-27 16:44:05] 
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
#> ℹ [2026-04-27 16:44:08] Start find dynamic features
#> ℹ [2026-04-27 16:44:09] Data type is raw counts
#> ℹ [2026-04-27 16:44:10] Number of candidate features (union): 2
#> ℹ [2026-04-27 16:44:10] Data type is raw counts
#> ℹ [2026-04-27 16:44:10] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-27 16:44:10] Using 1 core
#> ⠙ [2026-04-27 16:44:10] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-27 16:44:10] Completed 2 tasks in 164ms
#> 
#> ℹ [2026-04-27 16:44:10] Building results
#> ✔ [2026-04-27 16:44:11] Find dynamic features done
#> ℹ [2026-04-27 16:44:11] Start find dynamic features
#> ℹ [2026-04-27 16:44:12] Data type is raw counts
#> ℹ [2026-04-27 16:44:12] Number of candidate features (union): 2
#> ℹ [2026-04-27 16:44:13] Data type is raw counts
#> ℹ [2026-04-27 16:44:13] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-27 16:44:13] Using 1 core
#> ⠙ [2026-04-27 16:44:13] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-27 16:44:13] Completed 2 tasks in 147ms
#> 
#> ℹ [2026-04-27 16:44:13] Building results
#> ✔ [2026-04-27 16:44:13] Find dynamic features done
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
#> ℹ [2026-04-27 16:44:15] Start find dynamic features
#> ℹ [2026-04-27 16:44:15] Data type is raw counts
#> ℹ [2026-04-27 16:44:17] Number of candidate features (union): 236
#> ℹ [2026-04-27 16:44:18] Data type is raw counts
#> ℹ [2026-04-27 16:44:18] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-27 16:44:18] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-04-27 16:44:18] Find dynamic features done
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
#> ℹ [2026-04-27 16:44:18] [1] 172 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Chgb,Gast,Lrpprc,Slc38a5,Cck,Cdkn1a...
#> ℹ [2026-04-27 16:44:19] 
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
#> ℹ [2026-04-27 16:44:22] Start find dynamic features
#> ℹ [2026-04-27 16:44:23] Data type is raw counts
#> ℹ [2026-04-27 16:44:24] Number of candidate features (union): 2
#> ℹ [2026-04-27 16:44:24] Data type is raw counts
#> ℹ [2026-04-27 16:44:24] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-27 16:44:24] Using 1 core
#> ⠙ [2026-04-27 16:44:24] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-27 16:44:24] Completed 2 tasks in 174ms
#> 
#> ℹ [2026-04-27 16:44:24] Building results
#> ✔ [2026-04-27 16:44:25] Find dynamic features done
#> ℹ [2026-04-27 16:44:25] Start find dynamic features
#> ℹ [2026-04-27 16:44:26] Data type is raw counts
#> ℹ [2026-04-27 16:44:26] Number of candidate features (union): 2
#> ℹ [2026-04-27 16:44:27] Data type is raw counts
#> ℹ [2026-04-27 16:44:27] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-27 16:44:27] Using 1 core
#> ⠙ [2026-04-27 16:44:27] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-27 16:44:27] Completed 2 tasks in 152ms
#> 
#> ℹ [2026-04-27 16:44:27] Building results
#> ✔ [2026-04-27 16:44:27] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
