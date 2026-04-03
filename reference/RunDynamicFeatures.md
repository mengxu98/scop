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
#> ℹ [2026-04-03 09:59:39] Start standard processing workflow...
#> ℹ [2026-04-03 09:59:40] Checking a list of <Seurat>...
#> ! [2026-04-03 09:59:40] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-03 09:59:40] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 09:59:42] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-03 09:59:43] Use the separate HVF from `srt_list`
#> ℹ [2026-04-03 09:59:43] Number of available HVF: 2000
#> ℹ [2026-04-03 09:59:43] Finished check
#> ℹ [2026-04-03 09:59:43] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-03 09:59:44] Perform pca linear dimension reduction
#> ℹ [2026-04-03 09:59:44] Use stored estimated dimensions 1:12 for Standardpca
#> ℹ [2026-04-03 09:59:45] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-03 09:59:45] Reorder clusters...
#> ℹ [2026-04-03 09:59:45] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-03 09:59:45] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-03 09:59:45] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ℹ [2026-04-03 09:59:50] Perform umap nonlinear dimension reduction using Standardpca (1:12)
#> ✔ [2026-04-03 09:59:55] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Warning: Removed 15 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 15 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_path()`).


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "gam"
)
#> ℹ [2026-04-03 09:59:56] Start find dynamic features
#> ℹ [2026-04-03 09:59:57] Data type is raw counts
#> ℹ [2026-04-03 09:59:59] Number of candidate features (union): 223
#> ℹ [2026-04-03 10:00:00] Data type is raw counts
#> ℹ [2026-04-03 10:00:00] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-03 10:00:00] Using 1 core
#> ⠙ [2026-04-03 10:00:00] Running for Ghrl [1/223]              0% | ETA: 10s
#> ⠹ [2026-04-03 10:00:00] Running for 8430408G22Rik [29/223] ■           13% | ET…
#> ⠸ [2026-04-03 10:00:00] Running for Ptn [108/223] ■■■■        48% | ETA:  4s
#> ⠼ [2026-04-03 10:00:00] Running for Kctd14 [188/223] ■■■■■■■■    84% | ETA:  1s
#> ✔ [2026-04-03 10:00:00] Completed 223 tasks in 8.5s
#> 
#> ℹ [2026-04-03 10:00:00] Building results
#> ℹ [2026-04-03 10:00:08] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-03 10:00:08] Using 1 core
#> ⠙ [2026-04-03 10:00:08] Running for Tagln [42/223] ■           19% | ETA:  6s
#> ⠹ [2026-04-03 10:00:08] Running for Rgs4 [113/223] ■■■■■       51% | ETA:  4s
#> ⠸ [2026-04-03 10:00:08] Running for Notch2 [184/223] ■■■■■■■■    83% | ETA:  2s
#> ✔ [2026-04-03 10:00:08] Completed 223 tasks in 8.8s
#> 
#> ℹ [2026-04-03 10:00:08] Building results
#> ✔ [2026-04-03 10:00:17] Find dynamic features done

names(
  pancreas_sub@tools$DynamicFeatures_Lineage1
)
#> [1] "DynamicFeatures" "raw_matrix"      "fitted_matrix"   "upr_matrix"     
#> [5] "lwr_matrix"      "libsize"         "lineages"        "family"         
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime  valleytime pvalue
#> Ghrl     Ghrl        159 0.2808128 0.6151289 25.25569  0.09701269      0
#> Gcg       Gcg        176 0.6329316 0.8037828 28.09067 15.22083623      0
#> Iapp     Iapp        270 0.3098506 0.7603657 28.09067  0.09701269      0
#> Pyy       Pyy        411 0.4122953 0.7180595 26.78263  0.09701269      0
#> Rbp4     Rbp4        388 0.4207574 0.7013749 26.23850 13.92893485      0
#> Chgb     Chgb        268 0.4884483 0.7313913 21.51067  0.09701269      0
#>      padjust
#> Ghrl       0
#> Gcg        0
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
#> ℹ [2026-04-03 10:00:17] [1] 173 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Ghrl,Gcg,Iapp,Pyy,Rbp4,Chgb,Lrpprc,Slc38a5,2810417H13Rik,Cdc20...
#> ℹ [2026-04-03 10:00:18] 
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
#> ℹ [2026-04-03 10:00:21] Start find dynamic features
#> ℹ [2026-04-03 10:00:22] Data type is raw counts
#> ℹ [2026-04-03 10:00:23] Number of candidate features (union): 2
#> ℹ [2026-04-03 10:00:23] Data type is raw counts
#> ℹ [2026-04-03 10:00:23] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-03 10:00:23] Using 1 core
#> ⠙ [2026-04-03 10:00:23] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-03 10:00:23] Completed 2 tasks in 142ms
#> 
#> ℹ [2026-04-03 10:00:23] Building results
#> ✔ [2026-04-03 10:00:23] Find dynamic features done
#> ℹ [2026-04-03 10:00:23] Start find dynamic features
#> ℹ [2026-04-03 10:00:25] Data type is raw counts
#> ℹ [2026-04-03 10:00:25] Number of candidate features (union): 2
#> ℹ [2026-04-03 10:00:26] Data type is raw counts
#> ℹ [2026-04-03 10:00:26] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-03 10:00:26] Using 1 core
#> ⠙ [2026-04-03 10:00:26] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-03 10:00:26] Completed 2 tasks in 137ms
#> 
#> ℹ [2026-04-03 10:00:26] Building results
#> ✔ [2026-04-03 10:00:26] Find dynamic features done
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
#> ℹ [2026-04-03 10:00:27] Start find dynamic features
#> ℹ [2026-04-03 10:00:28] Data type is raw counts
#> ℹ [2026-04-03 10:00:30] Number of candidate features (union): 223
#> ℹ [2026-04-03 10:00:30] Data type is raw counts
#> ℹ [2026-04-03 10:00:30] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-03 10:00:30] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-04-03 10:00:30] Find dynamic features done
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime  valleytime        pvalue
#> Ghrl     Ghrl        159 0.1196995 0.1196995 28.09067  7.44754441  1.125033e-19
#> Gcg       Gcg        176 0.5306269 0.5306269 28.09067  8.46663469 3.947382e-117
#> Iapp     Iapp        270 0.6378814 0.6378814 28.09067  0.09701269 2.000605e-157
#> Pyy       Pyy        411 0.6748659 0.6748659 28.09067 12.42896282 3.657064e-174
#> Rbp4     Rbp4        388 0.6412010 0.6412010 28.09067 11.01057150 7.420358e-159
#> Chgb     Chgb        268 0.4061405 0.4061405 24.64629  5.90494191  1.299702e-80
#>            padjust
#> Ghrl  1.475779e-19
#> Gcg  1.354256e-116
#> Iapp 1.538397e-156
#> Pyy  5.825180e-173
#> Rbp4 5.909785e-158
#> Chgb  2.611112e-80
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-04-03 10:00:31] [1] 168 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Chgb,Gast,Lrpprc,Ppy,Slc38a5,2810417H13Rik...
#> ℹ [2026-04-03 10:00:32] 
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
#> ℹ [2026-04-03 10:00:34] Start find dynamic features
#> ℹ [2026-04-03 10:00:36] Data type is raw counts
#> ℹ [2026-04-03 10:00:36] Number of candidate features (union): 2
#> ℹ [2026-04-03 10:00:37] Data type is raw counts
#> ℹ [2026-04-03 10:00:37] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-04-03 10:00:37] Using 1 core
#> ⠙ [2026-04-03 10:00:37] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-03 10:00:37] Completed 2 tasks in 144ms
#> 
#> ℹ [2026-04-03 10:00:37] Building results
#> ✔ [2026-04-03 10:00:37] Find dynamic features done
#> ℹ [2026-04-03 10:00:37] Start find dynamic features
#> ℹ [2026-04-03 10:00:40] Data type is raw counts
#> ℹ [2026-04-03 10:00:40] Number of candidate features (union): 2
#> ℹ [2026-04-03 10:00:41] Data type is raw counts
#> ℹ [2026-04-03 10:00:41] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-04-03 10:00:41] Using 1 core
#> ⠙ [2026-04-03 10:00:41] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-04-03 10:00:41] Completed 2 tasks in 140ms
#> 
#> ℹ [2026-04-03 10:00:41] Building results
#> ✔ [2026-04-03 10:00:41] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
