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
#> ℹ [2026-06-24 04:13:03] Start standard processing workflow...
#> ℹ [2026-06-24 04:13:04] Checking a list of <Seurat>...
#> ! [2026-06-24 04:13:04] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-06-24 04:13:04] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:13:04] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-06-24 04:13:04] Use the separate HVF from `srt_list`
#> ℹ [2026-06-24 04:13:04] Number of available HVF: 2000
#> ℹ [2026-06-24 04:13:04] Finished check
#> ℹ [2026-06-24 04:13:04] Perform `ScaleData()`
#> ℹ [2026-06-24 04:13:04] Perform pca linear dimension reduction
#> ℹ [2026-06-24 04:13:05] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-06-24 04:13:05] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-06-24 04:13:05] Reorder clusters...
#> ℹ [2026-06-24 04:13:07] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-06-24 04:13:07] Perform umap nonlinear dimension reduction
#> ✔ [2026-06-24 04:13:15] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 11 rows containing missing values or values outside the scale range
#> (`geom_path()`).


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "gam"
)
#> ℹ [2026-06-24 04:13:17] Start find dynamic features
#> ℹ [2026-06-24 04:13:18] Data type is raw counts
#> ℹ [2026-06-24 04:13:19] Number of candidate features (union): 242
#> ℹ [2026-06-24 04:13:20] Data type is raw counts
#> ℹ [2026-06-24 04:13:20] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-06-24 04:13:20] Using 1 core
#> ⠙ [2026-06-24 04:13:20] Running for Gcg [1/242]              0% | ETA:  9s
#> ⠹ [2026-06-24 04:13:20] Running for Rbp4 [5/242]              2% | ETA:  8s
#> ⠸ [2026-06-24 04:13:20] Running for Gip [102/242] ■■■■        42% | ETA:  4s
#> ⠼ [2026-06-24 04:13:20] Running for Mid1ip1 [196/242] ■■■■■■■■    81% | ETA:  1s
#> ✔ [2026-06-24 04:13:20] Completed 242 tasks in 7.7s
#> 
#> ℹ [2026-06-24 04:13:20] Building results
#> ℹ [2026-06-24 04:13:28] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-06-24 04:13:28] Using 1 core
#> ⠙ [2026-06-24 04:13:28] Running for Top2a [43/242] ■           18% | ETA:  7s
#> ⠹ [2026-06-24 04:13:28] Running for Slc30a8 [134/242] ■■■■■       55% | ETA:  4s
#> ⠸ [2026-06-24 04:13:28] Running for Sds [226/242] ■■■■■■■■■   93% | ETA:  1s
#> ✔ [2026-06-24 04:13:28] Completed 242 tasks in 8s
#> 
#> ℹ [2026-06-24 04:13:28] Building results
#> ✔ [2026-06-24 04:13:36] Find dynamic features done

names(
  pancreas_sub@tools$DynamicFeatures_Lineage1
)
#> [1] "DynamicFeatures" "raw_matrix"      "fitted_matrix"   "upr_matrix"     
#> [5] "lwr_matrix"      "libsize"         "lineages"        "family"         
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells       r.sq  dev.expl peaktime valleytime pvalue
#> Gcg       Gcg        174 0.58125774 0.7804150 22.76792  13.520617      0
#> Ghrl     Ghrl        162 0.36100387 0.6970627 19.26746  22.767919      0
#> Iapp     Iapp        265 0.28289448 0.7644470 22.76792   0.100214      0
#> Pyy       Pyy        404 0.43607110 0.7930835 20.83092   0.100214      0
#> Rbp4     Rbp4        363 0.45900832 0.7533428 19.78387  11.043897      0
#> Gast     Gast         87 0.07789092 0.7359666 20.60505   0.100214      0
#>      padjust
#> Gcg        0
#> Ghrl       0
#> Iapp       0
#> Pyy        0
#> Rbp4       0
#> Gast       0
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-06-24 04:13:36] [1] 182 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Lrpprc,Slc38a5,Cck,Cdkn1a,Chga...
#> ℹ [2026-06-24 04:13:37] 
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
#> ℹ [2026-06-24 04:13:39] Start find dynamic features
#> ℹ [2026-06-24 04:13:41] Data type is raw counts
#> ℹ [2026-06-24 04:13:41] Number of candidate features (union): 2
#> ℹ [2026-06-24 04:13:42] Data type is raw counts
#> ℹ [2026-06-24 04:13:42] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-06-24 04:13:42] Using 1 core
#> ⠙ [2026-06-24 04:13:42] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-06-24 04:13:42] Completed 2 tasks in 137ms
#> 
#> ℹ [2026-06-24 04:13:42] Building results
#> ✔ [2026-06-24 04:13:42] Find dynamic features done
#> ℹ [2026-06-24 04:13:42] Start find dynamic features
#> ℹ [2026-06-24 04:13:43] Data type is raw counts
#> ℹ [2026-06-24 04:13:44] Number of candidate features (union): 2
#> ℹ [2026-06-24 04:13:44] Data type is raw counts
#> ℹ [2026-06-24 04:13:44] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-06-24 04:13:44] Using 1 core
#> ⠙ [2026-06-24 04:13:44] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-06-24 04:13:44] Completed 2 tasks in 139ms
#> 
#> ℹ [2026-06-24 04:13:44] Building results
#> ✔ [2026-06-24 04:13:44] Find dynamic features done
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
#> ℹ [2026-06-24 04:13:46] Start find dynamic features
#> ℹ [2026-06-24 04:13:46] Data type is raw counts
#> ℹ [2026-06-24 04:13:48] Number of candidate features (union): 242
#> ℹ [2026-06-24 04:13:48] Data type is raw counts
#> ℹ [2026-06-24 04:13:48] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-06-24 04:13:49] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-06-24 04:13:49] Find dynamic features done
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime valleytime        pvalue
#> Gcg       Gcg        174 0.5614256 0.5614256 22.76792  13.421731 5.941533e-122
#> Ghrl     Ghrl        162 0.1488719 0.1488719 19.09239   5.283551  9.293807e-24
#> Iapp     Iapp        265 0.6606537 0.6606537 22.76792   0.100214 5.131126e-160
#> Pyy       Pyy        404 0.7076209 0.7076209 22.76792   8.000384 3.969429e-182
#> Rbp4     Rbp4        363 0.6462734 0.6462734 22.76792   7.756116 7.410789e-154
#> Gast     Gast         87 0.2620377 0.2620377 22.76792   8.000384  7.917567e-45
#>            padjust
#> Gcg  2.765098e-121
#> Ghrl  1.330829e-23
#> Iapp 5.173885e-159
#> Pyy  6.003762e-181
#> Rbp4 6.405039e-153
#> Gast  1.321415e-44
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-06-24 04:13:49] [1] 176 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
#> ℹ [2026-06-24 04:13:50] 
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
#> ℹ [2026-06-24 04:13:52] Start find dynamic features
#> ℹ [2026-06-24 04:13:54] Data type is raw counts
#> ℹ [2026-06-24 04:13:54] Number of candidate features (union): 2
#> ℹ [2026-06-24 04:13:55] Data type is raw counts
#> ℹ [2026-06-24 04:13:55] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-06-24 04:13:55] Using 1 core
#> ⠙ [2026-06-24 04:13:55] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-06-24 04:13:55] Completed 2 tasks in 150ms
#> 
#> ℹ [2026-06-24 04:13:55] Building results
#> ✔ [2026-06-24 04:13:55] Find dynamic features done
#> ℹ [2026-06-24 04:13:55] Start find dynamic features
#> ℹ [2026-06-24 04:13:56] Data type is raw counts
#> ℹ [2026-06-24 04:13:57] Number of candidate features (union): 2
#> ℹ [2026-06-24 04:13:57] Data type is raw counts
#> ℹ [2026-06-24 04:13:57] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-06-24 04:13:57] Using 1 core
#> ⠙ [2026-06-24 04:13:57] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-06-24 04:13:57] Completed 2 tasks in 136ms
#> 
#> ℹ [2026-06-24 04:13:57] Building results
#> ✔ [2026-06-24 04:13:57] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
