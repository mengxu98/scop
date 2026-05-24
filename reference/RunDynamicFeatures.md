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
#> ℹ [2026-05-24 16:08:17] Start standard processing workflow...
#> ℹ [2026-05-24 16:08:18] Checking a list of <Seurat>...
#> ! [2026-05-24 16:08:18] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-24 16:08:18] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-24 16:08:20] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-24 16:08:20] Use the separate HVF from `srt_list`
#> ℹ [2026-05-24 16:08:21] Number of available HVF: 2000
#> ℹ [2026-05-24 16:08:21] Finished check
#> ℹ [2026-05-24 16:08:21] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-24 16:08:21] Perform pca linear dimension reduction
#> ℹ [2026-05-24 16:08:21] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-24 16:08:22] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-24 16:08:22] Reorder clusters...
#> ℹ [2026-05-24 16:08:22] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-24 16:08:22] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-24 16:08:22] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-24 16:08:27] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-24 16:08:33] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "gam"
)
#> ℹ [2026-05-24 16:08:34] Start find dynamic features
#> ℹ [2026-05-24 16:08:36] Data type is raw counts
#> ℹ [2026-05-24 16:08:38] Number of candidate features (union): 244
#> ℹ [2026-05-24 16:08:38] Data type is raw counts
#> ℹ [2026-05-24 16:08:38] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-24 16:08:38] Using 1 core
#> ⠙ [2026-05-24 16:08:38] Running for Gcg [1/244]              0% | ETA: 23s
#> ⠹ [2026-05-24 16:08:38] Running for Krtap17-1 [21/244]              9% | ETA: 1…
#> ⠸ [2026-05-24 16:08:38] Running for Kirrel2 [64/244] ■■          26% | ETA: 13s
#> ⠼ [2026-05-24 16:08:38] Running for Mki67 [107/244] ■■■■        44% | ETA: 10s
#> ⠴ [2026-05-24 16:08:38] Running for H19 [150/244] ■■■■■■      61% | ETA:  7s
#> ⠦ [2026-05-24 16:08:38] Running for Notch2 [191/244] ■■■■■■■     78% | ETA:  4s
#> ⠧ [2026-05-24 16:08:38] Running for Ppp1r14c [233/244] ■■■■■■■■■   95% | ETA:  …
#> ✔ [2026-05-24 16:08:38] Completed 244 tasks in 17.5s
#> 
#> ℹ [2026-05-24 16:08:38] Building results
#> ℹ [2026-05-24 16:08:56] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-24 16:08:56] Using 1 core
#> ⠙ [2026-05-24 16:08:56] Running for Cyr61 [50/244] ■■          20% | ETA:  9s
#> ⠹ [2026-05-24 16:08:56] Running for Asb4 [123/244] ■■■■■       50% | ETA:  5s
#> ⠸ [2026-05-24 16:08:56] Running for Notch2 [191/244] ■■■■■■■     78% | ETA:  2s
#> ✔ [2026-05-24 16:08:56] Completed 244 tasks in 10.4s
#> 
#> ℹ [2026-05-24 16:08:56] Building results
#> ✔ [2026-05-24 16:09:06] Find dynamic features done

names(
  pancreas_sub@tools$DynamicFeatures_Lineage1
)
#> [1] "DynamicFeatures" "raw_matrix"      "fitted_matrix"   "upr_matrix"     
#> [5] "lwr_matrix"      "libsize"         "lineages"        "family"         
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime valleytime pvalue padjust
#> Gcg       Gcg        218 0.4096024 0.5628034 24.45175 0.06479825      0       0
#> Ghrl     Ghrl        211 0.2263349 0.5882084 19.32336 7.41692614      0       0
#> Ins1     Ins1        245 0.0000000 0.6766883 24.45175 2.79616696      0       0
#> Ins2     Ins2        151 0.3676227 0.7611310 23.52609 0.06479825      0       0
#> Nnat     Nnat        339 0.4519823 0.6394348 22.88144 0.06479825      0       0
#> Iapp     Iapp        460 0.4978671 0.7816684 24.45175 0.06479825      0       0
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-05-24 16:09:06] [1] 184 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Ins2,Nnat,Iapp,Lrpprc,Chgb,Pyy,Slc38a5,2810417H13Rik...
#> ℹ [2026-05-24 16:09:08] 
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
#> ℹ [2026-05-24 16:09:10] Start find dynamic features
#> ℹ [2026-05-24 16:09:12] Data type is raw counts
#> ℹ [2026-05-24 16:09:12] Number of candidate features (union): 2
#> ℹ [2026-05-24 16:09:13] Data type is raw counts
#> ℹ [2026-05-24 16:09:13] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-24 16:09:13] Using 1 core
#> ⠙ [2026-05-24 16:09:13] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-24 16:09:13] Completed 2 tasks in 211ms
#> 
#> ℹ [2026-05-24 16:09:13] Building results
#> ✔ [2026-05-24 16:09:13] Find dynamic features done
#> ℹ [2026-05-24 16:09:13] Start find dynamic features
#> ℹ [2026-05-24 16:09:15] Data type is raw counts
#> ℹ [2026-05-24 16:09:15] Number of candidate features (union): 2
#> ℹ [2026-05-24 16:09:16] Data type is raw counts
#> ℹ [2026-05-24 16:09:16] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-24 16:09:16] Using 1 core
#> ℹ [2026-05-24 16:09:16] Building results
#> ✔ [2026-05-24 16:09:16] Find dynamic features done
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
#> ℹ [2026-05-24 16:09:17] Start find dynamic features
#> ℹ [2026-05-24 16:09:19] Data type is raw counts
#> ℹ [2026-05-24 16:09:21] Number of candidate features (union): 244
#> ℹ [2026-05-24 16:09:22] Data type is raw counts
#> ℹ [2026-05-24 16:09:22] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-24 16:09:22] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-05-24 16:09:22] Find dynamic features done
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells       r.sq   dev.expl peaktime  valleytime
#> Gcg       Gcg        218 0.28707130 0.28707130 24.45175 14.87003138
#> Ghrl     Ghrl        211 0.06680993 0.06680993 18.43544  5.54019097
#> Ins1     Ins1        245 0.27413470 0.27413470 24.45175 14.20370086
#> Ins2     Ins2        151 0.34759326 0.34759326 24.45175 14.09979230
#> Nnat     Nnat        339 0.35599396 0.35599396 24.45175  0.06479825
#> Iapp     Iapp        460 0.73111421 0.73111421 24.45175 13.23265992
#>             pvalue       padjust
#> Gcg   2.591714e-71  4.649840e-71
#> Ghrl  1.448465e-14  1.691031e-14
#> Ins1  1.640201e-67  2.900066e-67
#> Ins2  4.501568e-90  9.153188e-90
#> Nnat  8.162278e-93  1.716893e-92
#> Iapp 9.083739e-278 1.847027e-276
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-05-24 16:09:22] [1] 167 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ins1,Ins2,Nnat,Iapp,Lrpprc,Chgb,Pyy,Slc38a5,2810417H13Rik...
#> ℹ [2026-05-24 16:09:23] 
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
#> ℹ [2026-05-24 16:09:26] Start find dynamic features
#> ℹ [2026-05-24 16:09:28] Data type is raw counts
#> ℹ [2026-05-24 16:09:28] Number of candidate features (union): 2
#> ℹ [2026-05-24 16:09:29] Data type is raw counts
#> ℹ [2026-05-24 16:09:29] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-24 16:09:29] Using 1 core
#> ⠙ [2026-05-24 16:09:29] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-24 16:09:29] Completed 2 tasks in 209ms
#> 
#> ℹ [2026-05-24 16:09:29] Building results
#> ✔ [2026-05-24 16:09:29] Find dynamic features done
#> ℹ [2026-05-24 16:09:29] Start find dynamic features
#> ℹ [2026-05-24 16:09:30] Data type is raw counts
#> ℹ [2026-05-24 16:09:31] Number of candidate features (union): 2
#> ℹ [2026-05-24 16:09:31] Data type is raw counts
#> ℹ [2026-05-24 16:09:31] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-24 16:09:31] Using 1 core
#> ⠙ [2026-05-24 16:09:31] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-24 16:09:31] Completed 2 tasks in 158ms
#> 
#> ℹ [2026-05-24 16:09:31] Building results
#> ✔ [2026-05-24 16:09:32] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
