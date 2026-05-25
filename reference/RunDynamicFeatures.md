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
#> ℹ [2026-05-25 05:13:54] Start standard processing workflow...
#> ℹ [2026-05-25 05:13:54] Checking a list of <Seurat>...
#> ! [2026-05-25 05:13:55] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:13:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:13:57] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:13:57] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 05:13:57] Number of available HVF: 2000
#> ℹ [2026-05-25 05:13:57] Finished check
#> ℹ [2026-05-25 05:13:57] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 05:13:58] Perform pca linear dimension reduction
#> ℹ [2026-05-25 05:13:58] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-25 05:13:59] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 05:13:59] Reorder clusters...
#> ℹ [2026-05-25 05:13:59] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 05:13:59] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 05:13:59] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-25 05:14:05] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-25 05:14:11] Standard processing workflow completed
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
#> ℹ [2026-05-25 05:14:13] Start find dynamic features
#> ℹ [2026-05-25 05:14:14] Data type is raw counts
#> ℹ [2026-05-25 05:14:16] Number of candidate features (union): 244
#> ℹ [2026-05-25 05:14:17] Data type is raw counts
#> ℹ [2026-05-25 05:14:17] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-25 05:14:17] Using 1 core
#> ⠙ [2026-05-25 05:14:17] Running for Gcg [1/244]              0% | ETA: 24s
#> ⠹ [2026-05-25 05:14:17] Running for Sdf2l1 [41/244] ■           17% | ETA: 16s
#> ⠸ [2026-05-25 05:14:17] Running for Aard [82/244] ■■■         34% | ETA: 12s
#> ⠼ [2026-05-25 05:14:17] Running for Id2 [126/244] ■■■■■       52% | ETA:  8s
#> ⠴ [2026-05-25 05:14:17] Running for Pou3f1 [169/244] ■■■■■■      69% | ETA:  5s
#> ⠦ [2026-05-25 05:14:17] Running for Xist [212/244] ■■■■■■■■    87% | ETA:  2s
#> ✔ [2026-05-25 05:14:17] Completed 244 tasks in 17.2s
#> 
#> ℹ [2026-05-25 05:14:17] Building results
#> ℹ [2026-05-25 05:14:34] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-25 05:14:34] Using 1 core
#> ⠙ [2026-05-25 05:14:34] Running for Gip [18/244]              7% | ETA:  9s
#> ⠹ [2026-05-25 05:14:34] Running for Hes1 [91/244] ■■■         37% | ETA:  6s
#> ⠸ [2026-05-25 05:14:34] Running for Foxd3 [161/244] ■■■■■■      66% | ETA:  3s
#> ⠼ [2026-05-25 05:14:34] Running for 1500009L16Rik [240/244] ■■■■■■■■■   98% | E…
#> ✔ [2026-05-25 05:14:34] Completed 244 tasks in 10s
#> 
#> ℹ [2026-05-25 05:14:34] Building results
#> ✔ [2026-05-25 05:14:44] Find dynamic features done

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
#> ℹ [2026-05-25 05:14:44] [1] 184 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Ins2,Nnat,Iapp,Lrpprc,Chgb,Pyy,Slc38a5,2810417H13Rik...
#> ℹ [2026-05-25 05:14:45] 
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
#> ℹ [2026-05-25 05:14:48] Start find dynamic features
#> ℹ [2026-05-25 05:14:49] Data type is raw counts
#> ℹ [2026-05-25 05:14:50] Number of candidate features (union): 2
#> ℹ [2026-05-25 05:14:50] Data type is raw counts
#> ℹ [2026-05-25 05:14:50] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-25 05:14:50] Using 1 core
#> ⠙ [2026-05-25 05:14:50] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-25 05:14:50] Completed 2 tasks in 204ms
#> 
#> ℹ [2026-05-25 05:14:50] Building results
#> ✔ [2026-05-25 05:14:51] Find dynamic features done
#> ℹ [2026-05-25 05:14:51] Start find dynamic features
#> ℹ [2026-05-25 05:14:52] Data type is raw counts
#> ℹ [2026-05-25 05:14:53] Number of candidate features (union): 2
#> ℹ [2026-05-25 05:14:53] Data type is raw counts
#> ℹ [2026-05-25 05:14:53] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-25 05:14:53] Using 1 core
#> ⠙ [2026-05-25 05:14:53] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-25 05:14:53] Completed 2 tasks in 152ms
#> 
#> ℹ [2026-05-25 05:14:53] Building results
#> ✔ [2026-05-25 05:14:54] Find dynamic features done
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
#> ℹ [2026-05-25 05:14:55] Start find dynamic features
#> ℹ [2026-05-25 05:14:56] Data type is raw counts
#> ℹ [2026-05-25 05:14:58] Number of candidate features (union): 244
#> ℹ [2026-05-25 05:14:58] Data type is raw counts
#> ℹ [2026-05-25 05:14:58] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-25 05:14:58] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-05-25 05:14:59] Find dynamic features done
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
#> ℹ [2026-05-25 05:14:59] [1] 167 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ins1,Ins2,Nnat,Iapp,Lrpprc,Chgb,Pyy,Slc38a5,2810417H13Rik...
#> ℹ [2026-05-25 05:15:00] 
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
#> ℹ [2026-05-25 05:15:03] Start find dynamic features
#> ℹ [2026-05-25 05:15:06] Data type is raw counts
#> ℹ [2026-05-25 05:15:06] Number of candidate features (union): 2
#> ℹ [2026-05-25 05:15:07] Data type is raw counts
#> ℹ [2026-05-25 05:15:07] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-25 05:15:07] Using 1 core
#> ⠙ [2026-05-25 05:15:07] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-25 05:15:07] Completed 2 tasks in 248ms
#> 
#> ℹ [2026-05-25 05:15:07] Building results
#> ✔ [2026-05-25 05:15:07] Find dynamic features done
#> ℹ [2026-05-25 05:15:07] Start find dynamic features
#> ℹ [2026-05-25 05:15:09] Data type is raw counts
#> ℹ [2026-05-25 05:15:09] Number of candidate features (union): 2
#> ℹ [2026-05-25 05:15:10] Data type is raw counts
#> ℹ [2026-05-25 05:15:10] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-25 05:15:10] Using 1 core
#> ⠙ [2026-05-25 05:15:10] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-05-25 05:15:10] Completed 2 tasks in 172ms
#> 
#> ℹ [2026-05-25 05:15:10] Building results
#> ✔ [2026-05-25 05:15:10] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
