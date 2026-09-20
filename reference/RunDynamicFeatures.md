# Calculates dynamic features for lineages

Calculates dynamic features for lineages

## Usage

``` r
RunDynamicFeatures(
  object,
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
  pretsa_backend = c("cpp", "r"),
  padjust_method = "fdr",
  cores = 1,
  verbose = TRUE,
  seed = 11,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- lineages:

  Lineage names for which dynamic features should be calculated.

- features:

  Features to use. If `NULL`, n_candidates must be provided.

- suffix:

  Suffix to append to the output layer names for each lineage. Default
  is the lineage names.

- n_candidates:

  A number of candidate features to select when features is `NULL`.

- minfreq:

  Minimum frequency threshold for candidate features. Features with a
  frequency less than minfreq will be excluded.

- family:

  A character or character vector specifying the family of distributions
  to use for the GAM. If family is set to NULL, the appropriate family
  will be automatically determined based on the data. If length(family)
  is 1, the same family will be used for all features. Otherwise, family
  must have the same length as features.

- layer:

  Assay layer to use.

- assay:

  Assay to use. `NULL` uses the default assay.

- libsize:

  A numeric or numeric vector specifying the library size correction
  factors for each cell. If NULL, the library size correction factors
  will be calculated based on the expression matrix. If length(libsize)
  is 1, the same value will be used for all cells. Otherwise, libsize
  must have the same length as the number of cells in srt.

- fit_method:

  The method used for fitting features. Either `"gam"` (generalized
  additive models) or `"pretsa"` (Pattern recognition in Temporal and
  Spatial Analyses).

- knot:

  For `fit_method = "pretsa"`: B-spline knots. `0` or `"auto"`.

- max_knot_allowed:

  For `fit_method = "pretsa"` when `knot = "auto"`: max knots.

- pretsa_backend:

  PreTSA fitting backend. `"cpp"` batches model selection and fitting;
  `"r"` retains the reference implementation.

- padjust_method:

  The method used for p-value adjustment.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Optional integer seed. When supplied, every input receives a
  deterministic independent L'Ecuyer-CMRG random-number stream, making
  results reproducible across worker counts and scheduling order. The
  caller's random number state is restored when the call finishes.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

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
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-20 22:16:56] Start standard processing workflow...
#> ℹ [2026-09-20 22:16:56] Checking a list of <Seurat>...
#> ! [2026-09-20 22:16:56] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-20 22:16:56] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:16:56] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-20 22:16:57] Use the separate HVF from `srt_list`
#> ℹ [2026-09-20 22:16:57] Number of available HVF: 2000
#> ℹ [2026-09-20 22:16:57] Finished check
#> ℹ [2026-09-20 22:16:57] Perform `ScaleData()`
#> ℹ [2026-09-20 22:16:57] Perform pca linear dimension reduction
#> ℹ [2026-09-20 22:16:57] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-20 22:16:57] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-20 22:16:57] Reorder clusters...
#> ℹ [2026-09-20 22:16:58] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-20 22:16:58] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-20 22:17:04] Standard processing workflow completed
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
#> ℹ [2026-09-20 22:17:07] Start find dynamic features
#> ℹ [2026-09-20 22:17:07] Data type is raw counts
#> ℹ [2026-09-20 22:17:08] Number of candidate features (union): 242
#> ℹ [2026-09-20 22:17:08] Data type is raw counts
#> ℹ [2026-09-20 22:17:08] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 22:17:08] Using 1 core
#> ⠙ [2026-09-20 22:17:08] Running for Gcg [1/242]              0% | ETA: 10s
#> ⠹ [2026-09-20 22:17:08] Running for Gm8773 [91/242] ■■■         38% | ETA:  5s
#> ⠸ [2026-09-20 22:17:08] Running for Psat1 [177/242] ■■■■■■■     73% | ETA:  2s
#> ✔ [2026-09-20 22:17:08] Completed 242 tasks in 8.5s
#> 
#> ℹ [2026-09-20 22:17:08] Building results
#> ℹ [2026-09-20 22:17:17] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-09-20 22:17:17] Using 1 core
#> ⠙ [2026-09-20 22:17:17] Running for Krtap17-1 [19/242]              8% | ETA:  …
#> ⠹ [2026-09-20 22:17:17] Running for Plk1 [103/242] ■■■■        43% | ETA:  5s
#> ⠸ [2026-09-20 22:17:17] Running for Nptx2 [184/242] ■■■■■■■     76% | ETA:  2s
#> ✔ [2026-09-20 22:17:17] Completed 242 tasks in 8.5s
#> 
#> ℹ [2026-09-20 22:17:17] Building results
#> ✔ [2026-09-20 22:17:25] Find dynamic features done

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
#> ℹ [2026-09-20 22:17:25] [1] 182 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Lrpprc,Slc38a5,Cck,Cdkn1a,Chga...
#> ℹ [2026-09-20 22:17:26] 
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
#> ℹ [2026-09-20 22:17:28] Start find dynamic features
#> ℹ [2026-09-20 22:17:29] Data type is raw counts
#> ℹ [2026-09-20 22:17:29] Number of candidate features (union): 2
#> ℹ [2026-09-20 22:17:29] Data type is raw counts
#> ℹ [2026-09-20 22:17:29] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 22:17:29] Using 1 core
#> ⠙ [2026-09-20 22:17:29] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 22:17:29] Completed 2 tasks in 137ms
#> 
#> ℹ [2026-09-20 22:17:29] Building results
#> ✔ [2026-09-20 22:17:30] Find dynamic features done
#> ℹ [2026-09-20 22:17:30] Start find dynamic features
#> ℹ [2026-09-20 22:17:30] Data type is raw counts
#> ℹ [2026-09-20 22:17:30] Number of candidate features (union): 2
#> ℹ [2026-09-20 22:17:31] Data type is raw counts
#> ℹ [2026-09-20 22:17:31] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-09-20 22:17:31] Using 1 core
#> ℹ [2026-09-20 22:17:31] Building results
#> ✔ [2026-09-20 22:17:31] Find dynamic features done


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "pretsa"
)
#> ℹ [2026-09-20 22:17:32] Start find dynamic features
#> ℹ [2026-09-20 22:17:32] Data type is raw counts
#> ℹ [2026-09-20 22:17:33] Number of candidate features (union): 242
#> ℹ [2026-09-20 22:17:33] Data type is raw counts
#> ℹ [2026-09-20 22:17:33] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 22:17:33] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-09-20 22:17:34] Find dynamic features done
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
#> ℹ [2026-09-20 22:17:34] [1] 176 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Ppy,Cck...
#> ℹ [2026-09-20 22:17:34] 
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
#> ℹ [2026-09-20 22:17:37] Start find dynamic features
#> ℹ [2026-09-20 22:17:37] Data type is raw counts
#> ℹ [2026-09-20 22:17:37] Number of candidate features (union): 2
#> ℹ [2026-09-20 22:17:38] Data type is raw counts
#> ℹ [2026-09-20 22:17:38] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-09-20 22:17:38] Using 1 core
#> ⠙ [2026-09-20 22:17:38] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 22:17:38] Completed 2 tasks in 139ms
#> 
#> ℹ [2026-09-20 22:17:38] Building results
#> ✔ [2026-09-20 22:17:38] Find dynamic features done
#> ℹ [2026-09-20 22:17:38] Start find dynamic features
#> ℹ [2026-09-20 22:17:38] Data type is raw counts
#> ℹ [2026-09-20 22:17:39] Number of candidate features (union): 2
#> ℹ [2026-09-20 22:17:39] Data type is raw counts
#> ℹ [2026-09-20 22:17:39] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-09-20 22:17:39] Using 1 core
#> ⠙ [2026-09-20 22:17:39] Running for Arxes1 [1/2] ■■■■■       50% | ETA:  0s
#> ✔ [2026-09-20 22:17:39] Completed 2 tasks in 134ms
#> 
#> ℹ [2026-09-20 22:17:39] Building results
#> ✔ [2026-09-20 22:17:39] Find dynamic features done
```
