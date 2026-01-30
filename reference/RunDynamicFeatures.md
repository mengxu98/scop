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
  to use for the generalized additive models (GAMs). If family is set to
  NULL, the appropriate family will be automatically determined based on
  the data. If length(family) is 1, the same family will be used for all
  features. Otherwise, family must have the same length as features.

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

## See also

[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md),
[DynamicPlot](https://mengxu98.github.io/scop/reference/DynamicPlot.md),
[RunDynamicEnrichment](https://mengxu98.github.io/scop/reference/RunDynamicEnrichment.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-30 17:13:48] Start standard scop workflow...
#> ℹ [2026-01-30 17:13:48] Checking a list of <Seurat>...
#> ! [2026-01-30 17:13:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-30 17:13:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:13:51] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-30 17:13:51] Use the separate HVF from srt_list
#> ℹ [2026-01-30 17:13:52] Number of available HVF: 2000
#> ℹ [2026-01-30 17:13:52] Finished check
#> ℹ [2026-01-30 17:13:52] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-30 17:13:52] Perform pca linear dimension reduction
#> ℹ [2026-01-30 17:13:53] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-30 17:13:53] Reorder clusters...
#> ℹ [2026-01-30 17:13:54] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-30 17:13:54] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-30 17:13:59] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-30 17:14:04] Run scop standard workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200
)
#> ℹ [2026-01-30 17:14:06] Start find dynamic features
#> ℹ [2026-01-30 17:14:06] Data type is raw counts
#> ℹ [2026-01-30 17:14:08] Number of candidate features (union): 231
#> ℹ [2026-01-30 17:14:09] Data type is raw counts
#> ℹ [2026-01-30 17:14:09] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-30 17:14:09] Using 1 core
#> ⠙ [2026-01-30 17:14:09] Running for Gcg [1/231] ■                              …
#> ⠹ [2026-01-30 17:14:09] Running for Nnat [41/231] ■■■■■■                       …
#> ⠸ [2026-01-30 17:14:09] Running for Arx [127/231] ■■■■■■■■■■■■■■■■■            …
#> ⠼ [2026-01-30 17:14:09] Running for Dut [198/231] ■■■■■■■■■■■■■■■■■■■■■■■■■■■  …
#> ✔ [2026-01-30 17:14:09] Completed 231 tasks in 8.7s
#> 
#> ℹ [2026-01-30 17:14:09] Building results
#> ℹ [2026-01-30 17:14:17] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-01-30 17:14:17] Using 1 core
#> ⠙ [2026-01-30 17:14:17] Running for Top2a [43/231] ■■■■■■■                     …
#> ⠹ [2026-01-30 17:14:17] Running for Rgs4 [113/231] ■■■■■■■■■■■■■■■■            …
#> ⠸ [2026-01-30 17:14:17] Running for Hist1h1c [188/231] ■■■■■■■■■■■■■■■■■■■■■■■■…
#> ✔ [2026-01-30 17:14:17] Completed 231 tasks in 9.4s
#> 
#> ℹ [2026-01-30 17:14:17] Building results
#> ✔ [2026-01-30 17:14:27] Find dynamic features done

names(
  pancreas_sub@tools$DynamicFeatures_Lineage1
)
#> [1] "DynamicFeatures" "raw_matrix"      "fitted_matrix"   "upr_matrix"     
#> [5] "lwr_matrix"      "libsize"         "lineages"        "family"         
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime valleytime pvalue padjust
#> Gcg       Gcg        182 0.6401140 0.7670504 21.33789  0.1164741      0       0
#> Ghrl     Ghrl        167 0.2985260 0.6484728 18.01054 12.8984585      0       0
#> Iapp     Iapp        279 0.2241776 0.7503415 21.33789  0.1164741      0       0
#> Pyy       Pyy        434 0.4169149 0.7697879 19.21200  9.8743469      0       0
#> Rbp4     Rbp4        396 0.4403772 0.7285632 18.83423 10.1022997      0       0
#> Gast     Gast         92 0.0659722 0.7077578 19.84490  0.1164741      0       0
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 6,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-01-30 17:14:27] [1] 180 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Lrpprc,Slc38a5,Cdkn1a,2810417H13Rik,Chga...
#> ℹ [2026-01-30 17:14:28] 
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
#> ℹ [2026-01-30 17:14:32] Start find dynamic features
#> ℹ [2026-01-30 17:14:33] Data type is raw counts
#> ℹ [2026-01-30 17:14:33] Number of candidate features (union): 2
#> ℹ [2026-01-30 17:14:34] Data type is raw counts
#> ℹ [2026-01-30 17:14:34] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-30 17:14:34] Using 1 core
#> ⠙ [2026-01-30 17:14:34] Running for Arxes1 [1/2] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-01-30 17:14:34] Completed 2 tasks in 170ms
#> 
#> ℹ [2026-01-30 17:14:34] Building results
#> ✔ [2026-01-30 17:14:34] Find dynamic features done
#> ℹ [2026-01-30 17:14:34] Start find dynamic features
#> ℹ [2026-01-30 17:14:35] Data type is raw counts
#> ℹ [2026-01-30 17:14:35] Number of candidate features (union): 2
#> ℹ [2026-01-30 17:14:36] Data type is raw counts
#> ℹ [2026-01-30 17:14:36] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-01-30 17:14:36] Using 1 core
#> ℹ [2026-01-30 17:14:36] Building results
#> ✔ [2026-01-30 17:14:36] Find dynamic features done
```
