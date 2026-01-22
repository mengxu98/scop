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
#> ℹ [2026-01-22 04:01:55] Start standard scop workflow...
#> ℹ [2026-01-22 04:01:57] Checking a list of <Seurat>...
#> ! [2026-01-22 04:01:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-22 04:01:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:01:59] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 04:02:00] Use the separate HVF from srt_list
#> ℹ [2026-01-22 04:02:00] Number of available HVF: 2000
#> ℹ [2026-01-22 04:02:00] Finished check
#> ℹ [2026-01-22 04:02:00] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-22 04:02:01] Perform pca linear dimension reduction
#> ℹ [2026-01-22 04:02:01] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-22 04:02:02] Reorder clusters...
#> ℹ [2026-01-22 04:02:02] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-22 04:02:02] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-22 04:02:06] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-22 04:02:11] Run scop standard workflow completed
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
#> ℹ [2026-01-22 04:02:13] Start find dynamic features
#> ℹ [2026-01-22 04:02:14] Data type is raw counts
#> ℹ [2026-01-22 04:02:15] Number of candidate features (union): 231
#> ℹ [2026-01-22 04:02:16] Data type is raw counts
#> ℹ [2026-01-22 04:02:16] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-22 04:02:16] Using 1 core
#> ⠙ [2026-01-22 04:02:16] Running for Gcg [1/231] ■                              …
#> ⠹ [2026-01-22 04:02:16] Running for Cpe [40/231] ■■■■■■                        …
#> ⠸ [2026-01-22 04:02:16] Running for Pde10a [129/231] ■■■■■■■■■■■■■■■■■■        …
#> ⠼ [2026-01-22 04:02:16] Running for Sparc [208/231] ■■■■■■■■■■■■■■■■■■■■■■■■■■■…
#> ✔ [2026-01-22 04:02:16] Completed 231 tasks in 8.2s
#> 
#> ℹ [2026-01-22 04:02:16] Building results
#> ℹ [2026-01-22 04:02:24] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-01-22 04:02:24] Using 1 core
#> ⠙ [2026-01-22 04:02:24] Running for Irs4 [60/231] ■■■■■■■■■                    …
#> ⠹ [2026-01-22 04:02:24] Running for Ascl1 [133/231] ■■■■■■■■■■■■■■■■■■         …
#> ⠸ [2026-01-22 04:02:24] Running for H19 [217/231] ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■…
#> ✔ [2026-01-22 04:02:24] Completed 231 tasks in 8.7s
#> 
#> ℹ [2026-01-22 04:02:24] Building results
#> ✔ [2026-01-22 04:02:33] Find dynamic features done

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
#> ℹ [2026-01-22 04:02:33] [1] 180 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Lrpprc,Slc38a5,Cdkn1a,2810417H13Rik,Chga...
#> ℹ [2026-01-22 04:02:35] 
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
#> ℹ [2026-01-22 04:02:38] Start find dynamic features
#> ℹ [2026-01-22 04:02:39] Data type is raw counts
#> ℹ [2026-01-22 04:02:39] Number of candidate features (union): 2
#> ℹ [2026-01-22 04:02:40] Data type is raw counts
#> ℹ [2026-01-22 04:02:40] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-22 04:02:40] Using 1 core
#> ⠙ [2026-01-22 04:02:40] Running for Arxes1 [1/2] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-01-22 04:02:40] Completed 2 tasks in 162ms
#> 
#> ℹ [2026-01-22 04:02:40] Building results
#> ✔ [2026-01-22 04:02:40] Find dynamic features done
#> ℹ [2026-01-22 04:02:40] Start find dynamic features
#> ℹ [2026-01-22 04:02:40] Data type is raw counts
#> ℹ [2026-01-22 04:02:41] Number of candidate features (union): 2
#> ℹ [2026-01-22 04:02:41] Data type is raw counts
#> ℹ [2026-01-22 04:02:41] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-01-22 04:02:41] Using 1 core
#> ⠙ [2026-01-22 04:02:41] Running for Arxes1 [1/2] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-01-22 04:02:41] Completed 2 tasks in 132ms
#> 
#> ℹ [2026-01-22 04:02:41] Building results
#> ✔ [2026-01-22 04:02:42] Find dynamic features done
```
