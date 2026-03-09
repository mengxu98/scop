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
#> ℹ [2026-03-09 08:30:24] Start standard scop workflow...
#> ℹ [2026-03-09 08:30:24] Checking a list of <Seurat>...
#> ! [2026-03-09 08:30:25] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-03-09 08:30:25] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:30:27] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-03-09 08:30:27] Use the separate HVF from `srt_list`
#> ℹ [2026-03-09 08:30:27] Number of available HVF: 2000
#> ℹ [2026-03-09 08:30:28] Finished check
#> ℹ [2026-03-09 08:30:28] Perform `Seurat::ScaleData()`
#> ℹ [2026-03-09 08:30:28] Perform pca linear dimension reduction
#> ℹ [2026-03-09 08:30:29] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-03-09 08:30:29] Reorder clusters...
#> ℹ [2026-03-09 08:30:29] Perform umap nonlinear dimension reduction
#> ℹ [2026-03-09 08:30:29] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-03-09 08:30:34] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-03-09 08:30:39] Run scop standard workflow completed
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
#> ℹ [2026-03-09 08:30:41] Start find dynamic features
#> ℹ [2026-03-09 08:30:41] Data type is raw counts
#> ℹ [2026-03-09 08:30:43] Number of candidate features (union): 231
#> ℹ [2026-03-09 08:30:44] Data type is raw counts
#> ℹ [2026-03-09 08:30:44] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-03-09 08:30:44] Using 1 core
#> ⠙ [2026-03-09 08:30:44] Running for Gcg [1/231] ■                              …
#> ⠹ [2026-03-09 08:30:44] Running for Cks2 [37/231] ■■■■■■                       …
#> ⠸ [2026-03-09 08:30:44] Running for Ass1 [118/231] ■■■■■■■■■■■■■■■■            …
#> ⠼ [2026-03-09 08:30:44] Running for Ins2 [200/231] ■■■■■■■■■■■■■■■■■■■■■■■■■■■ …
#> ✔ [2026-03-09 08:30:44] Completed 231 tasks in 8.7s
#> 
#> ℹ [2026-03-09 08:30:44] Building results
#> ℹ [2026-03-09 08:30:52] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-03-09 08:30:52] Using 1 core
#> ⠙ [2026-03-09 08:30:52] Running for Peg10 [44/231] ■■■■■■■                     …
#> ⠹ [2026-03-09 08:30:52] Running for Smc4 [117/231] ■■■■■■■■■■■■■■■■            …
#> ⠸ [2026-03-09 08:30:52] Running for Ckb [192/231] ■■■■■■■■■■■■■■■■■■■■■■■■■■   …
#> ✔ [2026-03-09 08:30:52] Completed 231 tasks in 9.1s
#> 
#> ℹ [2026-03-09 08:30:52] Building results
#> ✔ [2026-03-09 08:31:01] Find dynamic features done

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
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-03-09 08:31:02] [1] 180 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Ghrl,Iapp,Pyy,Rbp4,Lrpprc,Slc38a5,Cdkn1a,2810417H13Rik,Chga...
#> ℹ [2026-03-09 08:31:03] 
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
#> ℹ [2026-03-09 08:31:05] Start find dynamic features
#> ℹ [2026-03-09 08:31:06] Data type is raw counts
#> ℹ [2026-03-09 08:31:07] Number of candidate features (union): 2
#> ℹ [2026-03-09 08:31:07] Data type is raw counts
#> ℹ [2026-03-09 08:31:07] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-03-09 08:31:07] Using 1 core
#> ⠙ [2026-03-09 08:31:07] Running for Arxes1 [1/2] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-03-09 08:31:07] Completed 2 tasks in 164ms
#> 
#> ℹ [2026-03-09 08:31:07] Building results
#> ✔ [2026-03-09 08:31:07] Find dynamic features done
#> ℹ [2026-03-09 08:31:07] Start find dynamic features
#> ℹ [2026-03-09 08:31:08] Data type is raw counts
#> ℹ [2026-03-09 08:31:08] Number of candidate features (union): 2
#> ℹ [2026-03-09 08:31:09] Data type is raw counts
#> ℹ [2026-03-09 08:31:09] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-03-09 08:31:09] Using 1 core
#> ⠙ [2026-03-09 08:31:09] Running for Arxes1 [1/2] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-03-09 08:31:09] Completed 2 tasks in 237ms
#> 
#> ℹ [2026-03-09 08:31:09] Building results
#> ✔ [2026-03-09 08:31:09] Find dynamic features done


pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200,
  fit_method = "pretsa"
)
#> ℹ [2026-03-09 08:31:10] Start find dynamic features
#> ℹ [2026-03-09 08:31:11] Data type is raw counts
#> ℹ [2026-03-09 08:31:13] Number of candidate features (union): 231
#> ℹ [2026-03-09 08:31:14] Data type is raw counts
#> ℹ [2026-03-09 08:31:14] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-03-09 08:31:14] Calculating dynamic features for "Lineage2"...
#> ✔ [2026-03-09 08:31:14] Find dynamic features done
head(
  pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
)
#>      features exp_ncells      r.sq  dev.expl peaktime valleytime        pvalue
#> Gcg       Gcg        182 0.5551612 0.5551612 21.33789 12.6658731 6.319785e-128
#> Ghrl     Ghrl        167 0.1285034 0.1285034 18.32857  5.0984726  1.228737e-21
#> Iapp     Iapp        279 0.6477017 0.6477017 21.33789  0.1164741 7.284924e-165
#> Pyy       Pyy        434 0.6877808 0.6877808 21.33789  7.2132681 5.380374e-184
#> Rbp4     Rbp4        396 0.6305632 0.6305632 21.33789  6.8807704 2.434637e-157
#> Gast     Gast         92 0.2618876 0.2618876 21.33789  8.3787085  8.118511e-48
#>            padjust
#> Gcg  2.517018e-127
#> Ghrl  1.669636e-21
#> Iapp 6.232657e-164
#> Pyy  7.767915e-183
#> Rbp4 1.654121e-156
#> Gast  1.330054e-47
ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  cell_annotation = "SubCellType",
  n_split = 3,
  reverse_ht = "Lineage1"
)
#> ℹ [2026-03-09 08:31:14] [1] 168 features from Lineage1,Lineage2 passed the threshold (exp_ncells>[1] 20 & r.sq>[1] 0.2 & dev.expl>[1] 0.2 & padjust<[1] 0.05): 
#> ℹ                       Gcg,Iapp,Pyy,Rbp4,Gast,Chgb,Lrpprc,Slc38a5,Cdkn1a,2810417H13Rik...
#> ℹ [2026-03-09 08:31:15] 
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
#> ℹ [2026-03-09 08:31:18] Start find dynamic features
#> ℹ [2026-03-09 08:31:18] Data type is raw counts
#> ℹ [2026-03-09 08:31:20] Number of candidate features (union): 2
#> ℹ [2026-03-09 08:31:21] Data type is raw counts
#> ℹ [2026-03-09 08:31:21] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-03-09 08:31:21] Using 1 core
#> ⠙ [2026-03-09 08:31:21] Running for Arxes1 [1/2] ■■■■■■■■■■■■■■■■              …
#> ✔ [2026-03-09 08:31:21] Completed 2 tasks in 161ms
#> 
#> ℹ [2026-03-09 08:31:21] Building results
#> ✔ [2026-03-09 08:31:21] Find dynamic features done
#> ℹ [2026-03-09 08:31:21] Start find dynamic features
#> ℹ [2026-03-09 08:31:22] Data type is raw counts
#> ℹ [2026-03-09 08:31:22] Number of candidate features (union): 2
#> ℹ [2026-03-09 08:31:23] Data type is raw counts
#> ℹ [2026-03-09 08:31:23] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-03-09 08:31:23] Using 1 core
#> ℹ [2026-03-09 08:31:23] Building results
#> ✔ [2026-03-09 08:31:23] Find dynamic features done
```
