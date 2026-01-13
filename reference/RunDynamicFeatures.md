# RunDynamicFeatures

Calculates dynamic features for lineages in a single-cell RNA-seq
dataset.

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
#> ⠙ [2026-01-13 09:17:32] Running for 1 [1/231] ■                                …
#> ⠹ [2026-01-13 09:17:32] Running for 36 [36/231] ■■■■■■                         …
#> ⠸ [2026-01-13 09:17:32] Running for 89 [89/231] ■■■■■■■■■■■■■                  …
#> ⠼ [2026-01-13 09:17:32] Running for 178 [178/231] ■■■■■■■■■■■■■■■■■■■■■■■■     …
#> ✔ [2026-01-13 09:17:32] Completed 231 tasks in 9s
#> 
#> ⠙ [2026-01-13 09:17:42] Running for 34 [34/231] ■■■■■                          …
#> ⠹ [2026-01-13 09:17:42] Running for 114 [114/231] ■■■■■■■■■■■■■■■■             …
#> ⠸ [2026-01-13 09:17:42] Running for 194 [194/231] ■■■■■■■■■■■■■■■■■■■■■■■■■■   …
#> ✔ [2026-01-13 09:17:42] Completed 231 tasks in 8.3s
#> 

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
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.

ht$plot


DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ⠙ [2026-01-13 09:17:55] Running for 1 [1/2] ■■■■■■■■■■■■■■■■                  5…
#> ✔ [2026-01-13 09:17:55] Completed 2 tasks in 135ms
#> 
#> ⠙ [2026-01-13 09:17:57] Running for 1 [1/2] ■■■■■■■■■■■■■■■■                  5…
#> ✔ [2026-01-13 09:17:57] Completed 2 tasks in 110ms
#> 
```
