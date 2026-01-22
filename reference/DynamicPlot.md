# Plot dynamic features across pseudotime

Plot dynamic features across pseudotime

## Usage

``` r
DynamicPlot(
  srt,
  lineages,
  features,
  group.by = NULL,
  cells = NULL,
  layer = "counts",
  assay = NULL,
  family = NULL,
  exp_method = c("log1p", "raw", "zscore", "fc", "log2fc"),
  lib_normalize = identical(layer, "counts"),
  libsize = NULL,
  compare_lineages = TRUE,
  compare_features = FALSE,
  add_line = TRUE,
  add_interval = TRUE,
  line.size = 1,
  line_palette = "Dark2",
  line_palcolor = NULL,
  add_point = TRUE,
  pt.size = 1,
  point_palette = "Paired",
  point_palcolor = NULL,
  add_rug = TRUE,
  flip = FALSE,
  reverse = FALSE,
  x_order = c("value", "rank"),
  aspect.ratio = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- lineages:

  A character vector specifying the lineages to plot.

- features:

  A character vector of features to use.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- cells:

  A character vector of cell names to use. Default is `NULL`.

- layer:

  Which layer to use. Default is `"counts"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- family:

  A character specifying the model used to calculate the dynamic
  features if needed. By default, this parameter is set to `NULL`, and
  the appropriate family will be automatically determined.

- exp_method:

  A character specifying the method to transform the expression values.
  Default is `"log1p"` with options "log1p", "raw", "zscore", "fc",
  "log2fc".

- lib_normalize:

  A boolean specifying whether to normalize the expression values using
  library size. Default the `layer` is counts, this parameter is set to
  `TRUE`. Otherwise, it is set to `FALSE`.

- libsize:

  A numeric vector specifying the library size for each cell. Default is
  `NULL`.

- compare_lineages:

  A boolean specifying whether to compare the lineages in the plot.
  Default is `TRUE`.

- compare_features:

  A boolean specifying whether to compare the features in the plot.
  Default is `FALSE`.

- add_line:

  A boolean specifying whether to add lines to the plot. Default is
  `TRUE`.

- add_interval:

  A boolean specifying whether to add confidence intervals to the plot.
  Default is `TRUE`.

- line.size:

  A numeric specifying the size of the lines. Default is `1`.

- line_palette:

  A character string specifying the name of the palette to use for the
  line colors. Default is `"Dark2"`.

- line_palcolor:

  A vector specifying the colors to use for the line palette. Default is
  `NULL`.

- add_point:

  A boolean specifying whether to add points to the plot. Default is
  `TRUE`.

- pt.size:

  A numeric specifying the size of the points. Default is `1`.

- point_palette:

  A character string specifying the name of the palette to use for the
  point colors. Default is `"Paired"`.

- point_palcolor:

  A vector specifying the colors to use for the point palette. Default
  is `NULL`.

- add_rug:

  A boolean specifying whether to add rugs to the plot. Default is
  `TRUE`.

- flip:

  A boolean specifying whether to flip the x-axis. Default is `FALSE`.

- reverse:

  A boolean specifying whether to reverse the x-axis. Default is
  `FALSE`.

- x_order:

  A character specifying the order of the x-axis values. Default is
  `c("value", "rank")`.

- aspect.ratio:

  Aspect ratio of the panel. Default is `NULL`.

- legend.position:

  The position of legends, one of `"none"`, `"left"`, `"right"`,
  `"bottom"`, `"top"`. Default is `"right"`.

- legend.direction:

  The direction of the legend in the plot. Can be one of `"vertical"` or
  `"horizontal"`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Combine plots into a single `patchwork` object. If `FALSE`, return a
  list of ggplot objects.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-01-22 03:33:47] Start standard scop workflow...
#> ℹ [2026-01-22 03:33:48] Checking a list of <Seurat>...
#> ! [2026-01-22 03:33:48] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-01-22 03:33:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 03:33:50] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2026-01-22 03:33:51] Use the separate HVF from srt_list
#> ℹ [2026-01-22 03:33:51] Number of available HVF: 2000
#> ℹ [2026-01-22 03:33:51] Finished check
#> ℹ [2026-01-22 03:33:51] Perform `Seurat::ScaleData()`
#> ℹ [2026-01-22 03:33:51] Perform pca linear dimension reduction
#> ℹ [2026-01-22 03:33:52] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-01-22 03:33:52] Reorder clusters...
#> ℹ [2026-01-22 03:33:52] Perform umap nonlinear dimension reduction
#> ℹ [2026-01-22 03:33:52] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2026-01-22 03:33:56] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ✔ [2026-01-22 03:34:00] Run scop standard workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_path()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1
)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


DynamicPlot(
  pancreas_sub,
  lineages = "Lineage1",
  features = c("Arxes1", "Ncoa2", "G2M_score"),
  group.by = "SubCellType",
  compare_features = TRUE
)
#> ℹ [2026-01-22 03:34:02] Start find dynamic features
#> ℹ [2026-01-22 03:34:03] Data type is raw counts
#> ℹ [2026-01-22 03:34:03] Number of candidate features (union): 3
#> ℹ [2026-01-22 03:34:03] Data type is raw counts
#> ! [2026-01-22 03:34:03] Negative values detected
#> ℹ [2026-01-22 03:34:03] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-22 03:34:03] Using 1 core
#> ⠙ [2026-01-22 03:34:03] Running for Arxes1 [1/3] ■■■■■■■■■■■                   …
#> ✔ [2026-01-22 03:34:03] Completed 3 tasks in 130ms
#> 
#> ℹ [2026-01-22 03:34:03] Building results
#> ✔ [2026-01-22 03:34:04] Find dynamic features done


DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2", "G2M_score"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ℹ [2026-01-22 03:34:04] Start find dynamic features
#> ℹ [2026-01-22 03:34:05] Data type is raw counts
#> ℹ [2026-01-22 03:34:06] Number of candidate features (union): 3
#> ℹ [2026-01-22 03:34:07] Data type is raw counts
#> ! [2026-01-22 03:34:07] Negative values detected
#> ℹ [2026-01-22 03:34:07] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-22 03:34:07] Using 1 core
#> ⠙ [2026-01-22 03:34:07] Running for Arxes1 [1/3] ■■■■■■■■■■■                   …
#> ✔ [2026-01-22 03:34:07] Completed 3 tasks in 130ms
#> 
#> ℹ [2026-01-22 03:34:07] Building results
#> ✔ [2026-01-22 03:34:07] Find dynamic features done
#> ℹ [2026-01-22 03:34:07] Start find dynamic features
#> ℹ [2026-01-22 03:34:07] Data type is raw counts
#> ℹ [2026-01-22 03:34:08] Number of candidate features (union): 3
#> ℹ [2026-01-22 03:34:08] Data type is raw counts
#> ! [2026-01-22 03:34:08] Negative values detected
#> ℹ [2026-01-22 03:34:08] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-01-22 03:34:08] Using 1 core
#> ⠙ [2026-01-22 03:34:08] Running for Arxes1 [1/3] ■■■■■■■■■■■                   …
#> ✔ [2026-01-22 03:34:08] Completed 3 tasks in 103ms
#> 
#> ℹ [2026-01-22 03:34:08] Building results
#> ✔ [2026-01-22 03:34:08] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2", "G2M_score"),
  group.by = "SubCellType",
  compare_lineages = FALSE,
  compare_features = FALSE
)
#> ℹ [2026-01-22 03:34:10] Start find dynamic features
#> ℹ [2026-01-22 03:34:10] Data type is raw counts
#> ℹ [2026-01-22 03:34:11] Number of candidate features (union): 3
#> ℹ [2026-01-22 03:34:11] Data type is raw counts
#> ! [2026-01-22 03:34:11] Negative values detected
#> ℹ [2026-01-22 03:34:11] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-01-22 03:34:11] Using 1 core
#> ⠙ [2026-01-22 03:34:11] Running for Arxes1 [1/3] ■■■■■■■■■■■                   …
#> ✔ [2026-01-22 03:34:11] Completed 3 tasks in 133ms
#> 
#> ℹ [2026-01-22 03:34:11] Building results
#> ✔ [2026-01-22 03:34:11] Find dynamic features done
#> ℹ [2026-01-22 03:34:11] Start find dynamic features
#> ℹ [2026-01-22 03:34:12] Data type is raw counts
#> ℹ [2026-01-22 03:34:12] Number of candidate features (union): 3
#> ℹ [2026-01-22 03:34:13] Data type is raw counts
#> ! [2026-01-22 03:34:13] Negative values detected
#> ℹ [2026-01-22 03:34:13] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-01-22 03:34:13] Using 1 core
#> ℹ [2026-01-22 03:34:13] Building results
#> ✔ [2026-01-22 03:34:13] Find dynamic features done
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
