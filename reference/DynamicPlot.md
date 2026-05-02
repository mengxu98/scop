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
  point_palette = "Chinese",
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
  cores = 1,
  verbose = TRUE,
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
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

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
  point colors. Default is `"Chinese"`.

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

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-02 04:37:04] Start standard processing workflow...
#> ℹ [2026-05-02 04:37:05] Checking a list of <Seurat>...
#> ! [2026-05-02 04:37:05] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-02 04:37:05] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-02 04:37:07] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-02 04:37:08] Use the separate HVF from `srt_list`
#> ℹ [2026-05-02 04:37:08] Number of available HVF: 2000
#> ℹ [2026-05-02 04:37:08] Finished check
#> ℹ [2026-05-02 04:37:08] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-02 04:37:08] Perform pca linear dimension reduction
#> ℹ [2026-05-02 04:37:09] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-02 04:37:09] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-02 04:37:09] Reorder clusters...
#> ℹ [2026-05-02 04:37:10] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-02 04:37:10] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-02 04:37:10] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-02 04:37:14] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-02 04:37:17] Standard processing workflow completed
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_path()`).
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_path()`).


CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1
)


DynamicPlot(
  pancreas_sub,
  lineages = "Lineage1",
  features = c("Arxes1", "Ncoa2", "G2M_score"),
  group.by = "SubCellType",
  compare_features = TRUE
)
#> ℹ [2026-05-02 04:37:20] Start find dynamic features
#> ℹ [2026-05-02 04:37:21] Data type is raw counts
#> ℹ [2026-05-02 04:37:22] Number of candidate features (union): 3
#> ℹ [2026-05-02 04:37:22] Data type is raw counts
#> ! [2026-05-02 04:37:22] Negative values detected
#> ℹ [2026-05-02 04:37:22] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-02 04:37:22] Using 1 core
#> ⠙ [2026-05-02 04:37:22] Running for Arxes1 [1/3] ■■■         33% | ETA:  0s
#> ✔ [2026-05-02 04:37:22] Completed 3 tasks in 151ms
#> 
#> ℹ [2026-05-02 04:37:22] Building results
#> ✔ [2026-05-02 04:37:22] Find dynamic features done


DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  features = c("Arxes1", "Ncoa2", "G2M_score"),
  group.by = "SubCellType",
  compare_lineages = TRUE,
  compare_features = FALSE
)
#> ℹ [2026-05-02 04:37:23] Start find dynamic features
#> ℹ [2026-05-02 04:37:24] Data type is raw counts
#> ℹ [2026-05-02 04:37:25] Number of candidate features (union): 3
#> ℹ [2026-05-02 04:37:25] Data type is raw counts
#> ! [2026-05-02 04:37:25] Negative values detected
#> ℹ [2026-05-02 04:37:25] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-02 04:37:25] Using 1 core
#> ⠙ [2026-05-02 04:37:25] Running for Arxes1 [1/3] ■■■         33% | ETA:  0s
#> ✔ [2026-05-02 04:37:25] Completed 3 tasks in 150ms
#> 
#> ℹ [2026-05-02 04:37:25] Building results
#> ✔ [2026-05-02 04:37:25] Find dynamic features done
#> ℹ [2026-05-02 04:37:25] Start find dynamic features
#> ℹ [2026-05-02 04:37:27] Data type is raw counts
#> ℹ [2026-05-02 04:37:27] Number of candidate features (union): 3
#> ℹ [2026-05-02 04:37:27] Data type is raw counts
#> ! [2026-05-02 04:37:27] Negative values detected
#> ℹ [2026-05-02 04:37:27] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-02 04:37:27] Using 1 core
#> ℹ [2026-05-02 04:37:27] Building results
#> ✔ [2026-05-02 04:37:28] Find dynamic features done
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
#> ℹ [2026-05-02 04:37:29] Start find dynamic features
#> ℹ [2026-05-02 04:37:30] Data type is raw counts
#> ℹ [2026-05-02 04:37:31] Number of candidate features (union): 3
#> ℹ [2026-05-02 04:37:31] Data type is raw counts
#> ! [2026-05-02 04:37:31] Negative values detected
#> ℹ [2026-05-02 04:37:31] Calculating dynamic features for "Lineage1"...
#> ℹ [2026-05-02 04:37:31] Using 1 core
#> ⠙ [2026-05-02 04:37:31] Running for Arxes1 [1/3] ■■■         33% | ETA:  0s
#> ✔ [2026-05-02 04:37:31] Completed 3 tasks in 353ms
#> 
#> ℹ [2026-05-02 04:37:31] Building results
#> ✔ [2026-05-02 04:37:32] Find dynamic features done
#> ℹ [2026-05-02 04:37:32] Start find dynamic features
#> ℹ [2026-05-02 04:37:33] Data type is raw counts
#> ℹ [2026-05-02 04:37:33] Number of candidate features (union): 3
#> ℹ [2026-05-02 04:37:34] Data type is raw counts
#> ! [2026-05-02 04:37:34] Negative values detected
#> ℹ [2026-05-02 04:37:34] Calculating dynamic features for "Lineage2"...
#> ℹ [2026-05-02 04:37:34] Using 1 core
#> ⠙ [2026-05-02 04:37:34] Running for Arxes1 [1/3] ■■■         33% | ETA:  0s
#> ✔ [2026-05-02 04:37:34] Completed 3 tasks in 182ms
#> 
#> ℹ [2026-05-02 04:37:34] Building results
#> ✔ [2026-05-02 04:37:34] Find dynamic features done
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
