# Features correlation plot

Pairwise feature correlations in a `Seurat` object.

## Usage

``` r
FeatureCorPlot(
  srt,
  features,
  group.by = NULL,
  split.by = NULL,
  cells = NULL,
  layer = "data",
  assay = NULL,
  cor_method = "pearson",
  adjust = 1,
  margin = 1,
  reverse = FALSE,
  add_equation = FALSE,
  add_r2 = TRUE,
  add_pvalue = TRUE,
  add_smooth = TRUE,
  palette = "Chinese",
  palcolor = NULL,
  cor_palette = "RdBu",
  cor_palcolor = NULL,
  cor_range = c(-1, 1),
  pt.size = NULL,
  pt.alpha = 1,
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  calculate_coexp = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  aspect.ratio = 1,
  title = NULL,
  subtitle = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  theme_use = "theme_scop",
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- features:

  Features to compare. Should be present in both the assay data and the
  metadata of the Seurat object.

- group.by:

  Metadata column(s) used to color cells.

- split.by:

  Metadata column to facet by.

- cells:

  Cell names to include.

- layer:

  Assay layer to use.

- assay:

  Assay to use. `NULL` uses the default assay.

- cor_method:

  `"pearson"` or `"spearman"`.

- adjust:

  The adjustment factor for the width of the violin plots.

- margin:

  The margin size for the plot.

- reverse:

  Whether to reverse the order of the features in the plot.

- add_equation:

  Whether to add the equation of the linear regression line to each
  scatter plot.

- add_r2:

  Whether to add the R-squared value of the linear regression line to
  each scatter plot.

- add_pvalue:

  Whether to add the p-value of the linear regression line to each
  scatter plot.

- add_smooth:

  Whether to add a smoothed line to each scatter plot.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- cor_palette:

  Name of the color palette to use for the correlation.

- cor_palcolor:

  Color for the correlation.

- cor_range:

  A two-length numeric vector specifying the range for the correlation.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- cells.highlight, cols.highlight, sizes.highlight, alpha.highlight,
  stroke.highlight:

  Cells to highlight and their appearance. `TRUE` highlights all cells.

- calculate_coexp:

  Whether to calculate the co-expression of selected features.

- raster, raster.dpi:

  Rasterize points. `raster = NULL` rasterizes when there are more than
  100,000 cells.

- aspect.ratio:

  Panel aspect ratio.

- title:

  Plot title. `NULL` hides the title for merged/single panels. When
  multiple lineages are plotted and `title` is `NULL`, each panel is
  titled with its lineage column.

- subtitle:

  Plot subtitle.

- legend.position:

  Legend placement (`"none"`, `"left"`, `"right"`, `"bottom"`, `"top"`),
  direction, and title. `legend.title = NULL` uses the group name.

- legend.direction:

  Legend direction: `"horizontal"` or `"vertical"`.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- force:

  Whether to force the creation of the plot, even if it contains more
  than 50 subplots.

- seed:

  Random seed.

- verbose:

  Whether to print messages.

## See also

[FeatureStatPlot](https://mengxu98.github.io/scop/reference/FeatureStatPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:28:25] Start standard processing workflow...
#> ℹ [2026-09-06 21:28:26] Checking a list of <Seurat>...
#> ! [2026-09-06 21:28:26] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:28:26] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:28:26] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:28:26] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:28:26] Number of available HVF: 2000
#> ℹ [2026-09-06 21:28:26] Finished check
#> ℹ [2026-09-06 21:28:26] Perform `ScaleData()`
#> ℹ [2026-09-06 21:28:26] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:28:26] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:28:27] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:28:27] Reorder clusters...
#> ℹ [2026-09-06 21:28:27] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:28:27] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:28:33] Standard processing workflow completed
FeatureCorPlot(
  pancreas_sub,
  features = rownames(pancreas_sub)[1:5],
  group.by = "SubCellType"
)


FeatureCorPlot(
  pancreas_sub,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "nFeature_spliced",
    "nCount_spliced",
    "nFeature_unspliced",
    "nCount_unspliced"
  ),
  group.by = "SubCellType",
  cor_palette = "Greys",
  cor_range = c(0, 1)
)


FeatureCorPlot(
  pancreas_sub,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "SubCellType",
  add_equation = TRUE
)
```
