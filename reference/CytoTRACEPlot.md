# Plot CytoTRACE 2 Results

Plot CytoTRACE 2 Results

## Usage

``` r
CytoTRACEPlot(
  srt,
  reduction = NULL,
  group.by = NULL,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- reduction:

  Reduction to plot. `NULL` uses
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- group.by:

  Metadata column(s) used to color cells.

- combine, nrow, ncol, byrow:

  Combine plots with
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html).
  `combine = FALSE` returns a list of ggplots.

- pt.size, pt.alpha:

  Point size and transparency. `pt.size = NULL` scales with `sqrt(n)`
  (minimum `0.3`). Rasterized points keep at least a two-pixel radius at
  `raster.dpi = c(512, 512)` and scale with `raster.dpi`.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- theme_use, theme_args:

  Theme name or function, plus extra theme arguments.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to
  [CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
  and
  [FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md).

## Value

If `combine = TRUE`, returns a `patchwork` object combining all plots.
If `combine = FALSE`, returns a named list of ggplot objects:

- `Score`: UMAP plot colored by score computed by CytoTRACE2;

- `Potency`: UMAP plot colored by potency category computed by
  CytoTRACE2;

- `Relative`: UMAP plot colored by relative score computed by
  CytoTRACE2;

- `Phenotype`: UMAP plot colored by phenotype (if `group.by` is
  provided);

- `Boxplot`: Boxplot of score computed by CytoTRACE2 corresponding to
  phenotype (if `group.by` is provided).

## See also

[RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 21:05:48] Start standard processing workflow...
#> ℹ [2026-09-06 21:05:49] Checking a list of <Seurat>...
#> ! [2026-09-06 21:05:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 21:05:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:05:49] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 21:05:49] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 21:05:49] Number of available HVF: 2000
#> ℹ [2026-09-06 21:05:49] Finished check
#> ℹ [2026-09-06 21:05:49] Perform `ScaleData()`
#> ℹ [2026-09-06 21:05:49] Perform pca linear dimension reduction
#> ℹ [2026-09-06 21:05:49] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 21:05:50] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 21:05:50] Reorder clusters...
#> ℹ [2026-09-06 21:05:50] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 21:05:50] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 21:05:56] Standard processing workflow completed
pancreas_sub <- RunCytoTRACE(
  pancreas_sub,
  species = "Mus_musculus",
  backend = "cpp"
)
#> ◌ [2026-09-06 21:05:56] Running CytoTRACE2
#> ℹ [2026-09-06 21:05:56] Extracting expression matrix from `assay = RNA, layer = counts`
#> ◌ [2026-09-06 21:05:56] Running CytoTRACE2 with `backend = cpp`
#> ℹ [2026-09-06 21:05:56] Preparing CytoTRACE2 database
#> ℹ [2026-09-06 21:05:56] Downloading CytoTRACE2 model data from datasets GitHub repository...
#> ℹ [2026-09-06 21:05:56]   Downloading model_parameters.rds ...
#> ℹ [2026-09-06 21:05:57]   Downloading features_model_training_17.csv ...
#> ℹ [2026-09-06 21:05:57]   Downloading mt_dict_human_to_mouse.csv ...
#> ℹ [2026-09-06 21:05:58]   Downloading mt_human_alias.csv ...
#> ℹ [2026-09-06 21:05:58]   Downloading mt_mouse_alias.csv ...
#> ✔ [2026-09-06 21:05:58] CytoTRACE2 data cached at /home/runner/.local/share/R/scop/CytoTRACE2
#> ℹ [2026-09-06 21:05:58] Species: "Homo_sapiens"
#> ℹ [2026-09-06 21:05:58] Species: "Mus_musculus"
#> ℹ [2026-09-06 21:05:58] Loading model from /home/runner/.local/share/R/scop/CytoTRACE2
#> ℹ [2026-09-06 21:06:01] Dataset contains 15998 genes and 1000 cells.
#> ℹ [2026-09-06 21:06:01] Running on 1 subsample
#> ℹ [2026-09-06 21:06:01] Using 1 core
#> ℹ [2026-09-06 21:06:01] 12486 input genes mapped to model genes.
#> ℹ [2026-09-06 21:06:01] Building results
#> ✔ [2026-09-06 21:06:20] CytoTRACE2 computed successfully

CytoTRACEPlot(
  pancreas_sub,
  group.by = "CellType",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)


plots <- CytoTRACEPlot(
  pancreas_sub,
  group.by = "CellType",
  combine = FALSE
)
plots$Boxplot
```
