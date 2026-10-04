# Plot CytoTRACE 2 Results

Plot CytoTRACE 2 Results

## Usage

``` r
CytoTRACEPlot(
  object,
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
  ...,
  srt = NULL
)
```

## Arguments

- object:

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

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

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
#> ℹ [2026-10-04 21:42:43] Start standard processing workflow...
#> ℹ [2026-10-04 21:42:43] Checking a list of <Seurat>...
#> ! [2026-10-04 21:42:43] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-04 21:42:43] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 21:42:43] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 21:42:43] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 21:42:43] Number of available HVF: 2000
#> ℹ [2026-10-04 21:42:43] Finished check
#> ℹ [2026-10-04 21:42:43] Perform `ScaleData()`
#> ℹ [2026-10-04 21:42:44] Perform pca linear dimension reduction
#> ℹ [2026-10-04 21:42:44] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-04 21:42:44] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 21:42:44] Reorder clusters...
#> ℹ [2026-10-04 21:42:44] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 21:42:44] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-04 21:42:50] Standard processing workflow completed
pancreas_sub <- RunCytoTRACE(
  pancreas_sub,
  species = "Mus_musculus",
  backend = "cpp"
)
#> ◌ [2026-10-04 21:42:50] Running CytoTRACE2
#> ℹ [2026-10-04 21:42:50] Extracting expression matrix from `assay = RNA, layer = counts`
#> ◌ [2026-10-04 21:42:50] Running CytoTRACE2 with `backend = cpp`
#> ℹ [2026-10-04 21:42:50] Preparing CytoTRACE2 database
#> ℹ [2026-10-04 21:42:50] Downloading CytoTRACE2 model data from datasets GitHub repository...
#> ℹ [2026-10-04 21:42:50]   Downloading model_parameters.rds ...
#> ℹ [2026-10-04 21:42:52]   Downloading features_model_training_17.csv ...
#> ℹ [2026-10-04 21:42:52]   Downloading mt_dict_human_to_mouse.csv ...
#> ℹ [2026-10-04 21:42:52]   Downloading mt_human_alias.csv ...
#> ℹ [2026-10-04 21:42:52]   Downloading mt_mouse_alias.csv ...
#> ✔ [2026-10-04 21:42:52] CytoTRACE2 data cached at /home/runner/.local/share/R/scop/CytoTRACE2
#> ℹ [2026-10-04 21:42:52] Loading model from /home/runner/.local/share/R/scop/CytoTRACE2
#> ℹ [2026-10-04 21:42:55] Dataset contains 15998 genes and 1000 cells.
#> ℹ [2026-10-04 21:42:55] Running on 1 subsample
#> ℹ [2026-10-04 21:42:55] Using 1 core
#> ℹ [2026-10-04 21:42:55] 12486 input genes mapped to model genes.
#> ℹ [2026-10-04 21:42:55] Building results
#> ✔ [2026-10-04 21:43:04] CytoTRACE2 computed successfully

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
