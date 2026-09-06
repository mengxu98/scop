# Plot stored spatial deconvolution proportions

Plot spot-by-cell-type proportions stored by
[`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md),
[`RunCARD()`](https://mengxu98.github.io/scop/reference/RunCARD.md),
[`RunSPOTlight()`](https://mengxu98.github.io/scop/reference/RunSPOTlight.md),
or
[`RunSpatialDWLS()`](https://mengxu98.github.io/scop/reference/RunSpatialDWLS.md).
The plot reads the stored result directly from `srt@tools[[tool_name]]`
and never reruns a deconvolution backend.
[`RunCSIDE()`](https://mengxu98.github.io/scop/reference/RunCSIDE.md) is
intentionally excluded because its output represents differential or
context effects rather than cell-type proportions.

## Usage

``` r
SpatialDeconvolutionPlot(
  srt,
  tool_name = NULL,
  cell_types = NULL,
  plot_type = c("point", "dominant", "pie"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...,
  image.scale = c("lowres", "hires")
)
```

## Arguments

- srt:

  A spatial `Seurat` object containing a stored deconvolution result.

- tool_name:

  Explicit non-empty key in `srt@tools`. Results are never discovered
  implicitly.

- cell_types:

  Optional cell types to display. The default uses all stored cell
  types.

- plot_type:

  Plot proportions as separate point maps, one dominant-type map derived
  from the stored proportions, or one spot-level pie map.

- combine:

  Whether to combine point maps. If `FALSE`, return a named list.

- nrow, ncol, byrow:

  Point-map layout controls. When both dimensions are `NULL`, a
  near-square layout with at most three columns is used.

- ...:

  Additional arguments passed to
  [`SpatialSpotPlot()`](https://mengxu98.github.io/scop/reference/SpatialSpotPlot.md).

- image.scale:

  Image scale factor matching the selected raster.

## Value

A `ggplot`, `patchwork`, or named list of `ggplot` objects.

## Examples

``` r
data(visium_human_pancreas_sub)
data(panc8_sub)
keep_spots <- unique(round(seq(1, ncol(visium_human_pancreas_sub), length.out = 120)))
spatial <- visium_human_pancreas_sub[, keep_spots]
reference <- panc8_sub[, panc8_sub@meta.data[["celltype"]] %in%
  c("ductal", "alpha", "beta")]
reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
shared <- head(intersect(
  SeuratObject::VariableFeatures(reference),
  rownames(spatial)
), 300)
spatial <- RunSpatialDWLS(
  srt = spatial,
  reference = reference,
  reference_label = "celltype",
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  features = shared,
  coord.cols = c("x", "y"),
  min_cells = 2,
  verbose = FALSE
)
SpatialDeconvolutionPlot(
  spatial,
  tool_name = "SpatialDWLS",
  plot_type = "dominant",
  overlay_image = FALSE,
  coord.cols = c("x", "y")
)
```
