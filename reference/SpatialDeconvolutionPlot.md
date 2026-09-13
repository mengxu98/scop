# Plot stored spatial deconvolution proportions

Plot spot-by-cell-type proportions stored by
[`RunRCTD()`](https://mengxu98.github.io/scop/reference/RunRCTD.md),
[`RunCARD()`](https://mengxu98.github.io/scop/reference/RunCARD.md),
[`RunSPOTlight()`](https://mengxu98.github.io/scop/reference/RunSPOTlight.md).
The plot reads the stored result directly from `srt@tools[[tool_name]]`
and never reruns a deconvolution backend.
[`RunCSIDE()`](https://mengxu98.github.io/scop/reference/RunCSIDE.md) is
intentionally excluded because its output represents differential or
context effects rather than cell-type proportions.

## Usage

``` r
SpatialDeconvolutionPlot(
  object,
  tool_name = NULL,
  cell_types = NULL,
  plot_type = c("point", "dominant", "pie"),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...,
  image.scale = c("lowres", "hires"),
  srt = NULL
)
```

## Arguments

- object:

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

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `ggplot`, `patchwork`, or named list of `ggplot` objects.

## Examples

``` r
if (FALSE) { # \dontrun{
data(visium_human_pancreas_sub)
data(panc8_sub)
reference <- panc8_sub[, panc8_sub$celltype %in% c("ductal", "alpha", "beta")]
reference <- Seurat::FindVariableFeatures(reference, nfeatures = 300, verbose = FALSE)
features <- head(intersect(
  SeuratObject::VariableFeatures(reference),
  rownames(visium_human_pancreas_sub)
), 300)
spatial <- RunRCTD(
  visium_human_pancreas_sub,
  reference = reference,
  reference_label = "celltype",
  features = features,
  assay = "Spatial",
  reference_assay = "RNA",
  layer = "counts",
  reference_layer = "counts",
  rctd_mode = "full",
  max_cores = 1,
  verbose = FALSE
)
SpatialDeconvolutionPlot(spatial, tool_name = "RCTD", plot_type = "dominant")
} # }
```
