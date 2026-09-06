# Plot scPagwas Scores

Plot the score and adjusted p-value metadata produced by
[`RunscPagwas()`](https://mengxu98.github.io/scop/reference/RunscPagwas.md)
on an existing Seurat reduction. The plots are returned and can
optionally be saved as PDF files.

## Usage

``` r
PlotscPagwas(
  srt,
  reduction = c("umap", "tsne"),
  features = NULL,
  p_threshold = 0.05,
  output.dir = NULL,
  width = 7,
  height = 7,
  point_size = NULL,
  palette = "Spectral",
  palcolor = NULL,
  significance_palette = "Chinese",
  significance_palcolor = NULL,
  do_plot = TRUE
)
```

## Arguments

- srt:

  A Seurat object returned by
  [`RunscPagwas()`](https://mengxu98.github.io/scop/reference/RunscPagwas.md).

- reduction:

  Reduction used for plotting, either `"umap"` or `"tsne"`.

- features:

  Numeric scPagwas metadata columns to plot. By default, available gPAS,
  TRS, and down-TRS score columns are used.

- p_threshold:

  Adjusted p-value threshold used to identify significant cells. Set to
  `NULL` to omit the significance plot.

- output.dir:

  Optional directory in which to save PDF files.

- width, height:

  PDF dimensions in inches.

- point_size:

  Point size passed to
  [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)
  and
  [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md).

- palette, palcolor:

  Palette used for continuous scPagwas scores, passed to
  [`FeatureDimPlot()`](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md).

- significance_palette, significance_palcolor:

  Palette used for the significance groups, passed to
  [`CellDimPlot()`](https://mengxu98.github.io/scop/reference/CellDimPlot.md).

- do_plot:

  Whether to print each plot.

## Value

A named list of ggplot objects.
