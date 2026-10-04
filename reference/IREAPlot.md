# Plot IREA response or polarization results

Display enrichment effects and reference-cell significance.

## Usage

``` r
IREAPlot(
  object,
  plot_type = c("compass", "radar", "dotplot", "heatmap"),
  padjustCutoff = 0.05,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop"
)
```

## Arguments

- object:

  An `irea_result`, a Seurat object containing `object@tools$IREA`, or a
  named list of results for dotplots and heatmaps.

- plot_type:

  One of `"compass"`, `"radar"`, `"dotplot"`, or `"heatmap"`.

- padjustCutoff:

  Significance threshold used for compass colouring and the radar
  display. Radar scores are zero when no positive effect meets this
  threshold; otherwise all positive effects are divided by their
  maximum.

- palette, palcolor:

  Palette name and optional custom colours passed to
  [`thisplot::palette_colors()`](https://mengxu98.github.io/thisplot/reference/palette_colors.html).
  For compass plots, colours map to Positive, Negative and Not
  significant, in that order.

- theme_use:

  Theme name, function or ggplot theme, as in other package plots.

## Value

An editable ggplot object; plotting does not change result statistics.

## References

Cui, Ang; Huang, Teddy; Li, Shuqiang; Ma, Aileen; Perez, Jorge L.;
Sander, Chris; Keskin, Derin B.; Wu, Catherine J.; Fraenkel, Ernest;
Hacohen, Nir. Dictionary of immune responses to cytokines at single-cell
resolution. Nature 625, 377-384 (2024). doi:10.1038/s41586-023-06816-9.

## Examples

``` r
if (FALSE) { # \dontrun{
reference <- PrepareDB(
  db = "IREA_Macrophage", species = "Homo_sapiens"
)[["Homo_sapiens"]][["IREA_Macrophage"]]
result <- RunIREA(c("ISG15", "IFIT3", "BST2"), reference = reference)
IREAPlot(result, palette = "Chinese")
} # }
```
