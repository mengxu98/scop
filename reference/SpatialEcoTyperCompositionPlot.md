# SpatialEcoTyper composition plot

Plot SpatialEcoTyper composition across cell types, samples, or other
metadata groups using the default plotting theme and palettes.

## Usage

``` r
SpatialEcoTyperCompositionPlot(
  srt,
  se.by = "SpatialEcoTyper_SE",
  group.by = "CellType",
  sample.by = NULL,
  position = c("fill", "stack"),
  palette = "Paired",
  palcolor = NULL,
  legend.position = "right",
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- se.by:

  Metadata column containing SpatialEcoTyper labels.

- group.by:

  Metadata column used as the composition group.

- sample.by:

  Optional metadata column used for faceting.

- position:

  Bar position. `"fill"` shows fractions and `"stack"` shows counts.

- palette, palcolor:

  Palette passed to `palette_colors()`.

- legend.position:

  Legend position.

- theme_use:

  Theme function name.

- theme_args:

  Additional arguments passed to the theme function.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `ggplot` object.

## Examples

``` r
data(visium_human_pancreas_pair_sub)
# This fixture stores real PRECAST domains for a fast plot demonstration;
# use `se.by = "domain"` explicitly because it is not a new SE inference.
SpatialEcoTyperCompositionPlot(
  visium_human_pancreas_pair_sub,
  se.by = "domain",
  group.by = "coda_label",
  sample.by = "sample",
  position = "fill"
)
```
