# ggplot2 theme for spatial plots

Built on
[`theme_scop()`](https://mengxu98.github.io/scop/reference/theme_scop.md)
with axes, ticks, and panel grid hidden by default. Spatial plot helpers
default to `theme_use = "theme_spatial"`.

## Usage

``` r
theme_spatial(show_axes = FALSE, aspect.ratio = 1, base_size = 12, ...)
```

## Arguments

- show_axes:

  Whether to keep axis titles, text, ticks, and grid lines.

- aspect.ratio:

  Aspect ratio of the panel.

- base_size:

  Base font size

- ...:

  Arguments passed to the
  [ggplot2::theme](https://ggplot2.tidyverse.org/reference/theme.html).

## Value

A ggplot2 theme object (class `theme`, `gg`).

## See also

[`theme_scop()`](https://mengxu98.github.io/scop/reference/theme_scop.md)

## Examples

``` r
theme_spatial()
#> <theme> List of 22
#>  $ text                 : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : chr "black"
#>   ..@ size         : num 12
#>   ..@ hjust        : NULL
#>   ..@ vjust        : NULL
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : NULL
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ aspect.ratio         : num 1
#>  $ axis.title           : <ggplot2::element_blank>
#>  $ axis.text            : <ggplot2::element_blank>
#>  $ axis.line            : <ggplot2::element_blank>
#>  $ legend.background    : <ggplot2::element_blank>
#>  $ legend.key           : <ggplot2::element_rect>
#>   ..@ fill         : chr "transparent"
#>   ..@ colour       : chr "transparent"
#>   ..@ linewidth    : NULL
#>   ..@ linetype     : NULL
#>   ..@ linejoin     : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ legend.key.size      : 'simpleUnit' num 10points
#>   ..- attr(*, "unit")= int 8
#>  $ legend.text          : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : chr "black"
#>   ..@ size         : num 11
#>   ..@ hjust        : NULL
#>   ..@ vjust        : NULL
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : NULL
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ legend.title         : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : chr "black"
#>   ..@ size         : num 12
#>   ..@ hjust        : num 0
#>   ..@ vjust        : NULL
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : NULL
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ panel.background     : <ggplot2::element_rect>
#>   ..@ fill         : chr "white"
#>   ..@ colour       : chr "white"
#>   ..@ linewidth    : NULL
#>   ..@ linetype     : NULL
#>   ..@ linejoin     : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ panel.border         : <ggplot2::element_blank>
#>  $ plot.background      : <ggplot2::element_rect>
#>   ..@ fill         : chr "white"
#>   ..@ colour       : chr "white"
#>   ..@ linewidth    : NULL
#>   ..@ linetype     : NULL
#>   ..@ linejoin     : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ plot.title           : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : chr "black"
#>   ..@ size         : num 14
#>   ..@ hjust        : NULL
#>   ..@ vjust        : num 1
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : NULL
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ plot.subtitle        : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : NULL
#>   ..@ size         : num 13
#>   ..@ hjust        : num 0
#>   ..@ vjust        : NULL
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : <ggplot2::margin> num [1:4] 0 0 3 0
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ strip.background     : <ggplot2::element_rect>
#>   ..@ fill         : chr "transparent"
#>   ..@ colour       : NULL
#>   ..@ linewidth    : NULL
#>   ..@ linetype     : num 0
#>   ..@ linejoin     : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ strip.placement      : chr "outside"
#>  $ strip.text           : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : chr "black"
#>   ..@ size         : num 12.5
#>   ..@ hjust        : num 0.5
#>   ..@ vjust        : NULL
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : <ggplot2::margin> num [1:4] 3 3 3 3
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ strip.switch.pad.grid: 'simpleUnit' num -1points
#>   ..- attr(*, "unit")= int 8
#>  $ strip.switch.pad.wrap: 'simpleUnit' num -1points
#>   ..- attr(*, "unit")= int 8
#>  $ axis.ticks           : <ggplot2::element_blank>
#>  $ panel.grid           : <ggplot2::element_blank>
#>  @ complete: logi FALSE
#>  @ validate: logi TRUE
```
