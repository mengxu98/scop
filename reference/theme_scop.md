# Default ggplot2 theme for scop plots

Backward-compatible alias of
[`thisplot::theme_this()`](https://mengxu98.github.io/thisplot/reference/theme_this.html).
Plot helpers default to `theme_use = "theme_scop"`. Exporting this name
lets users call `theme_scop()` after
[`library(scop)`](https://mengxu98.github.io/scop/) without `:::`.

## Usage

``` r
theme_scop(aspect.ratio = NULL, base_size = 12, ...)
```

## Arguments

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

[`thisplot::theme_this()`](https://mengxu98.github.io/thisplot/reference/theme_this.html)

## Examples

``` r
theme_scop()
#> <theme> List of 20
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
#>  $ aspect.ratio         : NULL
#>  $ axis.title           : <ggplot2::element_text>
#>   ..@ family       : NULL
#>   ..@ face         : NULL
#>   ..@ italic       : chr NA
#>   ..@ fontweight   : num NA
#>   ..@ fontwidth    : num NA
#>   ..@ colour       : chr "black"
#>   ..@ size         : num 13
#>   ..@ hjust        : NULL
#>   ..@ vjust        : NULL
#>   ..@ angle        : NULL
#>   ..@ lineheight   : NULL
#>   ..@ margin       : NULL
#>   ..@ debug        : NULL
#>   ..@ inherit.blank: logi FALSE
#>  $ axis.text            : <ggplot2::element_text>
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
#>  $ panel.border         : <ggplot2::element_rect>
#>   ..@ fill         : chr "transparent"
#>   ..@ colour       : chr "black"
#>   ..@ linewidth    : num 1
#>   ..@ linetype     : NULL
#>   ..@ linejoin     : NULL
#>   ..@ inherit.blank: logi FALSE
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
#>  @ complete: logi FALSE
#>  @ validate: logi TRUE
```
