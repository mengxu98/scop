# Compact Xenium human pancreas cell boundaries

Real cell-segmentation polygon vertices extracted from the
TENxXeniumData human pancreas SpatialExperiment. The plain data frame is
directly suitable for SpatialCellPlot and keeps no Bioconductor geometry
dependency at runtime.

## Format

A data frame with 3900 polygon vertices for 300 cells.

## Source

Derived from TENxXeniumData::spe_human_pancreas and preserved by the
reproducible assets in
<https://github.com/mengxu98/datasets/tree/main/SCOPExamples>.

## Examples

``` r
data(xenium_human_pancreas_boundaries_sub)
c(
  cells = length(unique(xenium_human_pancreas_boundaries_sub$cell_id)),
  vertices = nrow(xenium_human_pancreas_boundaries_sub)
)
#>    cells vertices 
#>      300     3900 
head(xenium_human_pancreas_boundaries_sub)
#>      cell_id polygon_id ring_id vertex_order        x        y           image
#> 1 gncmhnff-1  polygon_1  ring_1            1 4286.337 1156.638 xenium_pancreas
#> 2 gncmhnff-1  polygon_1  ring_1            2 4283.363 1160.037 xenium_pancreas
#> 3 gncmhnff-1  polygon_1  ring_1            3 4283.150 1161.738 xenium_pancreas
#> 4 gncmhnff-1  polygon_1  ring_1            4 4285.062 1165.350 xenium_pancreas
#> 5 gncmhnff-1  polygon_1  ring_1            5 4286.550 1167.688 xenium_pancreas
#> 6 gncmhnff-1  polygon_1  ring_1            6 4287.613 1167.900 xenium_pancreas
#>   cell_area
#> 1  96.54501
#> 2  96.54501
#> 3  96.54501
#> 4  96.54501
#> 5  96.54501
#> 6  96.54501
```
