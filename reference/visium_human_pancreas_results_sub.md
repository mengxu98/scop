# Compact Visium result example for spatial documentation

A spatially covering subset of the bundled human PanIN Visium object. It
retains raw x/y coordinates and native SCOP network and neighborhood
results for short documentation examples. Optional backend results are
only included when generated and validated from the pinned source
workflow.

## Format

A Seurat object with 5000 genes and 400 spots.

## Source

Derived from GSE254829 GSM8058244 through the reproducible assets in
<https://github.com/mengxu98/datasets/tree/main/SCOPExamples>.

## Examples

``` r
data(visium_human_pancreas_results_sub)
c(
  cells = ncol(visium_human_pancreas_results_sub),
  features = nrow(visium_human_pancreas_results_sub)
)
#>    cells features 
#>      400     5000 
names(visium_human_pancreas_results_sub@tools)
#>  [1] "GSE254829_coda_table"    "SpatialNetwork"         
#>  [3] "SpatialNeighborhood"     "SpatialGradientFeatures"
#>  [5] "SCOPExamples"            "MistyR"                 
#>  [7] "StatialKontextual"       "STdeconvolve"           
#>  [9] "SpatialCellChat"         "CCC"                    
#> [11] "MERINGUE"               
```
