# A human pancreas Visium spatial example dataset

A compact gene-filtered version of a human pancreatic intraepithelial
neoplasia (PanIN) 10x Visium dataset from GSE254829. The object keeps
the 1986 non-background tissue spots from sample GSM8058244 (PanIN-LG2),
with a `Spatial` assay, a `slice1` Visium image, and tissue coordinates
in metadata columns `x` and `y`. Metadata column `coda_label` stores the
dominant CODA microanatomical component for each spot, and `coda_score`
stores its percentage. Component percentage columns are stored with the
`coda_` prefix, and the matched CODA table is stored in
`@tools$GSE254829_coda_table`. To keep the package data small and
directly usable with the bundled `panc8_sub` reference, the object
retains the top 5000 genes shared with `panc8_sub`, ranked by total
spatial counts.

## Format

A `Seurat` object with 5000 genes, 1986 spots, and one Visium image
named `slice1`.

## Source

Derived from the GSE254829 human PanIN 10x Visium dataset:
[GSE254829](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254829).
The package object uses the GEO supplementary files
`GSM8058244_PanIN-LG2.tar.gz` and
`GSE254829_codatable_may202024.csv.gz`.

## Examples

``` r
data(visium_human_pancreas_sub)
SeuratObject::Images(visium_human_pancreas_sub)
#> [1] "slice1"
head(visium_human_pancreas_sub@meta.data[, c("x", "y")])
#>                       x   y
#> TGGTATCGGTCTGTAT-1 3859 662
#> ATTATCTCGACAGATC-1 3914 662
#> TGAGATCAAATACTCA-1 3969 662
#> CTGGTCCTAACTTGGC-1 4024 662
#> ATAGTCTTTGACGTGC-1 4079 662
#> GGGTGGTCCAGCCTGT-1 4134 662
SpatialSpotPlot(visium_human_pancreas_sub, group.by = "coda_label")

```
