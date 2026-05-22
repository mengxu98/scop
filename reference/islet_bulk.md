# Human pancreatic islet bulk RNA-seq example dataset

A full human pancreatic islet bulk RNA-seq `SummarizedExperiment`
derived from a brefeldin A perturbation study. The object keeps all
samples from the published islet arm and stores a symbol-level count
matrix that can be used directly in bulk DE and deconvolution examples
together with the bundled `panc8_sub` reference.

## Format

A `SummarizedExperiment` object with 19876 genes and 8 bulk RNA-seq
samples.

## Source

Derived from
[GSE152615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152615).
The bundled object is built from the supplementary matrix
`GSE152615_Rawcounts_filtered.txt.gz`. The published non-integer count
values are rounded to the nearest integer for count-based example
workflows.

## Examples

``` r
data(islet_bulk)
SummarizedExperiment::assayNames(islet_bulk)
#> [1] "counts"
head(rownames(islet_bulk))
#> [1] "A1BG"     "A1BG-AS1" "A1CF"     "A2M"      "A2M-AS1"  "A2ML1"   
table(SummarizedExperiment::colData(islet_bulk)$condition)
#> 
#>     bfa control 
#>       4       4 
```
