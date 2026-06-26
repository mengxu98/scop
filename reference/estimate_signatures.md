# ESTIMATE gene signatures

Gene signatures and common-gene universe used by
[`RunESTIMATE()`](https://mengxu98.github.io/scop/reference/RunESTIMATE.md)
to compute stromal, immune, combined ESTIMATE, and tumor-purity scores
without requiring the external `estimate` package.

## Usage

``` r
data(estimate_signatures)
```

## Format

A `list` with five entries:

- stromal_signature:

  Character vector of stromal signature genes.

- immune_signature:

  Character vector of immune signature genes.

- common_genes:

  Character vector of common genes used for filtering.

- common_gene_aliases:

  Data frame mapping common gene symbols to aliases.

- source:

  List with source method, source data, and DOI metadata.

## Source

Derived from `tidyestimate` 1.1.1 CRAN data files, which are derived
from the MD Anderson ESTIMATE implementation. The scoring references are
[Yoshihara et al. (2013)](https://doi.org/10.1038/ncomms3612) and
[Barbie et al. (2009)](https://doi.org/10.1038/nature08460).

## Examples

``` r
data(estimate_signatures)
names(estimate_signatures)
#> [1] "stromal_signature"   "immune_signature"    "common_genes"       
#> [4] "common_gene_aliases" "source"             
lengths(estimate_signatures[c("stromal_signature", "immune_signature", "common_genes")])
#> stromal_signature  immune_signature      common_genes 
#>               141               141             10391 
```
