# Prefetch cell cycle genes

Based on the human cell cycle genes, the cell cycle genes of the
corresponding species were captured by homologous gene conversion.

## Usage

``` r
CycGenePrefetch(
  species = "Homo_sapiens",
  Ensembl_version = NULL,
  mirror = NULL,
  max_tries = 5,
  use_cached_gene = TRUE,
  verbose = TRUE
)
```

## Arguments

- species:

  Latin names for animals, i.e., `"Homo_sapiens"`, `"Mus_musculus"`

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- use_cached_gene:

  Whether to use previously cached cell cycle gene conversion results
  for the species. Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list of S-phase and G2M-phase genes.

## See also

[GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)

## Examples

``` r
ccgenes <- CycGenePrefetch("Homo_sapiens")
#> ℹ [2026-05-24 15:01:44] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-05-24 15:01:44] Cell cycle gene prefetching completed "Homo_sapiens"
str(ccgenes)
#> List of 5
#>  $ res         : NULL
#>  $ S           : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M         : chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
#>  $ cc_S_genes  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ cc_G2M_genes: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...

ccgenes <- CycGenePrefetch("Mus_musculus")
#> ℹ [2026-05-24 15:01:44] Prefetching cell cycle genes for "Mus_musculus" ...
#> ℹ [2026-05-24 15:01:44] Connect to the Ensembl archives...
#> ℹ [2026-05-24 15:01:45] Using the 115 version of ensembl database...
#> ℹ [2026-05-24 15:01:45] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-24 15:02:06] Searching the dataset hsapiens ...
#> ! [2026-05-24 15:02:06] No matching datasets found
#> ! [2026-05-24 15:02:06] Can not find the dataset for the species: Homo_sapiens (hsapiens)
#> ℹ [2026-05-24 15:02:06] Cached conversion results for "Mus_musculus"
#> ✔ [2026-05-24 15:02:06] Cell cycle gene prefetching completed "Mus_musculus"
str(ccgenes)
#> List of 5
#>  $ res         :List of 6
#>   ..$ geneID_res     : NULL
#>   ..$ geneID_collapse: NULL
#>   ..$ geneID_expand  : NULL
#>   ..$ Ensembl_version: chr "115"
#>   ..$ Datasets       :'data.frame':  0 obs. of  3 variables:
#>   .. ..$ dataset    : 'AsIs' chr(0) 
#>   .. ..$ description: 'AsIs' chr(0) 
#>   .. ..$ version    : 'AsIs' chr(0) 
#>   ..$ Attributes     : NULL
#>  $ S           : NULL
#>  $ G2M         : NULL
#>  $ cc_S_genes  : NULL
#>  $ cc_G2M_genes: NULL
```
