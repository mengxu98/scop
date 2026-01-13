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
str(ccgenes)
#> List of 3
#>  $ res: NULL
#>  $ S  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...

ccgenes <- CycGenePrefetch("Mus_musculus")
#> Error in loadNamespace(x): there is no package called ‘R.cache’
str(ccgenes)
#> List of 3
#>  $ res: NULL
#>  $ S  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
```
