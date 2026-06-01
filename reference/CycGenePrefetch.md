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
#> ℹ [2026-06-01 08:56:31] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-06-01 08:56:31] Cell cycle gene prefetching completed "Homo_sapiens"
str(ccgenes)
#> List of 5
#>  $ res         : NULL
#>  $ S           : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M         : chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
#>  $ cc_S_genes  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ cc_G2M_genes: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...

ccgenes <- CycGenePrefetch("Mus_musculus")
#> ℹ [2026-06-01 08:56:31] Prefetching cell cycle genes for "Mus_musculus" ...
#> ℹ [2026-06-01 08:56:31] Connect to the Ensembl archives...
#> ℹ [2026-06-01 08:56:31] Using the 115 version of ensembl database...
#> ℹ [2026-06-01 08:56:31] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> ! [2026-06-01 08:56:34] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-06-01 08:56:34] Get errors when connecting with ensembl database...
#> ! [2026-06-01 08:56:35] Retrying...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-06-01 08:56:38] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-06-01 08:56:38] Get errors when connecting with ensembl database...
#> ! [2026-06-01 08:56:39] Retrying...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-06-01 08:56:42] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-06-01 08:56:42] Get errors when connecting with ensembl database...
#> ! [2026-06-01 08:56:43] Retrying...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> ! [2026-06-01 08:56:46] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-06-01 08:56:46] Get errors when connecting with ensembl database...
#> ! [2026-06-01 08:56:47] Retrying...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-06-01 08:56:50] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-06-01 08:56:50] Get errors when connecting with ensembl database...
#> Error in try_get(expr = {    if (!is.null(mirror)) {        biomaRt::useEnsembl(biomart = "ensembl", mirror = mirror)    }    else {        mart_try <- tryCatch(biomaRt::useMart(biomart = "ensembl",             host = url), error = function(e) e)        if (inherits(mart_try, "error")) {            mirror_candidates <- c("useast", "uswest", "asia",                 "www")            for (mirror_i in mirror_candidates) {                mart_mirror <- tryCatch(biomaRt::useEnsembl(biomart = "ensembl",                   mirror = mirror_i), error = function(e) e)                if (!inherits(mart_mirror, "error")) {                  mart_try <- mart_mirror                  break                }            }        }        if (inherits(mart_try, "error")) {            stop(mart_try)        }        mart_try    }}, max_tries = max_tries, error_message = "Get errors when connecting with ensembl database..."): <simpleError: Your query has been redirected to
#> https://status.ensembl.org indicating this Ensembl service is currently
#> unavailable. Look at ?useEnsembl for details on how to try a mirror site.>
str(ccgenes)
#> List of 5
#>  $ res         : NULL
#>  $ S           : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M         : chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
#>  $ cc_S_genes  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ cc_G2M_genes: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
```
