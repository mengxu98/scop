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
#> ℹ [2026-01-29 12:38:27] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-01-29 12:38:27] Cell cycle gene prefetching completed "Homo_sapiens"
str(ccgenes)
#> List of 3
#>  $ res: NULL
#>  $ S  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...

ccgenes <- CycGenePrefetch("Mus_musculus")
#> ℹ [2026-01-29 12:38:27] Prefetching cell cycle genes for "Mus_musculus" ...
#> ℹ [2026-01-29 12:38:48] Connect to the Ensembl archives...
#> ℹ [2026-01-29 12:38:49] Using the 115 version of ensembl database...
#> ℹ [2026-01-29 12:38:49] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-01-29 12:38:51] Searching the dataset hsapiens ...
#> ℹ [2026-01-29 12:38:51] Connecting to the dataset hsapiens_gene_ensembl ...
#> ℹ [2026-01-29 12:38:52] Converting the geneIDs...
#> ℹ [2026-01-29 12:38:54] 97 genes mapped with "ensembl_symbol"
#> ℹ [2026-01-29 12:38:54] ==============================
#> ℹ                       97 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-01-29 12:38:54] Cached conversion results for "Mus_musculus"
#> ✔ [2026-01-29 12:38:54] Cell cycle gene prefetching completed "Mus_musculus"
str(ccgenes)
#> List of 3
#>  $ res:List of 7
#>   ..$ geneID_res     :'data.frame':  100 obs. of  4 variables:
#>   .. ..$ from_IDtype: chr [1:100] "ensembl_symbol" "ensembl_symbol" "ensembl_symbol" "ensembl_symbol" ...
#>   .. ..$ from_geneID: chr [1:100] "NCAPD2" "ANLN" "UBR7" "TACC3" ...
#>   .. ..$ to_IDtype  : chr [1:100] "symbol" "symbol" "symbol" "symbol" ...
#>   .. ..$ to_geneID  : chr [1:100] "Ncapd2" "Anln" "Ubr7" "Tacc3" ...
#>   ..$ geneID_collapse:'data.frame':  96 obs. of  2 variables:
#>   .. ..$ from_geneID: chr [1:96] "ANLN" "ANP32E" "ATAD2" "AURKA" ...
#>   .. ..$ symbol     :List of 96
#>   .. .. ..$ : chr "Anln"
#>   .. .. ..$ : chr "Anp32e"
#>   .. .. ..$ : chr "Atad2"
#>   .. .. ..$ : chr "Aurka"
#>   .. .. ..$ : chr "Aurkb"
#>   .. .. ..$ : chr "Birc5"
#>   .. .. ..$ : chr "Blm"
#>   .. .. ..$ : chr "Bub1"
#>   .. .. ..$ : chr "Casp8ap2"
#>   .. .. ..$ : chr "Cbx5"
#>   .. .. ..$ : chr "Ccnb2"
#>   .. .. ..$ : chr "Ccne2"
#>   .. .. ..$ : chr "Cdc20"
#>   .. .. ..$ : chr "Cdc25c"
#>   .. .. ..$ : chr "Cdc45"
#>   .. .. ..$ : chr "Cdc6"
#>   .. .. ..$ : chr "Cdca2"
#>   .. .. ..$ : chr "Cdca3"
#>   .. .. ..$ : chr "Cdca7"
#>   .. .. ..$ : chr "Cdca8"
#>   .. .. ..$ : chr "Cdk1"
#>   .. .. ..$ : chr "Cenpa"
#>   .. .. ..$ : chr "Cenpe"
#>   .. .. ..$ : chr "Cenpf"
#>   .. .. ..$ : chr "Cenpu"
#>   .. .. ..$ : chr "Chaf1b"
#>   .. .. ..$ : chr "Ckap2"
#>   .. .. ..$ : chr "Ckap2l"
#>   .. .. ..$ : chr "Ckap5"
#>   .. .. ..$ : chr [1:2] "Cks1brt" "Cks1b"
#>   .. .. ..$ : chr "Cks2"
#>   .. .. ..$ : chr "Clspn"
#>   .. .. ..$ : chr "Ctcf"
#>   .. .. ..$ : chr "Dlgap5"
#>   .. .. ..$ : chr "Dscc1"
#>   .. .. ..$ : chr "Dtl"
#>   .. .. ..$ : chr "E2f8"
#>   .. .. ..$ : chr "Ect2"
#>   .. .. ..$ : chr "Exo1"
#>   .. .. ..$ : chr "Fen1"
#>   .. .. ..$ : chr "G2e3"
#>   .. .. ..$ : chr "Gas2l3"
#>   .. .. ..$ : chr "Gins2"
#>   .. .. ..$ : chr "Gmnn"
#>   .. .. ..$ : chr "Gtse1"
#>   .. .. ..$ : chr "Hells"
#>   .. .. ..$ : chr "Hjurp"
#>   .. .. ..$ : chr "Hmgb2"
#>   .. .. ..$ : chr "Hmmr"
#>   .. .. ..$ : chr "Jpt1"
#>   .. .. ..$ : chr "Kif11"
#>   .. .. ..$ : chr "Kif20b"
#>   .. .. ..$ : chr "Kif23"
#>   .. .. ..$ : chr "Kif2c"
#>   .. .. ..$ : chr "Lbr"
#>   .. .. ..$ : chr "Mcm4"
#>   .. .. ..$ : chr "Mcm5"
#>   .. .. ..$ : chr "Mcm6"
#>   .. .. ..$ : chr "Mcm7"
#>   .. .. ..$ : chr "Mki67"
#>   .. .. ..$ : chr "Mrpl36"
#>   .. .. ..$ : chr "Msh2"
#>   .. .. ..$ : chr "Nasp"
#>   .. .. ..$ : chr "Ncapd2"
#>   .. .. ..$ : chr "Ndc80"
#>   .. .. ..$ : chr "Nek2"
#>   .. .. ..$ : chr "Nuf2"
#>   .. .. ..$ : chr "Nusap1"
#>   .. .. ..$ : chr "Pcna"
#>   .. .. ..$ : chr "Pimreg"
#>   .. .. ..$ : chr "Pola1"
#>   .. .. ..$ : chr "Polr1b"
#>   .. .. ..$ : chr "Prim1"
#>   .. .. ..$ : chr "Psrc1"
#>   .. .. ..$ : chr "Rad51"
#>   .. .. ..$ : chr "Rad51ap1"
#>   .. .. ..$ : chr "Rangap1"
#>   .. .. ..$ : chr "Rfc2"
#>   .. .. ..$ : chr "Rrm1"
#>   .. .. ..$ : chr "Rrm2"
#>   .. .. ..$ : chr "Slbp"
#>   .. .. ..$ : chr "Smc4"
#>   .. .. ..$ : chr "Tacc3"
#>   .. .. ..$ : chr "Tipin"
#>   .. .. ..$ : chr "Tmpo"
#>   .. .. ..$ : chr "Top2a"
#>   .. .. ..$ : chr "Tpx2"
#>   .. .. ..$ : chr "Ttk"
#>   .. .. ..$ : chr "Tubb4b"
#>   .. .. ..$ : chr "Tyms"
#>   .. .. ..$ : chr "Ube2c"
#>   .. .. ..$ : chr "Ubr7"
#>   .. .. ..$ : chr "Uhrf1"
#>   .. .. ..$ : chr "Ung"
#>   .. .. ..$ : chr "Usp1"
#>   .. .. ..$ : chr "Wdr76"
#>   .. .. ..- attr(*, "class")= chr "AsIs"
#>   ..$ geneID_expand  :'data.frame':  97 obs. of  2 variables:
#>   .. ..$ from_geneID: chr [1:97] "ANLN" "ANP32E" "ATAD2" "AURKA" ...
#>   .. ..$ symbol     : chr [1:97] "Anln" "Anp32e" "Atad2" "Aurka" ...
#>   ..$ Ensembl_version: chr "115"
#>   ..$ Datasets       :'data.frame':  213 obs. of  3 variables:
#>   .. ..$ dataset    : 'AsIs' chr [1:213] "abrachyrhynchus_gene_ensembl" "acalliptera_gene_ensembl" "acarolinensis_gene_ensembl" "acchrysaetos_gene_ensembl" ...
#>   .. ..$ description: 'AsIs' chr [1:213] "Pink-footed goose genes (ASM259213v1)" "Eastern happy genes (fAstCal1.3)" "Green anole genes (AnoCar2.0v2)" "Golden eagle genes (bAquChr1.2)" ...
#>   .. ..$ version    : 'AsIs' chr [1:213] "ASM259213v1" "fAstCal1.3" "AnoCar2.0v2" "bAquChr1.2" ...
#>   ..$ Attributes     :'data.frame':  3170 obs. of  3 variables:
#>   .. ..$ name       : chr [1:3170] "ensembl_gene_id" "ensembl_gene_id_version" "ensembl_transcript_id" "ensembl_transcript_id_version" ...
#>   .. ..$ description: chr [1:3170] "Gene stable ID" "Gene stable ID version" "Transcript stable ID" "Transcript stable ID version" ...
#>   .. ..$ page       : chr [1:3170] "feature_page" "feature_page" "feature_page" "feature_page" ...
#>   ..$ geneID_unmapped: chr(0) 
#>  $ S  : chr [1:42] "Mcm5" "Pcna" "Tyms" "Fen1" ...
#>  $ G2M: chr [1:55] "Hmgb2" "Cdk1" "Nusap1" "Ube2c" ...
```
