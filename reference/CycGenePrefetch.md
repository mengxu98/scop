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

  Ensembl database version. If NULL, use the current release version.

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
#> ℹ [2025-11-13 11:42:52] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2025-11-13 11:42:52] Cell cycle gene prefetching completed "Homo_sapiens"
str(ccgenes)
#> List of 3
#>  $ res: NULL
#>  $ S  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...

ccgenes <- CycGenePrefetch("Mus_musculus")
#> ℹ [2025-11-13 11:42:52] Prefetching cell cycle genes for "Mus_musculus" ...
#> ◌ [2025-11-13 11:42:52] Installing: biomaRt...
#>  
#> → Will install 19 packages.
#> → All 19 packages (0 B) are cached.
#> + AnnotationDbi   1.72.0  [bld]
#> + BiocFileCache   3.0.0   [bld]
#> + Biostrings      2.78.0  [bld][cmp]
#> + DBI             1.2.3   
#> + KEGGREST        1.50.0  [bld]
#> + RSQLite         2.4.4   
#> + Seqinfo         1.0.0   [bld]
#> + XVector         0.50.0  [bld][cmp]
#> + biomaRt         2.66.0  [bld]
#> + bit             4.6.0   
#> + bit64           4.6.0-1 
#> + blob            1.2.4   
#> + dbplyr          2.5.1   
#> + filelock        1.0.3   
#> + hms             1.1.4   
#> + httr2           1.2.1   
#> + prettyunits     1.2.0   
#> + progress        1.2.3   
#> + xml2            1.4.1    + ✔ libxml2-dev
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 19 pkgs are cached
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libxml2-dev libcurl4-openssl-dev libssl-dev libpng-dev libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 23 not upgraded.
#> ℹ Building Seqinfo 1.0.0
#> ℹ Building XVector 0.50.0
#> ✔ Installed bit 4.6.0  (96ms)
#> ✔ Installed bit64 4.6.0-1  (117ms)
#> ✔ Installed blob 1.2.4  (176ms)
#> ✔ Installed DBI 1.2.3  (181ms)
#> ✔ Installed dbplyr 2.5.1  (124ms)
#> ✔ Installed filelock 1.0.3  (109ms)
#> ✔ Installed hms 1.1.4  (109ms)
#> ✔ Installed httr2 1.2.1  (106ms)
#> ✔ Installed prettyunits 1.2.0  (108ms)
#> ✔ Installed progress 1.2.3  (104ms)
#> ✔ Installed RSQLite 2.4.4  (170ms)
#> ℹ Building BiocFileCache 3.0.0
#> ✔ Installed xml2 1.4.1  (110ms)
#> ✔ Built BiocFileCache 3.0.0 (5.3s)
#> ✔ Installed BiocFileCache 3.0.0  (34ms)
#> ✔ Built Seqinfo 1.0.0 (6.8s)
#> ✔ Installed Seqinfo 1.0.0  (1.1s)
#> ✔ Built XVector 0.50.0 (10.9s)
#> ✔ Installed XVector 0.50.0  (1s)
#> ℹ Building Biostrings 2.78.0
#> ✔ Built Biostrings 2.78.0 (18.4s)
#> ✔ Installed Biostrings 2.78.0  (1.1s)
#> ℹ Building KEGGREST 1.50.0
#> ✔ Built KEGGREST 1.50.0 (4.6s)
#> ✔ Installed KEGGREST 1.50.0  (1s)
#> ℹ Building AnnotationDbi 1.72.0
#> ✔ Built AnnotationDbi 1.72.0 (11.8s)
#> ✔ Installed AnnotationDbi 1.72.0  (62ms)
#> ℹ Building biomaRt 2.66.0
#> ✔ Built biomaRt 2.66.0 (7.7s)
#> ✔ Installed biomaRt 2.66.0  (33ms)
#> ✔ 1 pkg + 54 deps: kept 36, added 19 [1m 1s]
#> ✔ [2025-11-13 11:43:53] biomaRt installed successfully
#> ℹ [2025-11-13 11:43:53] Connect to the Ensembl archives...
#> ℹ [2025-11-13 11:43:55] Using the 115 version of ensembl database...
#> ℹ [2025-11-13 11:43:55] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2025-11-13 11:43:56] Searching the dataset hsapiens ...
#> ℹ [2025-11-13 11:43:57] Connecting to the dataset hsapiens_gene_ensembl ...
#> ℹ [2025-11-13 11:43:57] Converting the geneIDs...
#> ℹ [2025-11-13 11:43:59] 97 genes mapped with "ensembl_symbol"
#> ℹ [2025-11-13 11:43:59] ==============================
#> ℹ                       97 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
#> ℹ [2025-11-13 11:43:59] Cached conversion results for "Mus_musculus"
#> ✔ [2025-11-13 11:43:59] Cell cycle gene prefetching completed "Mus_musculus"
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
