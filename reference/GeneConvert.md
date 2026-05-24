# Gene ID conversion function using biomart

This function can convert different gene ID types within one species or
between two species using the biomart service.

## Usage

``` r
GeneConvert(
  geneID,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "entrez_id",
  species_from = "Homo_sapiens",
  species_to = NULL,
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)
```

## Arguments

- geneID:

  A vector of the geneID character.

- geneID_from_IDtype:

  Gene ID type of the input `geneID`. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`

- geneID_to_IDtype:

  Gene ID type(s) to convert to. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`.

- species_from:

  Latin names for animals of the input geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- species_to:

  Latin names for animals of the output geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- biomart:

  The name of the BioMart database that you want to connect to. Possible
  options include `"ensembl"`, `"protists_mart"`, `"fungi_mart"`, and
  `"plants_mart"`.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list with the following elements:

- `geneID_res`: A data.frame contains the all gene IDs mapped in the
  database with columns: `"from_IDtype"`, `"from_geneID"`,
  `"to_IDtype"`, `"to_geneID"`.

- `geneID_collapse`: The data.frame contains all the successfully
  converted gene IDs, and the output gene IDs are collapsed into a list.
  As a result, the `"from_geneID"` column (which is set as the row
  names) of the data.frame is unique.

- `geneID_expand`: The data.frame contains all the successfully
  converted gene IDs, and the output gene IDs are expanded.

- `Ensembl_version`: Ensembl database version.

- `Datasets`: Datasets available in the selected BioMart database.

- `Attributes`: Attributes available in the selected BioMart database.

- `geneID_unmapped`: A character vector of gene IDs that are unmapped in
  the database.

## See also

[AnnotateFeatures](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md),
[ConvertHomologs](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)

## Examples

``` r
res <- GeneConvert(
  geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"),
  species_from = "Homo_sapiens",
  species_to = "Mus_musculus"
)
#> ℹ [2026-05-24 15:28:34] Connect to the Ensembl archives...
#> ℹ [2026-05-24 15:28:34] Using the 115 version of ensembl database...
#> ℹ [2026-05-24 15:28:34] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-24 15:28:42] Searching the dataset hsapiens ...
#> ℹ [2026-05-24 15:29:05] Connecting to the dataset hsapiens_gene_ensembl ...
#> ℹ [2026-05-24 15:29:09] Searching the dataset mmusculus ...
#> ℹ [2026-05-24 15:29:10] Connecting to the dataset mmusculus_gene_ensembl ...
#> ℹ [2026-05-24 15:29:10] Converting the geneIDs...
#> ℹ [2026-05-24 15:29:15] 5 genes mapped with "ensembl_symbol"
#> ℹ [2026-05-24 15:29:15] ==============================
#> ℹ                       5 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
str(res)
#> List of 7
#>  $ geneID_res     :'data.frame': 5 obs. of  4 variables:
#>   ..$ from_IDtype: chr [1:5] "ensembl_symbol" "ensembl_symbol" "ensembl_symbol" "ensembl_symbol" ...
#>   ..$ from_geneID: chr [1:5] "CTCF" "CDK1" "TOP2A" "AURKA" ...
#>   ..$ to_IDtype  : chr [1:5] "entrez_id" "entrez_id" "entrez_id" "entrez_id" ...
#>   ..$ to_geneID  : int [1:5] 13018 12534 21973 20878 17345
#>  $ geneID_collapse:'data.frame': 5 obs. of  2 variables:
#>   ..$ from_geneID: chr [1:5] "AURKA" "CDK1" "CTCF" "MKI67" ...
#>   ..$ entrez_id  :List of 5
#>   .. ..$ : int 20878
#>   .. ..$ : int 12534
#>   .. ..$ : int 13018
#>   .. ..$ : int 17345
#>   .. ..$ : int 21973
#>   .. ..- attr(*, "class")= chr "AsIs"
#>  $ geneID_expand  :'data.frame': 5 obs. of  2 variables:
#>   ..$ from_geneID: chr [1:5] "AURKA" "CDK1" "CTCF" "MKI67" ...
#>   ..$ entrez_id  : int [1:5] 20878 12534 13018 17345 21973
#>  $ Ensembl_version: chr "115"
#>  $ Datasets       :'data.frame': 213 obs. of  3 variables:
#>   ..$ dataset    : 'AsIs' chr [1:213] "abrachyrhynchus_gene_ensembl" "acalliptera_gene_ensembl" "acarolinensis_gene_ensembl" "acchrysaetos_gene_ensembl" ...
#>   ..$ description: 'AsIs' chr [1:213] "Pink-footed goose genes (ASM259213v1)" "Eastern happy genes (fAstCal1.3)" "Green anole genes (AnoCar2.0v2)" "Golden eagle genes (bAquChr1.2)" ...
#>   ..$ version    : 'AsIs' chr [1:213] "ASM259213v1" "fAstCal1.3" "AnoCar2.0v2" "bAquChr1.2" ...
#>  $ Attributes     :'data.frame': 3170 obs. of  3 variables:
#>   ..$ name       : chr [1:3170] "ensembl_gene_id" "ensembl_gene_id_version" "ensembl_transcript_id" "ensembl_transcript_id_version" ...
#>   ..$ description: chr [1:3170] "Gene stable ID" "Gene stable ID version" "Transcript stable ID" "Transcript stable ID version" ...
#>   ..$ page       : chr [1:3170] "feature_page" "feature_page" "feature_page" "feature_page" ...
#>  $ geneID_unmapped: chr(0) 
```
