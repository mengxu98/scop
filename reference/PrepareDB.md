# Prepare gene annotation databases

Build TERM2GENE / TERM2NAME (and GO semantic-similarity) databases for a
species from annotation packages, cached downloads, or local files.

## Usage

``` r
PrepareDB(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = c("GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM",
    "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink",
    "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat",
    "Chromosome", "GeneType", "Enzyme", "TF", "CytoTRACE2"),
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest",
  db_update = FALSE,
  data_dir = NULL,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  biomart = NULL,
  max_tries = 5,
  custom_TERM2GENE = NULL,
  custom_TERM2NAME = NULL,
  custom_species = NULL,
  custom_IDtype = NULL,
  custom_version = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- species:

  `"Homo_sapiens"` or `"Mus_musculus"`.

- db:

  Annotation sources. One or more of `"GO"`, `"GO_BP"`, `"GO_CC"`,
  `"GO_MF"`, `"KEGG"`, `"WikiPathway"`, `"Reactome"`, `"CORUM"`, `"MP"`,
  `"DO"`, `"HPO"`, `"PFAM"`, `"CSPA"`, `"Surfaceome"`, `"SPRomeDB"`,
  `"VerSeDa"`, `"TFLink"`, `"hTFtarget"`, `"TRRUST"`, `"JASPAR"`,
  `"ENCODE"`, `"MSigDB"`, `"CellTalk"`, `"CellChat"`, `"Chromosome"`,
  `"GeneType"`, `"Enzyme"`, `"TF"`, `"CytoTRACE2"`. MSigDB
  subcollections use `"MSigDB_<collection>"` (e.g. `"MSigDB_H"`).
  `"CytoTRACE2"` is species-independent and is required by
  [RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md).

- db_IDtypes:

  Gene ID types to include.

- db_version:

  Database version to retrieve.

- db_update:

  Force a refresh. `FALSE` loads the cache when available.

- data_dir:

  Directory or named list of local source files. Searches
  `data_dir/<db>/` then `data_dir`. Named lists override a path, e.g.
  `list(MSigDB = "~/db/msigdb")`.

- convert_species:

  Use a species-converted database when the annotation is missing for
  `species`.

- Ensembl_version:

  Ensembl version. `NULL` uses the latest.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- biomart:

  BioMart database that you want to connect to. Possible options include
  `"ensembl"`, `"protists_mart"`, `"fungi_mart"`, and `"plants_mart"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- custom_TERM2GENE, custom_TERM2NAME:

  Custom mappings for `custom_species`.

- custom_species, custom_IDtype, custom_version:

  Metadata for a custom database.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to helper functions.

## Value

A list containing the prepared gene annotation databases:

- `TERM2GENE`: mapping of gene identifiers to terms.

- `TERM2NAME`: mapping of terms to their names.

- `semData`: semantic similarity data for gene sets (only for Gene
  Ontology terms).

## See also

[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
db_list <- PrepareDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
#> ℹ [2026-09-06 21:44:35] Species: "Homo_sapiens"
#> 
#> ℹ [2026-09-06 21:48:54] Preparing database: GO_BP
#> ℹ [2026-09-06 21:49:02] Convert ID types for the GO_BP database
#> ℹ [2026-09-06 21:49:03] Converted ID types using local annotation package org.Hs.eg.db
ListDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
#>   Database      Species            Version                       Date
#> 1    GO_BP Homo_sapiens 3.23.1 nterm:14209 2026-09-06 21:49:04.743761
head(
  db_list[["Homo_sapiens"]][["GO_BP"]][["TERM2GENE"]]
)
#>         Term entrez_id symbol      ensembl_id
#> 1 GO:0000012      2074  ERCC6 ENSG00000225830
#> 2 GO:0000012      7515  XRCC1 ENSG00000073050
#> 3 GO:0000012       142  PARP1 ENSG00000143799
#> 4 GO:0000012      1161  ERCC8 ENSG00000049167
#> 5 GO:0000012     11284   PNKP ENSG00000039650
#> 6 GO:0000012     55775   TDP1 ENSG00000042088

# Based on homologous gene conversion,
# prepare a gene annotation database that originally does not exist in the species.
db_list <- PrepareDB(
  species = "Homo_sapiens",
  db = "MP"
)
#> ℹ [2026-09-06 21:49:05] Species: "Homo_sapiens"
#> ! [2026-09-06 21:49:05] Use the mouse annotation to create the MP database for "Homo_sapiens"
#> ℹ [2026-09-06 21:49:05] Preparing MP database
#> ℹ [2026-09-06 21:49:18] Convert species for the MP database
#> ℹ [2026-09-06 21:49:18] Connect to the Ensembl archives...
#> ℹ [2026-09-06 21:49:18] Using the 116 version of ensembl database...
#> ℹ [2026-09-06 21:49:18] Downloading the ensembl database from https://jun2026.archive.ensembl.org...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> ℹ [2026-09-06 21:49:22] Searching the dataset mmusculus ...
#> ℹ [2026-09-06 21:49:22] Connecting to the dataset mmusculus_gene_ensembl ...
#> ℹ [2026-09-06 21:49:24] Converting the geneIDs...
#> ℹ [2026-09-06 21:49:32] 14164 genes mapped with "ensembl_symbol"
#> ℹ [2026-09-06 21:49:34] 3 genes mapped with "entrez_symbol"
#> ℹ [2026-09-06 21:49:36] 12 genes mapped with "uniprot_symbol"
#> ℹ [2026-09-06 21:49:38] ==============================
#> ℹ                       14179 genes mapped
#> ℹ                       48 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-09-06 21:49:44] Convert ID types for the MP database
#> ℹ [2026-09-06 21:49:44] Converted ID types using local annotation package org.Hs.eg.db
ListDB(
  species = "Homo_sapiens",
  db = "MP"
)
#>   Database      Species                                             Version
#> 1       MP Homo_sapiens 2026-09-06(converted from Mus_musculus) nterm:10803
#>                         Date
#> 1 2026-09-06 21:49:46.990192
head(
  db_list[["Homo_sapiens"]][["MP"]][["TERM2GENE"]]
)
#>         Term      ensembl_id symbol entrez_id
#> 1 MP:0000600 ENSG00000139687    RB1      5925
#> 2 MP:0001716 ENSG00000139687    RB1      5925
#> 3 MP:0001698 ENSG00000139687    RB1      5925
#> 4 MP:0001092 ENSG00000139687    RB1      5925
#> 5 MP:0000961 ENSG00000139687    RB1      5925
#> 6 MP:0000828 ENSG00000139687    RB1      5925

# You can also build a custom database based on the gene sets you have
ccgenes <- CycGenePrefetch("Homo_sapiens")
#> ℹ [2026-09-06 21:49:47] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-09-06 21:49:47] Cell cycle gene prefetching completed "Homo_sapiens"
custom_TERM2GENE <- rbind(
  data.frame(
    term = "S_genes",
    gene = ccgenes[["cc_S_genes"]]
  ),
  data.frame(
    term = "G2M_genes",
    gene = ccgenes[["cc_G2M_genes"]]
  )
)
str(custom_TERM2GENE)
#> 'data.frame':    97 obs. of  2 variables:
#>  $ term: chr  "S_genes" "S_genes" "S_genes" "S_genes" ...
#>  $ gene: chr  "MCM5" "PCNA" "TYMS" "FEN1" ...

# Set convert_species = TRUE to build a custom database for both species,
# with the name "CellCycle"
db_list <- PrepareDB(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "CellCycle",
  convert_species = TRUE,
  custom_TERM2GENE = custom_TERM2GENE,
  custom_species = "Homo_sapiens",
  custom_IDtype = "symbol",
  custom_version = "Seurat_v5"
)
#> ℹ [2026-09-06 21:49:47] Species: "Homo_sapiens"
#> ℹ [2026-09-06 21:49:47] Convert ID types for the CellCycle database
#> ℹ [2026-09-06 21:49:47] Converted ID types using local annotation package org.Hs.eg.db
#> ℹ [2026-09-06 21:49:47] Species: "Mus_musculus"
#> ! [2026-09-06 21:49:47] Use the "Homo_sapiens" annotation to create the "CellCycle" database for "Mus_musculus"
#> ℹ [2026-09-06 21:49:47] Convert species for the CellCycle database
#> ℹ [2026-09-06 21:49:47] Connect to the Ensembl archives...
#> ℹ [2026-09-06 21:49:47] Using the 116 version of ensembl database...
#> ℹ [2026-09-06 21:49:47] Downloading the ensembl database from https://jun2026.archive.ensembl.org...
#> Ensembl site unresponsive, trying asia mirror
#> ℹ [2026-09-06 21:49:48] Searching the dataset hsapiens ...
#> ℹ [2026-09-06 21:49:49] Connecting to the dataset hsapiens_gene_ensembl ...
#> ℹ [2026-09-06 21:49:51] Converting the geneIDs...
#> ℹ [2026-09-06 21:49:53] 97 genes mapped with "ensembl_symbol"
#> ℹ [2026-09-06 21:49:53] ==============================
#> ℹ                       97 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-09-06 21:49:53] Convert ID types for the CellCycle database
#> ℹ [2026-09-06 21:49:53] Converted ID types using local annotation package org.Mm.eg.db
ListDB(db = "CellCycle")
#>    Database      Species                                        Version
#> 1 CellCycle Homo_sapiens                              Seurat_v5 nterm:2
#> 2 CellCycle Mus_musculus Seurat_v5(converted from Homo_sapiens) nterm:2
#>                         Date
#> 1 2026-09-06 21:49:47.371647
#> 2 2026-09-06 21:49:53.636954

db_list <- PrepareDB(species = "Mus_musculus", db = "CellCycle")
#> ℹ [2026-09-06 21:49:53] Species: "Mus_musculus"
#> ℹ [2026-09-06 21:49:53] Loading cached: CellCycle version: Seurat_v5(converted from Homo_sapiens) nterm:2 created: 2026-09-06 21:49:53
head(
  db_list[["Mus_musculus"]][["CellCycle"]][["TERM2GENE"]]
)
#>      Term         ensembl_id symbol entrez_id
#> 1 S_genes ENSMUSG00000005410   Mcm5     17218
#> 2 S_genes ENSMUSG00000027342   Pcna     18538
#> 3 S_genes ENSMUSG00000025747   Tyms     22171
#> 4 S_genes ENSMUSG00000024742   Fen1     14156
#> 5 S_genes ENSMUSG00000029730   Mcm7     17220
#> 6 S_genes ENSMUSG00000022673   Mcm4     17217
```
