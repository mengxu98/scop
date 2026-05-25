# Prepare the gene annotation databases

This function prepares the gene annotation databases for a given species
and set of annotation sources. It retrieves the necessary information
from various annotation packages or external resources and organizes it
into a list. The list contains the annotation data for each specified
annotation source.

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
  verbose = TRUE
)
```

## Arguments

- species:

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Can be one or more of
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF", "CytoTRACE2"`.
  Note: `"CytoTRACE2"` is species-independent and downloads pre-trained
  model data required by
  [RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md).

- db_IDtypes:

  A character vector specifying the desired ID types to be used for gene
  identifiers in the gene annotation databases. Default is
  `c("symbol", "entrez_id", "ensembl_id")`.

- db_version:

  A character vector specifying the version of the gene annotation
  databases to be retrieved. Default is `"latest"`.

- db_update:

  Whether the gene annotation databases should be forcefully updated. If
  set to FALSE, the function will attempt to load the cached databases
  instead. Default is `FALSE`.

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- biomart:

  The name of the BioMart database that you want to connect to. Possible
  options include `"ensembl"`, `"protists_mart"`, `"fungi_mart"`, and
  `"plants_mart"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- custom_TERM2GENE:

  A data frame containing a custom TERM2GENE mapping for the specified
  species and annotation source. Default is `NULL`.

- custom_TERM2NAME:

  A data frame containing a custom TERM2NAME mapping for the specified
  species and annotation source. Default is `NULL`.

- custom_species:

  A character vector specifying the species name to be used in a custom
  database. Default is `NULL`.

- custom_IDtype:

  A character vector specifying the ID type to be used in a custom
  database. Default is `NULL`.

- custom_version:

  A character vector specifying the version to be used in a custom
  database. Default is `NULL`.

- verbose:

  Whether to print the message. Default is `TRUE`.

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
#> ℹ [2026-05-25 10:20:27] Species: "Homo_sapiens"
#> ℹ [2026-05-25 10:20:27] Loading cached: GO_BP version: 3.23.1 nterm:14209 created: 2026-05-25 09:30:03
#> ℹ [2026-05-25 10:20:28] Convert ID types for the GO_BP database
#> ℹ [2026-05-25 10:20:28] Converted ID types using local annotation package org.Hs.eg.db
ListDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.23.1 nterm:14209|Homo_sapiens-GO_BP 1779704429 2026-05-25 10:20:29.395787
#>           db_version            db_name
#> 1 3.23.1 nterm:14209 Homo_sapiens-GO_BP
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/82886eae0621aa67a1839db9a8d85cc5.Rcache
#>        Species    DB
#> 1 Homo_sapiens GO_BP
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
if (interactive()) {
  db_list <- PrepareDB(
    species = "Homo_sapiens",
    db = "MP"
  )
  ListDB(
    species = "Homo_sapiens",
    db = "MP"
  )
  head(
    db_list[["Homo_sapiens"]][["MP"]][["TERM2GENE"]]
  )
}

# You can also build a custom database based on the gene sets you have
ccgenes <- CycGenePrefetch("Homo_sapiens")
#> ℹ [2026-05-25 10:20:29] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-05-25 10:20:29] Cell cycle gene prefetching completed "Homo_sapiens"
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
#> ℹ [2026-05-25 10:20:29] Species: "Homo_sapiens"
#> ℹ [2026-05-25 10:20:30] Convert ID types for the CellCycle database
#> ℹ [2026-05-25 10:20:30] Converted ID types using local annotation package org.Hs.eg.db
#> ℹ [2026-05-25 10:20:30] Species: "Mus_musculus"
#> ! [2026-05-25 10:20:30] Use the "Homo_sapiens" annotation to create the "CellCycle" database for "Mus_musculus"
#> ℹ [2026-05-25 10:20:30] Convert species for the CellCycle database
#> ℹ [2026-05-25 10:20:30] Connect to the Ensembl archives...
#> ℹ [2026-05-25 10:20:30] Using the 115 version of ensembl database...
#> ℹ [2026-05-25 10:20:30] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-25 10:21:02] Searching the dataset hsapiens ...
#> ℹ [2026-05-25 10:21:03] Connecting to the dataset hsapiens_gene_ensembl ...
#> ℹ [2026-05-25 10:21:04] Converting the geneIDs...
#> ℹ [2026-05-25 10:21:07] 97 genes mapped with "ensembl_symbol"
#> ℹ [2026-05-25 10:21:07] ==============================
#> ℹ                       97 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-25 10:21:07] Convert ID types for the CellCycle database
#> ℹ [2026-05-25 10:21:07] Converted ID types using local annotation package org.Mm.eg.db
ListDB(db = "CellCycle")
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#> 2 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                                                 comment
#> 1                              Seurat_v5 nterm:2|Homo_sapiens-CellCycle
#> 2 Seurat_v5(converted from Homo_sapiens) nterm:2|Mus_musculus-CellCycle
#>    timestamp                       date
#> 1 1779704430 2026-05-25 10:20:30.262067
#> 2 1779704467 2026-05-25 10:21:07.336615
#>                                       db_version                db_name
#> 1                              Seurat_v5 nterm:2 Homo_sapiens-CellCycle
#> 2 Seurat_v5(converted from Homo_sapiens) nterm:2 Mus_musculus-CellCycle
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/a6aa81007b9564b5bf1f3fa5dc7997fa.Rcache
#> 2 /home/runner/.cache/R/R.cache/c536960e834e01aaece391bfc36f44cf.Rcache
#>        Species        DB
#> 1 Homo_sapiens CellCycle
#> 2 Mus_musculus CellCycle

db_list <- PrepareDB(species = "Mus_musculus", db = "CellCycle")
#> ℹ [2026-05-25 10:21:07] Species: "Mus_musculus"
#> ℹ [2026-05-25 10:21:07] Loading cached: CellCycle version: Seurat_v5(converted from Homo_sapiens) nterm:2 created: 2026-05-25 10:21:07
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
