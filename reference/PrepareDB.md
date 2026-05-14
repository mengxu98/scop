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
#> ℹ [2026-05-14 06:36:42] Species: "Homo_sapiens"
#> ℹ [2026-05-14 06:36:42] Loading cached: GO_BP version: 3.23.1 nterm:14209 created: 2026-05-14 05:50:31
#> ℹ [2026-05-14 06:36:43] Convert ID types for the GO_BP database
#> ℹ [2026-05-14 06:36:43] Converted ID types using local annotation package org.Hs.eg.db
ListDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.23.1 nterm:14209|Homo_sapiens-GO_BP 1778740604 2026-05-14 06:36:44.344281
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

# Prepare databases for other species
db_list <- PrepareDB(
  species = "Macaca_fascicularis",
  db = "GO_BP",
  convert_species = TRUE
)
#> ℹ [2026-05-14 06:36:44] Species: "Macaca_fascicularis"
#> ! [2026-05-14 06:36:47] Annotation package org.Mf.eg.db does not exist
#> ! [2026-05-14 06:36:47] Use the human annotation to create the GO_BP database for "Macaca_fascicularis"
#> ✔ [2026-05-14 06:36:47] org.Hs.eg.db installed successfully
#> ℹ [2026-05-14 06:37:10] Preparing database: GO_BP
#> ℹ [2026-05-14 06:37:21] Convert species for the GO_BP database
#> ℹ [2026-05-14 06:37:21] Connect to the Ensembl archives...
#> ℹ [2026-05-14 06:37:21] Using the 115 version of ensembl database...
#> ℹ [2026-05-14 06:37:21] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-14 06:37:24] Searching the dataset hsapiens ...
#> ℹ [2026-05-14 06:37:24] Connecting to the dataset hsapiens_gene_ensembl ...
#> ℹ [2026-05-14 06:37:25] Converting the geneIDs...
#> ℹ [2026-05-14 06:37:44] 18209 genes mapped with "entrez_id"
#> ℹ [2026-05-14 06:37:44] ==============================
#> ℹ                       18209 genes mapped
#> ℹ                       633 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-14 06:38:04] Convert ID types for the GO_BP database
#> ℹ [2026-05-14 06:38:04] Connect to the Ensembl archives...
#> ℹ [2026-05-14 06:38:04] Using the 115 version of ensembl database...
#> ℹ [2026-05-14 06:38:04] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-14 06:38:05] Searching the dataset mfascicularis ...
#> ℹ [2026-05-14 06:38:05] Connecting to the dataset mfascicularis_gene_ensembl ...
#> ℹ [2026-05-14 06:38:05] Converting the geneIDs...
#> ℹ [2026-05-14 06:38:16] 15880 genes mapped with "ensembl_id"
#> ℹ [2026-05-14 06:38:16] ==============================
#> ℹ                       15880 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
ListDB(
  species = "Macaca_fascicularis",
  db = "GO_BP"
)
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                                                     comment
#> 1 3.23.1(converted from Homo_sapiens) nterm:13988|Macaca_fascicularis-GO_BP
#>    timestamp                       date
#> 1 1778740722 2026-05-14 06:38:42.131725
#>                                        db_version                   db_name
#> 1 3.23.1(converted from Homo_sapiens) nterm:13988 Macaca_fascicularis-GO_BP
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/144353694e2c2f86321358fd146ced75.Rcache
#>               Species    DB
#> 1 Macaca_fascicularis GO_BP
head(
  db_list[["Macaca_fascicularis"]][["GO_BP"]][["TERM2GENE"]]
)
#>         Term         ensembl_id symbol entrez_id
#> 1 GO:0000012 ENSMFAG00000035553  ERCC6      <NA>
#> 2 GO:0000012 ENSMFAG00000040166  XRCC1      <NA>
#> 3 GO:0000012 ENSMFAG00000002526  PARP1 101925579
#> 4 GO:0000012 ENSMFAG00000037556  ERCC8 102143286
#> 5 GO:0000012 ENSMFAG00000039622   PNKP 102140873
#> 6 GO:0000012 ENSMFAG00000000067   TDP1 101866378

db_list <- PrepareDB(
  species = "Saccharomyces_cerevisiae",
  db = "GO_BP"
)
#> ℹ [2026-05-14 06:38:42] Species: "Saccharomyces_cerevisiae"
#> 
#> ✔ [2026-05-14 06:39:31] org.Sc.sgd.db installed successfully
#> ℹ [2026-05-14 06:39:41] Preparing database: GO_BP
#> ℹ [2026-05-14 06:39:45] Convert ID types for the GO_BP database
#> ℹ [2026-05-14 06:39:45] Connect to the Ensembl archives...
#> ℹ [2026-05-14 06:39:45] Using the 115 version of ensembl database...
#> ℹ [2026-05-14 06:39:45] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-14 06:39:46] Searching the dataset scerevisiae ...
#> ℹ [2026-05-14 06:39:46] Connecting to the dataset scerevisiae_gene_ensembl ...
#> ℹ [2026-05-14 06:39:46] Converting the geneIDs...
#> ℹ [2026-05-14 06:39:49] 5793 genes mapped with "entrez_id"
#> ℹ [2026-05-14 06:39:49] ==============================
#> ℹ                       5793 genes mapped
#> ℹ                       656 genes unmapped
#> ℹ                       ==============================
ListDB(
  species = "Saccharomyces_cerevisiae",
  db = "GO_BP"
)
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                            comment  timestamp
#> 1 3.22.0 nterm:4474|Saccharomyces_cerevisiae-GO_BP 1778740795
#>                         date        db_version                        db_name
#> 1 2026-05-14 06:39:55.098693 3.22.0 nterm:4474 Saccharomyces_cerevisiae-GO_BP
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/526227780e9714f17fb7fabf51d9fd1b.Rcache
#>                    Species    DB
#> 1 Saccharomyces_cerevisiae GO_BP
head(
  db_list[["Saccharomyces_cerevisiae"]][["GO_BP"]][["TERM2GENE"]]
)
#>         Term entrez_id symbol ensembl_id
#> 1 GO:0000001    854867   MDM1    YML104C
#> 2 GO:0000001    850887   MMR1    YLR190W
#> 3 GO:0000001    855094   ABF2    YMR072W
#> 4 GO:0000001    856198  MDM36    YPR083W
#> 5 GO:0000001    855412  YPT11    YNL304W
#> 6 GO:0000001    851532   ARP2    YDL029W

# Prepare databases for Arabidopsis (plant)
db_list <- PrepareDB(
  species = "Arabidopsis_thaliana",
  db = c(
    "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway",
    "ENZYME", "Chromosome"
  ),
  biomart = "plants_mart"
)
#> ℹ [2026-05-14 06:39:55] Species: "Arabidopsis_thaliana"
#> 
#> ✔ [2026-05-14 06:42:11] org.At.tair.db installed successfully
#> ℹ [2026-05-14 06:42:21] Preparing database: GO_BP
#> ℹ [2026-05-14 06:42:27] Preparing database: GO_CC
#> ℹ [2026-05-14 06:42:28] Preparing database: GO_MF
#> ℹ [2026-05-14 06:42:35] Preparing KEGG database
#> ℹ [2026-05-14 06:42:50] Preparing WikiPathway database
#> ℹ [2026-05-14 06:42:51] Preparing Chromosome database
#> ℹ [2026-05-14 06:42:51] Convert ID types for the GO_BP database
#> ℹ [2026-05-14 06:42:51] Connecting to the ensembl database...
#> ℹ [2026-05-14 06:42:51] Searching the dataset athaliana ...
#> ℹ [2026-05-14 06:42:52] Connecting to the dataset athaliana_eg_gene ...
#> ℹ [2026-05-14 06:42:52] Converting the geneIDs...
#> ℹ [2026-05-14 06:42:58] 21184 genes mapped with "tair_locus"
#> ℹ [2026-05-14 06:42:58] ==============================
#> ℹ                       21184 genes mapped
#> ℹ                       8 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-14 06:43:23] Convert ID types for the GO_CC database
#> ℹ [2026-05-14 06:43:23] Connecting to the ensembl database...
#> ℹ [2026-05-14 06:43:23] Searching the dataset athaliana ...
#> ℹ [2026-05-14 06:43:23] Connecting to the dataset athaliana_eg_gene ...
#> ℹ [2026-05-14 06:43:23] Converting the geneIDs...
#> ℹ [2026-05-14 06:43:30] 26896 genes mapped with "tair_locus"
#> ℹ [2026-05-14 06:43:30] ==============================
#> ℹ                       26896 genes mapped
#> ℹ                       6 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-14 06:44:01] Convert ID types for the GO_MF database
#> ℹ [2026-05-14 06:44:01] Connecting to the ensembl database...
#> ℹ [2026-05-14 06:44:01] Searching the dataset athaliana ...
#> ℹ [2026-05-14 06:44:02] Connecting to the dataset athaliana_eg_gene ...
#> ℹ [2026-05-14 06:44:02] Converting the geneIDs...
#> ℹ [2026-05-14 06:44:08] 26716 genes mapped with "tair_locus"
#> ℹ [2026-05-14 06:44:08] ==============================
#> ℹ                       26716 genes mapped
#> ℹ                       6 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-14 06:44:39] Convert ID types for the KEGG database
#> ℹ [2026-05-14 06:44:39] Connecting to the ensembl database...
#> ℹ [2026-05-14 06:44:39] Searching the dataset athaliana ...
#> ℹ [2026-05-14 06:44:39] Connecting to the dataset athaliana_eg_gene ...
#> ℹ [2026-05-14 06:44:40] Converting the geneIDs...
#> ℹ [2026-05-14 06:44:41] 5651 genes mapped with "entrez_id"
#> ℹ [2026-05-14 06:44:41] ==============================
#> ℹ                       5651 genes mapped
#> ℹ                       699 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-14 06:44:45] Convert ID types for the WikiPathway database
#> ℹ [2026-05-14 06:44:45] Connecting to the ensembl database...
#> ℹ [2026-05-14 06:44:45] Searching the dataset athaliana ...
#> ℹ [2026-05-14 06:44:45] Connecting to the dataset athaliana_eg_gene ...
#> ℹ [2026-05-14 06:44:45] Converting the geneIDs...
#> ℹ [2026-05-14 06:44:46] 626 genes mapped with "entrez_id"
#> ℹ [2026-05-14 06:44:46] ==============================
#> ℹ                       626 genes mapped
#> ℹ                       0 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-14 06:44:46] Convert ID types for the Chromosome database
#> ℹ [2026-05-14 06:44:46] Connecting to the ensembl database...
#> ℹ [2026-05-14 06:44:46] Searching the dataset athaliana ...
#> ℹ [2026-05-14 06:44:47] Connecting to the dataset athaliana_eg_gene ...
#> ℹ [2026-05-14 06:44:47] Converting the geneIDs...
#> ℹ [2026-05-14 06:44:53] 26925 genes mapped with "tair_locus"
#> ℹ [2026-05-14 06:44:53] ==============================
#> ℹ                       26925 genes mapped
#> ℹ                       491 genes unmapped
#> ℹ                       ==============================
head(
  db_list[["Arabidopsis_thaliana"]][["KEGG"]][["TERM2GENE"]]
)
#>       Term entrez_id       symbol ensembl_id
#> 1 ath00190    839579         PPa1  AT1G01050
#> 2 ath04712    839341          LHY  AT1G01060
#> 3 ath00010    839429 PDH-E1 ALPHA  AT1G01090
#> 4 ath00620    839429 PDH-E1 ALPHA  AT1G01090
#> 5 ath01210    839429 PDH-E1 ALPHA  AT1G01090
#> 6 ath01200    839429 PDH-E1 ALPHA  AT1G01090

# You can also build a custom database based on the gene sets you have
ccgenes <- CycGenePrefetch("Homo_sapiens")
#> ℹ [2026-05-14 06:45:19] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-05-14 06:45:19] Cell cycle gene prefetching completed "Homo_sapiens"
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
#> Error in data.frame(term = "S_genes", gene = ccgenes[["cc_S_genes"]]): arguments imply differing number of rows: 1, 0
str(custom_TERM2GENE)
#> Error: object 'custom_TERM2GENE' not found

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
#> ℹ [2026-05-14 06:45:19] Species: "Homo_sapiens"
#> Error: object 'custom_TERM2GENE' not found
ListDB(db = "CellCycle")
#>  [1] identifier version    comment    timestamp  date       db_version
#>  [7] db_name    file       Species    DB        
#> <0 rows> (or 0-length row.names)

db_list <- PrepareDB(species = "Mus_musculus", db = "CellCycle")
#> ℹ [2026-05-14 06:45:19] Species: "Mus_musculus"
head(
  db_list[["Mus_musculus"]][["CellCycle"]][["TERM2GENE"]]
)
#> NULL
```
