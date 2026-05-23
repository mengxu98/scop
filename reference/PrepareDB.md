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
#> ℹ [2026-05-23 14:02:04] Species: "Homo_sapiens"
#> ℹ [2026-05-23 14:02:04] Loading cached: GO_BP version: 3.23.1 nterm:14209 created: 2026-05-23 13:07:27
#> ℹ [2026-05-23 14:02:05] Convert ID types for the GO_BP database
#> ℹ [2026-05-23 14:02:06] Converted ID types using local annotation package org.Hs.eg.db
ListDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.23.1 nterm:14209|Homo_sapiens-GO_BP 1779544927 2026-05-23 14:02:07.022842
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
#> ℹ [2026-05-23 14:02:07] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-05-23 14:02:07] Cell cycle gene prefetching completed "Homo_sapiens"
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
#> ℹ [2026-05-23 14:02:07] Species: "Homo_sapiens"
#> ℹ [2026-05-23 14:02:07] Convert ID types for the CellCycle database
#> ℹ [2026-05-23 14:02:07] Converted ID types using local annotation package org.Hs.eg.db
#> ℹ [2026-05-23 14:02:07] Species: "Mus_musculus"
#> ! [2026-05-23 14:02:07] Use the "Homo_sapiens" annotation to create the "CellCycle" database for "Mus_musculus"
#> ℹ [2026-05-23 14:02:07] Convert species for the CellCycle database
#> ℹ [2026-05-23 14:02:07] Connect to the Ensembl archives...
#> ℹ [2026-05-23 14:02:08] Using the 115 version of ensembl database...
#> ℹ [2026-05-23 14:02:08] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-05-23 14:02:12] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-23 14:02:12] Get errors when connecting with ensembl database...
#> ! [2026-05-23 14:02:13] Retrying...
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying www mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-05-23 14:02:16] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-23 14:02:16] Get errors when connecting with ensembl database...
#> ! [2026-05-23 14:02:17] Retrying...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> ! [2026-05-23 14:02:20] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-23 14:02:20] Get errors when connecting with ensembl database...
#> ! [2026-05-23 14:02:21] Retrying...
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Warning: Invalid mirror. Select a mirror from [www, useast, asia].
#> Default when no mirror is specified is to use www.ensembl.org which may be automatically redirected.
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> Ensembl site unresponsive, trying www mirror
#> Ensembl site unresponsive, trying asia mirror
#> Ensembl site unresponsive, trying useast mirror
#> ! [2026-05-23 14:02:24] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-23 14:02:24] Get errors when connecting with ensembl database...
#> ! [2026-05-23 14:02:25] Retrying...
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
#> ! [2026-05-23 14:02:28] <simpleError: Your query has been redirected to https://status.ensembl.org indicating this Ensembl service is currently unavailable.
#> !                       Look at ?useEnsembl for details on how to try a mirror site.>
#> ! [2026-05-23 14:02:28] Get errors when connecting with ensembl database...
#> Error in try_get(expr = {    if (!is.null(mirror)) {        biomaRt::useEnsembl(biomart = "ensembl", mirror = mirror)    }    else {        mart_try <- tryCatch(biomaRt::useMart(biomart = "ensembl",             host = url), error = function(e) e)        if (inherits(mart_try, "error")) {            mirror_candidates <- c("useast", "uswest", "asia",                 "www")            for (mirror_i in mirror_candidates) {                mart_mirror <- tryCatch(biomaRt::useEnsembl(biomart = "ensembl",                   mirror = mirror_i), error = function(e) e)                if (!inherits(mart_mirror, "error")) {                  mart_try <- mart_mirror                  break                }            }        }        if (inherits(mart_try, "error")) {            stop(mart_try)        }        mart_try    }}, max_tries = max_tries, error_message = "Get errors when connecting with ensembl database..."): <simpleError: Your query has been redirected to
#> https://status.ensembl.org indicating this Ensembl service is currently
#> unavailable. Look at ?useEnsembl for details on how to try a mirror site.>
ListDB(db = "CellCycle")
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                    comment  timestamp
#> 1 Seurat_v5 nterm:2|Homo_sapiens-CellCycle 1779544928
#>                         date        db_version                db_name
#> 1 2026-05-23 14:02:07.894984 Seurat_v5 nterm:2 Homo_sapiens-CellCycle
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/a6aa81007b9564b5bf1f3fa5dc7997fa.Rcache
#>        Species        DB
#> 1 Homo_sapiens CellCycle

db_list <- PrepareDB(species = "Mus_musculus", db = "CellCycle")
#> ℹ [2026-05-23 14:02:29] Species: "Mus_musculus"
head(
  db_list[["Mus_musculus"]][["CellCycle"]][["TERM2GENE"]]
)
#> NULL
```
