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
    "Chromosome", "GeneType", "Enzyme", "TF"),
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
  annotation databases should be prepared. Default is
  `c("Homo_sapiens", "Mus_musculus")`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Default is
  `c("GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF")`.

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

  Ensembl database version. If NULL, use the current release version.

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

## Details

The `PrepareDB` function prepares gene annotation databases for a given
species and set of annotation sources. It retrieves the necessary
information from various annotation packages or external resources and
organizes it into a list. The function also supports creating custom
databases based on user-provided gene sets.

## See also

[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
if (FALSE) { # \dontrun{
db_list <- PrepareDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
ListDB(
  species = "Homo_sapiens",
  db = "GO_BP"
)
head(
  db_list[["Homo_sapiens"]][["GO_BP"]][["TERM2GENE"]]
)

# Based on homologous gene conversion,
# prepare a gene annotation database that originally does not exist in the species.
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

# Prepare databases for other species
db_list <- PrepareDB(
  species = "Macaca_fascicularis",
  db = "GO_BP"
)
ListDB(
  species = "Macaca_fascicularis",
  db = "GO_BP"
)
head(
  db_list[["Macaca_fascicularis"]][["GO_BP"]][["TERM2GENE"]]
)

db_list <- PrepareDB(
  species = "Saccharomyces_cerevisiae",
  db = "GO_BP"
)
ListDB(
  species = "Saccharomyces_cerevisiae",
  db = "GO_BP"
)
head(
  db_list[["Saccharomyces_cerevisiae"]][["GO_BP"]][["TERM2GENE"]]
)

# Prepare databases for Arabidopsis (plant)
db_list <- PrepareDB(
  species = "Arabidopsis_thaliana",
  db = c(
    "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway",
    "ENZYME", "Chromosome"
  ),
  biomart = "plants_mart"
)
head(
  db_list[["Arabidopsis_thaliana"]][["KEGG"]][["TERM2GENE"]]
)

# You can also build a custom database based on the gene sets you have
ccgenes <- CycGenePrefetch("Homo_sapiens")
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
ListDB(db = "CellCycle")

db_list <- PrepareDB(species = "Mus_musculus", db = "CellCycle")
head(
  db_list[["Mus_musculus"]][["CellCycle"]][["TERM2GENE"]]
)
} # }
```
