# Prepare databases and reference resources

Prepare species-specific databases and reference resources from
annotation packages, cached downloads, or supplied files.

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

  Character vector of database or resource selectors: `"GO"`, `"GO_BP"`,
  `"GO_CC"`, `"GO_MF"`, `"KEGG"`, `"WikiPathway"`, `"Reactome"`,
  `"CORUM"`, `"MP"`, `"DO"`, `"HPO"`, `"PFAM"`, `"Chromosome"`,
  `"GeneType"`, `"Enzyme"`, `"TF"`, `"CSPA"`, `"Surfaceome"`,
  `"SPRomeDB"`, `"VerSeDa"`, `"TFLink"`, `"hTFtarget"`, `"TRRUST"`,
  `"JASPAR"`, `"ENCODE"`, `"MSigDB"`, `"CellTalk"`, `"CellChat"`,
  `"CytoTRACE2"` or `"MSigDB_<collection>"`. A vector may combine
  sources. See the corresponding preparation function for
  source-specific settings.

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

A named list of prepared resources. Gene annotation databases are nested
by species and database, with `TERM2GENE` (gene-to-term mappings),
`TERM2NAME` (term names) and `version`. Additional resource fields
follow the return contract of the corresponding preparation function.

## See also

[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md),
[PrepareGO](https://mengxu98.github.io/scop/reference/PrepareGO.md),
[PrepareKEGG](https://mengxu98.github.io/scop/reference/PrepareKEGG.md),
[PrepareMSigDB](https://mengxu98.github.io/scop/reference/PrepareMSigDB.md),
[PrepareIREA](https://mengxu98.github.io/scop/reference/PrepareIREA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
databases <- PrepareDB(species = "Homo_sapiens", db = c("GO_BP", "KEGG"))
names(databases[["Homo_sapiens"]])
ListDB(species = "Homo_sapiens", db = c("GO_BP", "KEGG"))
} # }
```
