# Prepare custom gene annotation databases

Build custom mappings using the common cache and species/identifier
conversion pipeline.

## Usage

``` r
PrepareCustomDB(
  species = c("Homo_sapiens", "Mus_musculus"),
  db,
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

A named species list containing the selected annotation database.

## See also

\[PrepareDB\], \[ListDB\]

## Examples

``` r
mappings <- data.frame(Term = c("Response", "Response"), symbol = c("Isg15", "Ifit3"))
databases <- PrepareCustomDB(
  species = "Mus_musculus", db = "Response", db_IDtypes = "symbol",
  custom_TERM2GENE = mappings, custom_species = "Mus_musculus",
  custom_IDtype = "symbol", custom_version = "v1", verbose = FALSE
)
databases[["Mus_musculus"]][["Response"]]$TERM2GENE
#>       Term symbol
#> 1 Response  Isg15
#> 2 Response  Ifit3
```
