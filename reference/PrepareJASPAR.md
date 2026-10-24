# Prepare JASPAR databases

Prepare JASPAR annotations, reusing the common cache and
species/identifier conversion pipeline.

## Usage

``` r
PrepareJASPAR(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "JASPAR",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest",
  db_update = FALSE,
  data_dir = NULL,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  biomart = NULL,
  max_tries = 5,
  verbose = TRUE,
  ...
)
```

## Arguments

- species:

  `"Homo_sapiens"` or `"Mus_musculus"`.

- db:

  Character selector; must be `"JASPAR"`.

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

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to helper functions.

## Value

A named species list containing the selected databases. Annotation
entries contain `TERM2GENE`, `TERM2NAME` and `version`.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
if (FALSE) { # \dontrun{
databases <- PrepareJASPAR(species = "Homo_sapiens", db = "JASPAR")
names(databases[["Homo_sapiens"]])
} # }
```
