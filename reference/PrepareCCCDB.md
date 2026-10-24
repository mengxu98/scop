# Prepare cell-cell communication databases

Prepares the ligand-receptor interaction databases used by the cell-cell
communication wrappers: the CellTalkDB ligand-receptor pairs
(`"CellTalk"`) and the CellChat curated ligand-receptor database
(`"CellChat"`). Each database is returned as a `TERM2GENE`/`TERM2NAME`
mapping (with a `"ligand_*"`/`"receptor_*"` term convention) and cached
with [R.cache](https://rdrr.io/pkg/R.cache/man/R.cache-package.html) so
that [ListDB](https://mengxu98.github.io/scop/reference/ListDB.md) can
list it like any other annotation database.

## Usage

``` r
PrepareCCCDB(
  species = c("Homo_sapiens", "Mus_musculus", "Danio_rerio"),
  db = c("CellChat", "CellTalk"),
  convert_species = TRUE,
  data_dir = NULL,
  db_version = "latest",
  db_update = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- species:

  `"Homo_sapiens"` or `"Mus_musculus"`.

- db:

  Cell-cell communication databases to prepare. Can be `"CellTalk"`
  and/or `"CellChat"`.

- convert_species:

  Use a species-converted database when the annotation is missing for
  `species`.

- data_dir:

  Directory or named list of local source files. Searches
  `data_dir/<db>/` then `data_dir`. Named lists override a path, e.g.
  `list(MSigDB = "~/db/msigdb")`.

- db_version:

  Database version to retrieve.

- db_update:

  Whether the databases should be forcefully updated. If `FALSE`, cached
  databases are reused when available.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to helper functions.

## Value

A list with the same structure as
[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md): for
each species a named list of databases, each with `TERM2GENE`,
`TERM2NAME`, and `version` entries.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListCCCDB](https://mengxu98.github.io/scop/reference/ListCCCDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ccc_db <- PrepareCCCDB(
  species = "Homo_sapiens",
  db = "CellChat"
)
head(ccc_db[["Homo_sapiens"]][["CellChat"]][["TERM2GENE"]])
} # }
```
