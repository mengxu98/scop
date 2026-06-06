# List cached databases

Retrieves information about databases based on a given species and
database name.

## Usage

``` r
ListDB(species = c("Homo_sapiens", "Mus_musculus"), db = NULL)
```

## Arguments

- species:

  A character vector of species for which to retrieve database
  information. Default is `c("Homo_sapiens", "Mus_musculus")`.

- db:

  The pattern to match against the database names. Default is `NULL`,
  which matches all databases.

## Value

A data frame containing information about the databases, including a
`Species` column and a `DB` column.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md)

## Examples

``` r
if (requireNamespace("R.cache", quietly = TRUE)) {
  ListDB(species = "Homo_sapiens")
  ListDB(species = c("Homo_sapiens", "Mus_musculus"))
  ListDB(species = "Mus_musculus", db = "GO_BP")
}
```
