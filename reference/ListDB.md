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
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.23.0 nterm:14957|Mus_musculus-GO_BP 1782271428 2026-06-24 03:23:47.768831
#>           db_version            db_name
#> 1 3.23.0 nterm:14957 Mus_musculus-GO_BP
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/a759db6d549ed9ffa83592264ac75787.Rcache
#>        Species    DB
#> 1 Mus_musculus GO_BP
```
