# List cached databases

Retrieves information about databases based on a given species and
database name.

## Usage

``` r
ListDB(species = "Homo_sapiens", db = NULL)
```

## Arguments

- species:

  The species for which to retrieve database information. Default is
  `"Homo_sapiens"`.

- db:

  The pattern to match against the database names. Default is `NULL`,
  which matches all databases.

## Value

A data frame containing information about the databases.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md)

## Examples

``` r
ListDB(species = "Homo_sapiens")
#> Error in loadNamespace(x): there is no package called ‘R.cache’
ListDB(species = "Mus_musculus", db = "GO_BP")
#> Error in loadNamespace(x): there is no package called ‘R.cache’
```
