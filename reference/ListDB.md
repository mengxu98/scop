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
#> [1] identifier version    comment    timestamp  date       db_version db_name   
#> [8] file      
#> <0 rows> (or 0-length row.names)
ListDB(species = "Mus_musculus", db = "GO_BP")
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.22.0 nterm:15169|Mus_musculus-GO_BP 1769052666 2026-01-22 03:31:05.647468
#>           db_version            db_name
#> 1 3.22.0 nterm:15169 Mus_musculus-GO_BP
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/4363ecdce2b08b4d38a5c290b4f4ae60.Rcache
```
