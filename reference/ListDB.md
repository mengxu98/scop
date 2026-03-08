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
ListDB(species = "Homo_sapiens")
#>  [1] identifier version    comment    timestamp  date       db_version
#>  [7] db_name    file       Species    DB        
#> <0 rows> (or 0-length row.names)
ListDB(species = c("Homo_sapiens", "Mus_musculus"))
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#> 2 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#> 3 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.22.0 nterm:15169|Mus_musculus-GO_BP 1772954150 2026-03-08 07:15:49.623996
#> 2        CSPA nterm:1|Mus_musculus-CSPA 1772952407 2026-03-08 06:46:46.514255
#> 3   AnimalTFDB4 nterm:2|Mus_musculus-TF 1772952396 2026-03-08 06:46:36.204431
#>            db_version            db_name
#> 1  3.22.0 nterm:15169 Mus_musculus-GO_BP
#> 2        CSPA nterm:1  Mus_musculus-CSPA
#> 3 AnimalTFDB4 nterm:2    Mus_musculus-TF
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/4363ecdce2b08b4d38a5c290b4f4ae60.Rcache
#> 2 /home/runner/.cache/R/R.cache/267624ab3d48b245ce01deac5165a7b1.Rcache
#> 3 /home/runner/.cache/R/R.cache/82f852e1d9ed5c0a5ad672b843c18742.Rcache
#>        Species    DB
#> 1 Mus_musculus GO_BP
#> 2 Mus_musculus  CSPA
#> 3 Mus_musculus    TF
ListDB(species = "Mus_musculus", db = "GO_BP")
#>                                                         identifier version
#> 1 Rcache v0.1.7 (R package R.cache by Henrik Bengtsson)              0.1.7
#>                                 comment  timestamp                       date
#> 1 3.22.0 nterm:15169|Mus_musculus-GO_BP 1772954150 2026-03-08 07:15:49.623996
#>           db_version            db_name
#> 1 3.22.0 nterm:15169 Mus_musculus-GO_BP
#>                                                                    file
#> 1 /home/runner/.cache/R/R.cache/4363ecdce2b08b4d38a5c290b4f4ae60.Rcache
#>        Species    DB
#> 1 Mus_musculus GO_BP
```
