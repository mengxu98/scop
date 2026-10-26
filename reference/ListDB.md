# List cached gene annotation databases

Lists the gene annotation databases cached by
[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md) and
[PrepareCCCDB](https://mengxu98.github.io/scop/reference/PrepareCCCDB.md)
in the R.cache root, optionally filtered by species and database name.
Each row describes one cached database with its version and creation
date.

## Usage

``` r
ListDB(species = c("Homo_sapiens", "Mus_musculus"), db = NULL)
```

## Arguments

- species:

  Species. Can be `"Homo_sapiens"` or `"Mus_musculus"`.

- db:

  Database names (for example `"GO_BP"`, `"KEGG"`, `"CellChat"`), or a
  regular expression. If `NULL`, all databases are listed.

## Value

A data frame with columns `Database`, `Species`, `Version`, and `Date`.

## Examples

``` r
ListDB(species = "Homo_sapiens")
#> [1] Database Species  Version  Date    
#> <0 rows> (or 0-length row.names)
ListDB(species = c("Homo_sapiens", "Mus_musculus"))
#>   Database      Species                Version                       Date
#> 1     CSPA Mus_musculus           CSPA nterm:1 2026-10-26 14:36:31.814562
#> 2       DO Mus_musculus       9.0.0 nterm:6595 2026-10-26 14:38:45.490315
#> 3    GO_BP Mus_musculus     3.23.0 nterm:14957 2026-10-26 14:38:07.894376
#> 4    GO_CC Mus_musculus      3.23.0 nterm:2065 2026-10-26 14:38:09.014323
#> 5       MP Mus_musculus 2026-10-26 nterm:10963 2026-10-26 14:38:40.990333
#> 6       TF Mus_musculus    AnimalTFDB4 nterm:2 2026-10-26 14:25:50.581074
ListDB(species = "Mus_musculus", db = "GO_BP")
#>   Database      Species            Version                       Date
#> 1    GO_BP Mus_musculus 3.23.0 nterm:14957 2026-10-26 14:38:07.894376
```
