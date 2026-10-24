# Prepare Immune Dictionary references

Prepare cell-type reference objects and signatures used by
[`RunIREA()`](https://mengxu98.github.io/scop/reference/RunIREA.md).

## Usage

``` r
PrepareIREA(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "IREA_NK_cell",
  db_update = FALSE,
  data_dir = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- species:

  `"Homo_sapiens"` or `"Mus_musculus"`.

- db:

  Reference selectors such as `"IREA_NK_cell"` or `"IREA_Macrophage"`.
  The suffix is the exact portal cell-type identifier.

- db_update:

  Force a refresh. `FALSE` loads the cache when available.

- data_dir:

  Directory or named list of local source files. Searches
  `data_dir/<db>/` then `data_dir`. Named lists override a path, e.g.
  `list(MSigDB = "~/db/msigdb")`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Passed to helper functions.

## Value

A named species list with an `irea_reference` under each selector. Each
reference contains `object`, `cytokine`, `polarization`, `cell_type`,
`species`, `paths`, `checksum_md5` and `provenance`. Expression rows are
mouse genes; columns are reference cells in their original order.

## Details

Missing files are downloaded from the official portal. Source assets are
unversioned; the cache manifest records their checksums and rejects
changed files unless `db_update = TRUE`. Named `data_dir` lists may use
the selector or `"IREA"`. Each file is searched under the selector
subdirectory before the parent directory. Human input uses partial
orthologue mappings; the expression reference remains mouse lymph-node
cells. Files are not bundled.

## References

Cui, Ang; Huang, Teddy; Li, Shuqiang; Ma, Aileen; Perez, Jorge L.;
Sander, Chris; Keskin, Derin B.; Wu, Catherine J.; Fraenkel, Ernest;
Hacohen, Nir. Dictionary of immune responses to cytokines at single-cell
resolution. Nature 625, 377-384 (2024). doi:10.1038/s41586-023-06816-9.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[RunIREA](https://mengxu98.github.io/scop/reference/RunIREA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
references <- PrepareIREA(db = "IREA_NK_cell", species = "Mus_musculus")
result <- RunIREA(c("Isg15", "Ifit3", "Bst2"),
  reference = references[["Mus_musculus"]][["IREA_NK_cell"]],
  analysis = "cell_polarization"
)
IREAPlot(result, plot_type = "radar")
} # }
```
