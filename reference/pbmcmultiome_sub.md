# A small human PBMC multiome example dataset

A near-balanced 500-cell subset of the PBMC multiome dataset from
SeuratData, containing paired `RNA` and `peaks` assays for package
examples and tests. The dataset keeps approximately equal numbers of
cells for each major PBMC cell type and retains the top 12000 accessible
peaks by total counts within the selected cells. When available, the
`peaks` assay stores a compact hg38 gene annotation derived from
`EnsDb.Hsapiens.v86` and collapsed to the longest transcript per gene.

## Format

A `Seurat` object.

## Source

Derived from the PBMC multiome reference data distributed through
[SeuratData](https://github.com/satijalab/seurat-data) /
`pbmcMultiome.SeuratData`, using the helper object
`test/data/pbmc_multiome_1k.rds` in this repository.

## Examples

``` r
if (interactive()) {
  source("test/data/create_pbmcmultiome_sub.R")
  pbmcmultiome_sub <- create_pbmcmultiome_sub()
  use_data <- get_namespace_fun("usethis", "use_data")
  use_data(
    pbmcmultiome_sub,
    compress = "xz",
    overwrite = TRUE
  )
}
```
