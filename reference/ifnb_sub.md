# A subsetted version of 'ifnb' datasets

Human PBMC control/IFNB-stimulated dataset

## Format

A `Seurat` object.

## Source

[paper_ifnb](https://www.nature.com/articles/nbt.4042)

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  if (!require("SeuratData", quietly = TRUE)) {
    pak::pak("satijalab/seurat-data")
  }
  library(SeuratData)
  library(Seurat)
  suppressWarnings(InstallData("ifnb"))
  data(ifnb)
  set.seed(11)
  cells_sub <- unlist(
    lapply(
      split(colnames(ifnb), ifnb$stim),
      function(x) sample(x, size = 1000)
    )
  )
  ifnb_sub <- subset(ifnb, cells = cells_sub)
  ifnb_sub <- ifnb_sub[Matrix::rowSums(
    GetAssayData5(
      ifnb_sub,
      assay = "RNA",
      layer = "counts"
    )
  ) > 0, ]
  ifnb_sub <- UpdateSeuratObject(ifnb_sub)
  usethis::use_data(ifnb_sub, compress = "xz")
}
} # }
```
