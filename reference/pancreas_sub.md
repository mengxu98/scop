# A subsetted version of mouse 'pancreas' datasets

Mouse pancreatic endocrinogenesis dataset from [Bastidas-Ponce et al.
(2019)](https://doi.org/10.1242/dev.173849). A total of 1000 cells were
downsampled to form the `pancreas_sub` dataset.

## Format

A `Seurat` object.

## Source

[scvelo.datasets.pancreas](https://scvelo.readthedocs.io/scvelo.datasets.pancreas/),
[endocrinogenesis_day15.h5ad](https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad)

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  library(Seurat)
  library(reticulate)
  PrepareEnv()
  check_python("scvelo")
  scv <- import("scvelo")
  adata <- scv$datasets$pancreas()
  pancreas <- adata_to_srt(adata)
  set.seed(98)
  cells <- sample(colnames(pancreas), size = 1000)
  pancreas_sub <- pancreas[, cells]
  pancreas_sub <- pancreas_sub[Matrix::rowSums(
    GetAssayData5(
      pancreas_sub,
      layer = "counts"
    )
  ) > 0, ]
  pancreas_sub[["CellType"]] <- pancreas_sub[["clusters_coarse"]]
  pancreas_sub[["SubCellType"]] <- pancreas_sub[["clusters"]]
  pancreas_sub[["clusters_coarse"]] <- pancreas_sub[["clusters"]] <- NULL
  pancreas_sub[["Phase"]] <- ifelse(
    pancreas_sub$S_score > pancreas_sub$G2M_score,
    "S",
    "G2M"
  )
  pancreas_sub[["Phase"]][apply(
    pancreas_sub[[]][, c("S_score", "G2M_score")],
    1,
    max
  ) < 0, ] <- "G1"
  pancreas_sub[["Phase", drop = TRUE]] <- factor(
    pancreas_sub[["Phase", drop = TRUE]],
    levels = c("G1", "S", "G2M")
  )
  pancreas_sub$CellType <- gsub("_", "-", pancreas_sub$CellType)
  pancreas_sub$CellType <- gsub(" ", "-", pancreas_sub$CellType)
  pancreas_sub$SubCellType <- gsub("_", "-", pancreas_sub$SubCellType)
  pancreas_sub$SubCellType <- gsub(" ", "-", pancreas_sub$SubCellType)
  pancreas_sub@reductions$X_pca <- NULL
  pancreas_sub@reductions$X_umap <- NULL
  usethis::use_data(
    pancreas_sub,
    compress = "xz",
    overwrite = TRUE
  )
}
} # }
```
