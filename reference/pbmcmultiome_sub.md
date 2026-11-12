# A small human PBMC multiome example dataset

A near-balanced 500-cell subset of the PBMC multiome dataset from
SeuratData, containing paired `RNA` and `peaks` assays for package
examples and tests. The dataset keeps approximately equal numbers of
cells for each major PBMC cell type, retains the top 12000 accessible
peaks by total counts within the selected cells, and keeps RNA genes
with at least 10 total counts in the subset (12405 genes). When
available, the `peaks` assay stores a compact hg38 gene annotation
derived from `EnsDb.Hsapiens.v86` and collapsed to the longest
transcript per gene.

## Format

A `Seurat` object with 12405 RNA genes, 12000 peaks, and 500 cells.

## Source

Derived from the PBMC multiome reference data distributed through
[SeuratData](https://github.com/satijalab/seurat-data).

## Examples

``` r
if (interactive()) {
  archive <- tempfile(fileext = ".tar.gz")
  utils::download.file(
    "https://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.4.tar.gz",
    archive,
    mode = "wb"
  )
  exdir <- tempfile()
  utils::untar(archive, files = "pbmcMultiome.SeuratData/data", exdir = exdir)
  load(file.path(exdir, "pbmcMultiome.SeuratData", "data", "pbmc.rna.rda"))
  load(file.path(exdir, "pbmcMultiome.SeuratData", "data", "pbmc.atac.rda"))
  pbmcMultiome <- SeuratObject::UpdateSeuratObject(pbmc.rna)
  pbmcMultiome[["peaks"]] <- SeuratObject::UpdateSeuratObject(pbmc.atac)[["ATAC"]]

  groups <- c(
    "Intermediate B" = "B", "Memory B" = "B", "Naive B" = "B",
    "CD4 Naive" = "CD4 T", "CD4 TCM" = "CD4 T", "CD4 TEM" = "CD4 T",
    Treg = "CD4 T", "CD8 Naive" = "CD8 T", "CD8 TEM_1" = "CD8 T",
    "CD8 TEM_2" = "CD8 T", gdT = "CD8 T", MAIT = "CD8 T",
    cDC = "DC", pDC = "DC", "CD14 Mono" = "Mono", "CD16 Mono" = "Mono",
    NK = "NK"
  )
  pbmcMultiome$SubCellType <- as.character(pbmcMultiome$seurat_annotations)
  pbmcMultiome$CellType <- unname(groups[pbmcMultiome$SubCellType])

  set.seed(98)
  cell_types <- split(colnames(pbmcMultiome), pbmcMultiome$CellType)
  n_per_type <- rep(500 %/% length(cell_types), length(cell_types))
  extra <- seq_len(500 %% length(cell_types))
  n_per_type[extra] <- n_per_type[extra] + 1L
  cells_sub <- unlist(Map(sample, cell_types, n_per_type), use.names = FALSE)
  pbmcmultiome_sub <- subset(pbmcMultiome, cells = cells_sub)

  peak_counts <- GetAssayData5(pbmcmultiome_sub, assay = "peaks", layer = "counts")
  peaks_keep <- names(sort(Matrix::rowSums(peak_counts), decreasing = TRUE))[seq_len(12000)]
  rna_counts <- GetAssayData5(pbmcmultiome_sub, assay = "RNA", layer = "counts")
  genes_keep <- names(Matrix::rowSums(rna_counts))[Matrix::rowSums(rna_counts) >= 10]
  pbmcmultiome_sub <- subset(pbmcmultiome_sub, features = c(genes_keep, peaks_keep))

  usethis::use_data(
    pbmcmultiome_sub,
    compress = "xz",
    overwrite = TRUE
  )
}
```
