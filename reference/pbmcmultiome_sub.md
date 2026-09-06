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
[SeuratData](https://github.com/satijalab/seurat-data) /
`pbmcMultiome.SeuratData`.

## Examples

``` r
if (interactive()) {
  thisutils::check_r("satijalab/seurat-data")
  InstallData <- thisutils::get_namespace_fun("SeuratData", "InstallData")
  InstallData("pbmcMultiome")
  data(pbmcMultiome)
  pbmcMultiome <- UpdateSeuratObject(pbmcMultiome)

  set.seed(98)
  n_per_type <- ceiling(500 / length(unique(pbmcMultiome$cell_type)))
  cells_sub <- unlist(
    lapply(
      split(colnames(pbmcMultiome), pbmcMultiome$cell_type),
      function(x) sample(x, size = min(n_per_type, length(x)))
    ),
    use.names = FALSE
  )
  pbmcmultiome_sub <- subset(pbmcMultiome, cells = cells_sub)

  peak_counts <- GetAssayData5(pbmcmultiome_sub, assay = "peaks", layer = "counts")
  peaks_keep <- names(sort(Matrix::rowSums(peak_counts), decreasing = TRUE))[seq_len(12000)]
  rna_counts <- GetAssayData5(pbmcmultiome_sub, assay = "RNA", layer = "counts")
  genes_keep <- names(Matrix::rowSums(rna_counts))[Matrix::rowSums(rna_counts) >= 10]
  pbmcmultiome_sub <- subset(pbmcmultiome_sub, features = c(genes_keep, peaks_keep))

  use_data <- thisutils::get_namespace_fun("usethis", "use_data")
  use_data(
    pbmcmultiome_sub,
    compress = "xz",
    overwrite = TRUE
  )
}
```
