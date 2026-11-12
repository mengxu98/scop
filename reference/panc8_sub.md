# A subsetted version of human 'panc8' datasets

Human pancreatic islet cell datasets produced across four technologies,
SMART-Seq2 (E-MTAB-5061), CelSeq (GSE81076), CelSeq2 (GSE85241), and
Fluidigm C1 (GSE86469), from
[SeuratData](https://github.com/satijalab/seurat-data) package. For each
data set in `panc8`, 200 cells were downsampled to form the `panc8_sub`
dataset.

## Format

A `Seurat` object.

## Source

[E-MTAB-5061](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/),
[GSE81076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076),
[GSE85241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241),
[GSE86469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469)

## Examples

``` r
if (interactive()) {
  data(pancreas_sub)

  archive <- tempfile(fileext = ".tar.gz")
  utils::download.file(
    "https://seurat.nygenome.org/src/contrib/panc8.SeuratData_3.0.2.tar.gz",
    archive,
    mode = "wb"
  )
  exdir <- tempfile()
  utils::untar(archive, files = "panc8.SeuratData/data/panc8.rda", exdir = exdir)
  load(file.path(exdir, "panc8.SeuratData", "data", "panc8.rda"))
  panc8 <- SeuratObject::UpdateSeuratObject(panc8)
  set.seed(98)
  cells_sub <- unlist(
    lapply(
      split(colnames(panc8), panc8$dataset),
      function(x) sample(x, size = 200)
    )
  )
  panc8_sub <- subset(panc8, cells = cells_sub)
  counts <- GetAssayData5(
    panc8_sub,
    layer = "counts"
  )
  panc8_sub <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = panc8_sub@meta.data
  )
  panc8_sub <- panc8_sub[Matrix::rowSums(counts) > 0, ]
  panc8_sub <- panc8_sub[toupper(
    rownames(panc8_sub)
  ) %in% toupper(
    rownames(pancreas_sub)
  ), ]
  panc8_sub$celltype <- gsub("_", "-", panc8_sub$celltype)
  panc8_sub$celltype <- gsub(" ", "-", panc8_sub$celltype)
  usethis::use_data(
    panc8_sub,
    compress = "xz",
    overwrite = TRUE
  )
}
```
