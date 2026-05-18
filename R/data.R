#' @title A subsetted version of mouse 'pancreas' datasets
#'
#' @description
#' Mouse pancreatic endocrinogenesis dataset from \href{https://doi.org/10.1242/dev.173849}{Bastidas-Ponce et al. (2019)}.
#' A total of 1000 cells were downsampled to form the `pancreas_sub` dataset.
#'
#' @md
#' @format A `Seurat` object.
#' @concept data
#' @source
#' \href{https://scvelo.readthedocs.io/scvelo.datasets.pancreas/}{scvelo.datasets.pancreas},
#' \href{https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad}{endocrinogenesis_day15.h5ad}
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(Seurat)
#'   library(reticulate)
#'   PrepareEnv()
#'   check_python("scvelo")
#'   scv <- import("scvelo")
#'   adata <- scv$datasets$pancreas()
#'   pancreas <- adata_to_srt(adata)
#'   set.seed(98)
#'   cells <- sample(colnames(pancreas), size = 1000)
#'   pancreas_sub <- pancreas[, cells]
#'   pancreas_sub <- pancreas_sub[Matrix::rowSums(
#'     GetAssayData5(
#'       pancreas_sub,
#'       layer = "counts"
#'     )
#'   ) > 0, ]
#'   pancreas_sub[["CellType"]] <- pancreas_sub[["clusters_coarse"]]
#'   pancreas_sub[["SubCellType"]] <- pancreas_sub[["clusters"]]
#'   pancreas_sub[["clusters_coarse"]] <- pancreas_sub[["clusters"]] <- NULL
#'   pancreas_sub[["Phase"]] <- ifelse(
#'     pancreas_sub$S_score > pancreas_sub$G2M_score,
#'     "S",
#'     "G2M"
#'   )
#'   pancreas_sub[["Phase"]][apply(
#'     pancreas_sub[[]][, c("S_score", "G2M_score")],
#'     1,
#'     max
#'   ) < 0, ] <- "G1"
#'   pancreas_sub[["Phase", drop = TRUE]] <- factor(
#'     pancreas_sub[["Phase", drop = TRUE]],
#'     levels = c("G1", "S", "G2M")
#'   )
#'   pancreas_sub$CellType <- gsub("_", "-", pancreas_sub$CellType)
#'   pancreas_sub$CellType <- gsub(" ", "-", pancreas_sub$CellType)
#'   pancreas_sub$SubCellType <- gsub("_", "-", pancreas_sub$SubCellType)
#'   pancreas_sub$SubCellType <- gsub(" ", "-", pancreas_sub$SubCellType)
#'   pancreas_sub@reductions$X_pca <- NULL
#'   pancreas_sub@reductions$X_umap <- NULL
#'   use_data <- get_namespace_fun("usethis", "use_data")
#'   use_data(
#'     pancreas_sub,
#'     compress = "xz",
#'     overwrite = TRUE
#'   )
#' }
#' }
#' @name pancreas_sub
NULL

#' @title A subsetted version of human 'panc8' datasets
#'
#' @description
#' Human pancreatic islet cell datasets produced across four technologies,
#' SMART-Seq2 (E-MTAB-5061), CelSeq (GSE81076), CelSeq2 (GSE85241), and Fluidigm C1 (GSE86469),
#' from \href{https://github.com/satijalab/seurat-data}{SeuratData} package.
#' For each data set in `panc8`, 200 cells were downsampled to form the `panc8_sub` dataset.
#'
#' @md
#' @format A `Seurat` object.
#' @concept data
#' @source
#' \href{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/}{E-MTAB-5061},
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076}{GSE81076},
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241}{GSE85241},
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469}{GSE86469}
#'
#' @examples
#' if (interactive()) {
#'   data(pancreas_sub)
#'   check_r("satijalab/seurat-data")
#'   library(Seurat)
#'   InstallData <- get_namespace_fun("SeuratData", "InstallData")
#'   InstallData("panc8")
#'   data(panc8)
#'   panc8 <- UpdateSeuratObject(panc8)
#'   set.seed(98)
#'   cells_sub <- unlist(
#'     lapply(
#'       split(colnames(panc8), panc8$dataset),
#'       function(x) sample(x, size = 200)
#'     )
#'   )
#'   panc8_sub <- subset(panc8, cells = cells_sub)
#'   counts <- GetAssayData5(
#'     panc8_sub,
#'     layer = "counts"
#'   )
#'   panc8_sub <- CreateSeuratObject(
#'     counts = counts,
#'     meta.data = panc8_sub@meta.data
#'   )
#'   panc8_sub <- panc8_sub[Matrix::rowSums(counts) > 0, ]
#'   panc8_sub <- panc8_sub[toupper(
#'     rownames(panc8_sub)
#'   ) %in% toupper(
#'     rownames(pancreas_sub)
#'   ), ]
#'   panc8_sub$celltype <- gsub("_", "-", panc8_sub$celltype)
#'   panc8_sub$celltype <- gsub(" ", "-", panc8_sub$celltype)
#'   use_data <- get_namespace_fun("usethis", "use_data")
#'   use_data(
#'     panc8_sub,
#'     compress = "xz",
#'     overwrite = TRUE
#'   )
#' }
#' @name panc8_sub
NULL

#' @title A small human PBMC multiome example dataset
#'
#' @description
#' A near-balanced 500-cell subset of the PBMC multiome dataset from SeuratData,
#' containing paired `RNA` and `peaks` assays for package examples and tests.
#' The dataset keeps approximately equal numbers of cells for each major PBMC cell type and retains the
#' top 12000 accessible peaks by total counts within the selected cells.
#' When available, the `peaks` assay stores a compact hg38 gene annotation
#' derived from `EnsDb.Hsapiens.v86` and collapsed to the longest transcript
#' per gene.
#'
#' @md
#' @format A `Seurat` object.
#' @concept data
#' @source
#' Derived from the PBMC multiome reference data distributed through
#' \href{https://github.com/satijalab/seurat-data}{SeuratData} /
#' `pbmcMultiome.SeuratData`, using the helper object
#' `test/data/pbmc_multiome_1k.rds` in this repository.
#'
#' @examples
#' if (interactive()) {
#'   source("test/data/create_pbmcmultiome_sub.R")
#'   pbmcmultiome_sub <- create_pbmcmultiome_sub()
#'   use_data <- get_namespace_fun("usethis", "use_data")
#'   use_data(
#'     pbmcmultiome_sub,
#'     compress = "xz",
#'     overwrite = TRUE
#'   )
#' }
#' @name pbmcmultiome_sub
NULL

#' @title A human pancreas Visium spatial example dataset
#'
#' @description
#' A compact gene-filtered version of a human pancreatic intraepithelial
#' neoplasia (PanIN) 10x Visium dataset from GSE254829. The object keeps the
#' 1986 non-background tissue spots from sample GSM8058244 (PanIN-LG2), with a
#' `Spatial` assay, a `slice1` Visium image, and tissue coordinates in metadata
#' columns `x` and `y`. Metadata column `coda_label` stores the dominant CODA
#' microanatomical component for each spot, and `coda_score` stores its
#' percentage. Component percentage columns are stored with the `coda_` prefix,
#' and the matched CODA table is stored in `@tools$GSE254829_coda_table`. To
#' keep the package data small and directly usable with the bundled `panc8_sub`
#' reference, the object retains the top 5000 genes shared with `panc8_sub`,
#' ranked by total spatial counts.
#'
#' @md
#' @format A `Seurat` object with 5000 genes, 1986 spots, and one Visium image
#' named `slice1`.
#' @concept data
#' @source
#' Derived from the GSE254829 human PanIN 10x
#' Visium dataset:
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254829}{GSE254829}.
#' The package object uses the GEO supplementary files
#' `GSM8058244_PanIN-LG2.tar.gz` and
#' `GSE254829_codatable_may202024.csv.gz`.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' SeuratObject::Images(visium_human_pancreas_sub)
#' head(visium_human_pancreas_sub@meta.data[, c("x", "y")])
#' SpatialSpotPlot(visium_human_pancreas_sub, group.by = "coda_label")
#'
#' @name visium_human_pancreas_sub
NULL

#' @title Human pancreatic islet bulk RNA-seq example dataset
#'
#' @description
#' A full human pancreatic islet bulk RNA-seq `SummarizedExperiment` derived
#' from a brefeldin A perturbation study. The object keeps all samples from the
#' published islet arm and stores a symbol-level count matrix that can be used
#' directly in bulk DE and deconvolution examples together with the
#' bundled `panc8_sub` reference.
#'
#' @md
#' @format A `SummarizedExperiment` object with 19876 genes and 8 bulk RNA-seq
#' samples.
#' @concept data
#' @source
#' Derived from
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152615}{GSE152615}.
#' The bundled object is built from the supplementary matrix
#' `GSE152615_Rawcounts_filtered.txt.gz`. The published non-integer count values
#' are rounded to the nearest integer for count-based example workflows.
#'
#' @examples
#' data(islet_bulk)
#' SummarizedExperiment::assayNames(islet_bulk)
#' head(rownames(islet_bulk))
#' table(SummarizedExperiment::colData(islet_bulk)$condition)
#'
#' @name islet_bulk
NULL

#' @title Excluded words in keyword enrichment analysis and extraction
#'
#' @md
#' @description
#' The variable `words_excluded` represents the words that are excluded during keyword enrichment analysis or keyword extraction process.
#' These mainly include words that are excessively redundant or of little value.
#'
#' @concept data
#' @examples
#' if (interactive()) {
#'   words_excluded <- c(
#'     "the", "is", "and", "or", "a",
#'     "in", "on", "under", "between", "of",
#'     "through", "via", "along", "that",
#'     "for", "with", "within", "without",
#'     "cell", "cellular", "dna", "rna",
#'     "protein", "peptide", "amino", "acid",
#'     "development", "involved", "organization", "system",
#'     "regulation", "regulated", "positive", "negative",
#'     "response", "process", "processing", "small", "large", "change"
#'   )
#'   use_data <- get_namespace_fun("usethis", "use_data")
#'   use_data(words_excluded, compress = "xz")
#' }
#' @name words_excluded
NULL

#' @title Reference datasets for cell type annotation in single-cell RNA data
#'
#' @concept data
#' @source
#' \href{https://github.com/ggjlab/scMCA}{scMCA}
#' @examples
#' if (interactive()) {
#'   library(Seurat)
#'   check_r(c("ggjlab/scMCA"))
#'   ref_scMCA <- NormalizeData(get("ref.expr", envir = asNamespace("scMCA")))
#'   Encoding(colnames(ref_scMCA)) <- "latin1"
#'   colnames(ref_scMCA) <- iconv(colnames(ref_scMCA), "latin1", "UTF-8")
#'   # get_namespace_fun("usethis", "use_data")(ref_scMCA, compress = "xz")
#' }
#' @name ref_scMCA
NULL
