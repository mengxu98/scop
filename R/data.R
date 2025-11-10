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
#'   usethis::use_data(
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
#' \dontrun{
#' if (interactive()) {
#'   data(pancreas_sub)
#'   if (!require("SeuratData", quietly = TRUE)) {
#'     pak::pak("satijalab/seurat-data")
#'   }
#'   library(SeuratData)
#'   library(Seurat)
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
#'   usethis::use_data(
#'     panc8_sub,
#'     compress = "xz",
#'     overwrite = TRUE
#'   )
#' }
#' }
#' @name panc8_sub
NULL

#' @title A subsetted version of 'ifnb' datasets
#'
#' @description
#' Human PBMC control/IFNB-stimulated dataset
#'
#' @md
#' @format A `Seurat` object.
#' @concept data
#' @source \href{https://www.nature.com/articles/nbt.4042}{paper_ifnb}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   if (!require("SeuratData", quietly = TRUE)) {
#'     pak::pak("satijalab/seurat-data")
#'   }
#'   library(SeuratData)
#'   library(Seurat)
#'   suppressWarnings(InstallData("ifnb"))
#'   data(ifnb)
#'   set.seed(11)
#'   cells_sub <- unlist(
#'     lapply(
#'       split(colnames(ifnb), ifnb$stim),
#'       function(x) sample(x, size = 1000)
#'     )
#'   )
#'   ifnb_sub <- subset(ifnb, cells = cells_sub)
#'   ifnb_sub <- ifnb_sub[Matrix::rowSums(
#'     GetAssayData5(
#'       ifnb_sub,
#'       assay = "RNA",
#'       layer = "counts"
#'     )
#'   ) > 0, ]
#'   ifnb_sub <- UpdateSeuratObject(ifnb_sub)
#'   usethis::use_data(ifnb_sub, compress = "xz")
#' }
#' }
#' @name ifnb_sub
NULL

#' Excluded words in keyword enrichment analysis and extraction
#'
#' The variable "words_excluded" represents the words that are excluded during keyword enrichment analysis or keyword extraction process.
#' These mainly include words that are excessively redundant or of little value.
#'
#' @concept data
#' @examples
#' \dontrun{
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
#'   usethis::use_data(words_excluded, compress = "xz")
#' }
#' }
#' @name words_excluded
NULL

#' Reference datasets for cell type annotation in single-cell RNA data
#'
#' @concept data
#' @source
#' \href{https://github.com/ggjlab/scHCL}{scHCL},
#' \href{https://github.com/ggjlab/scMCA}{scMCA},
#' \href{https://github.com/ggjlab/scZCL}{scZCL}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(Seurat)
#'   check_r(c("ggjlab/scZCL", "ggjlab/scHCL", "ggjlab/scMCA"))
#'   ref_scHCL <- NormalizeData(scHCL::ref.expr)
#'   ref_scMCA <- NormalizeData(scMCA::ref.expr)
#'   ref_scZCL <- NormalizeData(scZCL::ref.expr)
#'   Encoding(colnames(ref_scHCL)) <- "latin1"
#'   colnames(ref_scHCL) <- iconv(colnames(ref_scHCL), "latin1", "UTF-8")
#'   Encoding(colnames(ref_scMCA)) <- "latin1"
#'   colnames(ref_scMCA) <- iconv(colnames(ref_scMCA), "latin1", "UTF-8")
#'   Encoding(colnames(ref_scZCL)) <- "latin1"
#'   colnames(ref_scZCL) <- iconv(colnames(ref_scZCL), "latin1", "UTF-8")
#'   # usethis::use_data(ref_scHCL, compress = "xz")
#'   # usethis::use_data(ref_scMCA, compress = "xz")
#'   # usethis::use_data(ref_scZCL, compress = "xz")
#' }
#' }
#' @usage ref_scHCL
#' @usage ref_scMCA
#' @usage ref_scZCL
#' @name ref_scHCL
#' @aliases ref_scHCL ref_scMCA ref_scZCL
NULL
