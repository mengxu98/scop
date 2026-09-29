# Download only the source resources required for a local IREA evaluation.
# Run: Rscript inst/irea/irea_download.R D:/scop-dev/datasets/IREA
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Supply a reference cache directory.")
dest <- args[[1]]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
base <- "https://www.immune-dictionary.org/static/"
cells <- c("B_cell","cDC1","cDC2","Langerhans","Macrophage","MigDC",
           "Monocyte","Neutrophil","NK_cell","pDC","T_cell_CD4",
           "T_cell_CD8","T_cell_gd","Treg")
paths <- c(paste0("downloadableData/ligands-seurat-",cells,".RDS"),
           paste0("dataFiles/SuppTable",c("3_Cytokine_Signatures",
             "3_Cytokine_Signatures_Human","7_Polarization_Signatures",
             "7_Polarization_Signatures_Human"),".xlsx"),
           "exampleFiles/gene_test.xlsx")
manifest <- lapply(paths,function(path) {
  target <- file.path(dest,gsub("/","_",path,fixed=TRUE))
  if (!file.exists(target)) utils::download.file(paste0(base,path),target,mode="wb",quiet=TRUE)
  data.frame(source_url=paste0(base,path),local_file=target,
             downloaded_or_checked=as.character(Sys.Date()),
             bytes=unname(file.info(target)$size),
             checksum_md5=unname(tools::md5sum(target)))
})
utils::write.csv(do.call(rbind,manifest),file.path(dest,"reference_manifest.csv"),row.names=FALSE)
