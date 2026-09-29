# Run with: Rscript inst/irea/irea_demo.R D:/scop-dev/datasets/IREA output-directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply reference and output directories.")
reference_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pkgload::load_all(".", quiet = TRUE, compile = FALSE)
ref <- PrepareIREAReference(reference_dir, "NK_cell", "Mouse")
genes <- c("Ifitm3", "Isg15", "Ifit3", "Bst2", "Slfn5",
           "Isg20", "Phf11b", "Zbp1", "Rtp4")
example <- readxl::read_excel(file.path(reference_dir,"exampleFiles_gene_test.xlsx"))
for (analysis in c("cytokine_response", "cell_polarization")) {
  for (method in c("score", "hypergeometric")) {
    out <- RunIREA(ref, genes = genes, analysis = analysis, method = method)
    stem <- paste("NK", "genes", analysis, method, sep = "_")
    utils::write.csv(out$table, file.path(output_dir, paste0(stem,".csv")), row.names = FALSE)
    grDevices::pdf(file.path(output_dir,paste0(stem,".pdf")),width=11,height=8)
    print(IREAPlot(out, if (analysis == "cytokine_response") "compass" else "radar"))
    grDevices::dev.off()
  }
  for (column in names(example)[-1]) {
    out <- RunIREA(ref, matrix = example, contrast = column, analysis = analysis)
    stem <- paste("NK", column, analysis, sep = "_")
    utils::write.csv(out$table, file.path(output_dir,paste0(stem,".csv")),row.names=FALSE)
    grDevices::pdf(file.path(output_dir,paste0(stem,".pdf")),width=11,height=8)
    print(IREAPlot(out, if (analysis == "cytokine_response") "compass" else "radar"))
    grDevices::dev.off()
  }
}
