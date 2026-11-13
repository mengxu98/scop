run_case <- function(output) {
  data <- new.env(parent = emptyenv())
  utils::data("pancreas_sub", package = "scop", envir = data)
  object <- data$pancreas_sub[, seq_len(120L)]
  assay <- SeuratObject::DefaultAssay(object)
  object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
  expression <- SeuratObject::LayerData(object, assay = assay, layer = "data")
  variance <- Matrix::rowMeans(expression^2) - Matrix::rowMeans(expression)^2
  genes <- names(sort(variance, decreasing = TRUE))[seq_len(24L)]
  scop::RunGRN(object, assay = assay, layer = "data", methods = "ppcor",
    regulators = genes[1:8], targets = genes[9:24], return_data = TRUE,
    work_dir = file.path(output, "work"), verbose = FALSE)
}
