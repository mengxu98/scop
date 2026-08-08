# End-to-end consistency between the scop CNV wrappers and the original
# backend methods, run with real bundled example datasets (copykat::exp.rawdata
# breast-cancer scRNA-seq, infercnv example, numbat ATC2 example). No simulated
# data and no mocked backends: the wrappers and the original pipelines receive
# the same real input and must return the same scientific results.

cnv_real_breast_seurat <- function() {
  data(exp.rawdata, package = "copykat")
  srt <- Seurat::CreateSeuratObject(counts = as.matrix(exp.rawdata))
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt
}

cnv_real_infercnv_inputs <- function() {
  data(infercnv_data_example, package = "infercnv")
  data(infercnv_annots_example, package = "infercnv")
  data(infercnv_genes_example, package = "infercnv")
  gene_order <- data.frame(
    gene = rownames(infercnv_genes_example),
    chr = infercnv_genes_example$V2,
    start = infercnv_genes_example$V3,
    end = infercnv_genes_example$V4,
    stringsAsFactors = FALSE
  )
  srt <- Seurat::CreateSeuratObject(counts = as.matrix(infercnv_data_example))
  srt$celltype <- infercnv_annots_example$V2
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  list(srt = srt, gene_order = gene_order)
}

test_that("RunCNV copykat matches the original copykat pipeline", {
  skip_on_cran()
  skip_if_not_installed("copykat")

  srt <- cnv_real_breast_seurat()
  wrapped <- RunCNV(srt, method = "copykat", genome = "hg38", verbose = FALSE)
  wrapped_pred <- wrapped$CNV_prediction

  suppressPackageStartupMessages(library(copykat))
  oldwd <- getwd()
  setwd(tempdir())
  on.exit(setwd(oldwd), add = TRUE)
  set.seed(1)
  direct <- copykat::copykat(
    rawmat = as.matrix(copykat::exp.rawdata),
    id.type = "S",
    sam.name = "scop_consistency",
    genome = "hg20",
    n.cores = 1,
    plot.genes = FALSE,
    output.seg = FALSE
  )
  direct_pred <- stats::setNames(direct$prediction$copykat.pred, direct$prediction$cell.names)

  common <- intersect(names(wrapped_pred), names(direct_pred))
  expect_gt(length(common), 100)
  expect_equal(
    mean(as.character(wrapped_pred[common]) == as.character(direct_pred[common])),
    1
  )
  expect_true(all(c("aneuploid", "diploid") %in% wrapped_pred))
  mat <- wrapped@tools$CNV$methods$copykat$cnv_matrix
  expect_identical(colnames(mat), colnames(srt))
  expect_gt(nrow(mat), 1000)
})

test_that("RunCNV scevan matches the original SCEVAN pipeline", {
  skip_on_cran()
  skip_if_not_installed("SCEVAN")

  srt <- cnv_real_breast_seurat()
  wrapped <- RunCNV(srt, method = "scevan", genome = "hg38", verbose = FALSE)
  wrapped_pred <- wrapped$CNV_prediction

  oldwd <- getwd()
  setwd(tempdir())
  on.exit(setwd(oldwd), add = TRUE)
  set.seed(1)
  direct <- SCEVAN::pipelineCNA(
    count_mtx = as.matrix(copykat::exp.rawdata),
    sample = "scop_cnv",
    par_cores = 1,
    SUBCLONES = TRUE,
    plotTree = FALSE,
    organism = "human"
  )
  direct_class <- direct$class
  names(direct_class) <- colnames(copykat::exp.rawdata)

  expect_equal(
    mean(as.character(wrapped_pred) == as.character(direct_class)),
    1
  )
  expect_true(all(c("normal", "tumor") %in% wrapped_pred))
  mat <- wrapped@tools$CNV$methods$scevan$cnv_matrix
  expect_gt(nrow(mat), 1000)
  expect_identical(colnames(mat), colnames(srt))
})

test_that("RunCNV infercnv matches the original infercnv pipeline", {
  skip_on_cran()
  skip_if_not_installed("infercnv")

  inputs <- cnv_real_infercnv_inputs()
  srt <- inputs$srt
  gene_order <- inputs$gene_order

  wrapped <- RunCNV(
    srt,
    method = "infercnv",
    reference.by = "celltype",
    reference = "normal",
    gene_order = gene_order,
    verbose = FALSE
  )
  wrapped_mat <- wrapped@tools$CNV$methods$infercnv$cnv_matrix

  set.seed(2024)
  out_dir <- tempfile("infercnv_cons_")
  dir.create(out_dir)
  ann_file <- tempfile("ann_", fileext = ".txt")
  gene_file <- tempfile("gene_", fileext = ".txt")
  annotations <- data.frame(
    cell = colnames(srt),
    group = srt$celltype,
    stringsAsFactors = FALSE
  )
  utils::write.table(annotations, ann_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(gene_order, gene_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = as.matrix(infercnv::infercnv_data_example),
    annotations_file = ann_file,
    delim = "\t",
    gene_order_file = gene_file,
    ref_group_names = "normal"
  )
  direct <- infercnv::run(obj,
    cutoff = 0.1,
    out_dir = out_dir,
    cluster_by_groups = FALSE,
    cluster_references = FALSE,
    denoise = FALSE,
    HMM = FALSE,
    analysis_mode = "samples",
    tumor_subcluster_partition_method = "qnorm",
    plot_steps = FALSE,
    resume_mode = FALSE,
    no_plot = TRUE,
    no_prelim_plot = TRUE,
    save_rds = FALSE,
    save_final_rds = FALSE,
    plot_probabilities = FALSE,
    write_expr_matrix = FALSE,
    write_phylo = FALSE,
    inspect_subclusters = FALSE,
    num_threads = 1
  )
  direct_expr <- direct@expr.data

  common_bins <- intersect(rownames(wrapped_mat), rownames(direct_expr))
  common_cells <- intersect(colnames(wrapped_mat), colnames(direct_expr))
  expect_gt(length(common_bins), 1000)
  rho <- suppressWarnings(stats::cor(
    as.numeric(wrapped_mat[common_bins, common_cells]),
    as.numeric(direct_expr[common_bins, common_cells]),
    method = "spearman",
    use = "complete.obs"
  ))
  expect_gt(rho, 0.999)
})

test_that("RunCNV numbat matches the original numbat pipeline output tables", {
  skip_on_cran()
  skip_if_not_installed("numbat")

  cm <- numbat::count_mat_example
  srt <- Seurat::CreateSeuratObject(counts = cm)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  out_dir <- tempfile("numbat_cons_")
  wrapped <- RunCNV(
    srt,
    method = "numbat",
    allele_counts = numbat::df_allele_example,
    reference_counts = numbat::ref_hca,
    genome = "hg38",
    output_dir = out_dir,
    verbose = FALSE
  )
  wrapped_pred <- wrapped$CNV_prediction

  final_clone <- grep("^clone_post_2[.]tsv$", list.files(out_dir), value = TRUE)
  expect_length(final_clone, 1)
  clone_post <- utils::read.delim(file.path(out_dir, final_clone))
  direct_pred <- stats::setNames(as.character(clone_post$compartment_opt), clone_post$cell)
  common <- intersect(names(wrapped_pred), names(direct_pred))
  expect_gt(length(common), 100)
  expect_equal(
    mean(as.character(wrapped_pred[common]) == direct_pred[common]),
    1
  )
  expect_true(all(c("normal", "tumor") %in% wrapped_pred))
})

test_that("RunCNV fastCNV matches the original fastCNV pipeline", {
  skip_on_cran()
  skip_if_not_installed("fastCNV")

  srt <- cnv_real_breast_seurat()
  srt$celltype <- ifelse(seq_len(ncol(srt)) %% 5 == 0, "Normal", "Tumor")
  wrapped <- RunCNV(
    srt,
    method = "fastCNV",
    reference.by = "celltype",
    reference = "Normal",
    genome = "hg38",
    verbose = FALSE
  )
  wrapped_score <- wrapped$CNV_fastCNV_score

  suppressPackageStartupMessages(library(fastCNV))
  direct <- fastCNV::fastCNV(
    seuratObj = srt,
    sampleName = "scop_cnv",
    referenceVar = "celltype",
    referenceLabel = "Normal",
    assay = "RNA",
    prepareCounts = FALSE,
    doPlot = FALSE,
    getCNVClusters = FALSE,
    getCNVPerChromosomeArm = FALSE,
    savePath = NULL
  )
  direct_score <- direct$cnv_fraction

  common <- intersect(names(wrapped_score), names(direct_score))
  expect_gt(length(common), 100)
  rho <- suppressWarnings(stats::cor(
    wrapped_score[common],
    direct_score[common],
    method = "spearman",
    use = "complete.obs"
  ))
  expect_gt(rho, 0.999)
  expect_gt(nrow(wrapped@tools$CNV$methods$fastCNV$cnv_matrix), 50)
})
