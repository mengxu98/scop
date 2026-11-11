make_cnv_seurat <- function() {
  counts <- matrix(
    c(
      5, 4, 0,
      0, 2, 6,
      3, 1, 0,
      0, 1, 4
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$celltype <- c("Tumor", "Normal", "Tumor")
  srt$sample <- c("S1", "S1", "S2")
  emb <- matrix(
    c(0, 0, 1, 0, 0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(colnames(srt), c("UMAP_1", "UMAP_2"))
  )
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb,
    assay = "RNA",
    key = "UMAP_"
  )
  srt
}

mock_cnv_backend <- function(...) {
  list(
    cnv_matrix = matrix(
      c(
        -0.3, 0.2,
        0.4, -0.1,
        0.0, 0.8
      ),
      nrow = 3,
      dimnames = list(paste0("bin", 1:3), c("Cell2", "Cell1"))
    ),
    bin_info = data.frame(
      bin_id = paste0("bin", 1:3),
      chr = c("chr1", "chr1", "chr2"),
      start = c(1, 101, 1),
      end = c(100, 200, 100),
      gene = paste0("Gene", 1:3)
    ),
    cell_info = data.frame(
      cell = c("Cell2", "Cell1"),
      prediction = c("diploid", "aneuploid"),
      cluster = c("C0", "C1")
    )
  )
}

test_that("RunCNV stores a unified result bundle and aligned metadata", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_model = function(method, counts, ...) {
      expect_identical(method, "copykat")
      expect_s4_class(counts, "dgCMatrix")
      mock_cnv_backend()
    }
  )

  out <- RunCNV(srt, method = "copykat", verbose = FALSE)

  expect_true("CNV" %in% names(out@tools))
  expect_identical(out@tools$CNV$active_method, "copykat")
  bundle <- out@tools$CNV$methods$copykat
  expect_equal(dim(bundle$cnv_matrix), c(3, 3))
  expect_identical(colnames(bundle$cnv_matrix), colnames(out))
  expect_equal(bundle$bin_info$chr, c("chr1", "chr1", "chr2"))
  expect_equal(unname(out$CNV_prediction), c("aneuploid", "diploid", NA))
  expect_equal(unname(out$CNV_copykat_cluster), c("C1", "C0", NA))
  expect_true(is.na(out$CNV_score[3]))
})

test_that("RunCNV orients cell-by-bin backend matrices", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_model = function(...) {
      list(
        cnv_matrix = matrix(
          c(1, 2, 3, 4, 5, 6),
          nrow = 3,
          dimnames = list(colnames(srt), c("binA", "binB"))
        )
      )
    }
  )
  out <- RunCNV(srt, method = "scevan", verbose = FALSE)
  mat <- out@tools$CNV$methods$scevan$cnv_matrix
  expect_equal(dim(mat), c(2, 3))
  expect_identical(rownames(mat), c("binA", "binB"))
  expect_identical(colnames(mat), colnames(srt))
})

test_that("RunCNV validates reference and gene-order requirements before backend work", {
  srt <- make_cnv_seurat()
  gene_order <- data.frame(
    gene = rownames(srt),
    chr = "chr1",
    start = seq_len(nrow(srt)) * 100,
    end = seq_len(nrow(srt)) * 100 + 99
  )
  testthat::local_mocked_bindings(
    cnv_run_model = function(...) stop("backend should not run")
  )
  expect_error(
    RunCNV(srt, method = "missing_method", verbose = FALSE),
    "method"
  )
  expect_error(
    RunCNV(srt, method = "numbat", verbose = FALSE),
    "allele_counts"
  )
  expect_error(
    RunCNV(
      srt,
      method = "numbat",
      allele_counts = data.frame(cell = colnames(srt)),
      verbose = FALSE
    ),
    "reference_counts"
  )
  expect_error(
    RunCNV(srt, method = "fastCNV", verbose = FALSE),
    "reference.by"
  )
  expect_error(
    RunCNV(
      srt,
      method = "infercnv",
      reference.by = "celltype",
      reference = "Normal",
      verbose = FALSE
    ),
    "gene_order"
  )
  expect_error(
    RunCNV(
      srt,
      method = "fastCNV",
      reference.by = "celltype",
      reference = "Absent",
      verbose = FALSE
    ),
    "No reference cells"
  )
})

test_that("RunCNV accepts common method aliases", {
  srt <- make_cnv_seurat()
  gene_order <- data.frame(
    gene = rownames(srt),
    chr = "chr1",
    start = seq_len(nrow(srt)) * 100,
    end = seq_len(nrow(srt)) * 100 + 99
  )
  seen <- character()
  testthat::local_mocked_bindings(
    cnv_run_model = function(method, ...) {
      seen <<- c(seen, method)
      mock_cnv_backend()
    }
  )

  RunCNV(srt, method = "fastcnv", reference.by = "celltype", reference = "Normal", verbose = FALSE)
  RunCNV(srt, method = "SCEVAN", verbose = FALSE)
  RunCNV(
    srt,
    method = "Numbat",
    allele_counts = data.frame(cell = colnames(srt)),
    reference_counts = matrix(1, nrow = nrow(srt), ncol = 1, dimnames = list(rownames(srt), "ref")),
    verbose = FALSE
  )

  expect_equal(seen, c("fastCNV", "scevan", "numbat"))
})

test_that("RunCNV keeps method-specific metadata across multiple methods", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_model = function(method, ...) {
      out <- mock_cnv_backend()
      if (identical(method, "scevan")) {
        out$cell_info$prediction <- c("normal", "malignant")
        out$cell_info$cluster <- c("S0", "S1")
      }
      out
    }
  )

  out <- RunCNV(srt, method = "copykat", verbose = FALSE)
  out <- RunCNV(out, method = "scevan", verbose = FALSE)

  expect_identical(out@tools$CNV$active_method, "scevan")
  expect_true(all(c("copykat", "scevan") %in% names(out@tools$CNV$methods)))
  expect_equal(unname(out$CNV_copykat_prediction), c("aneuploid", "diploid", NA))
  expect_equal(unname(out$CNV_scevan_prediction), c("malignant", "normal", NA))
  expect_equal(unname(out$CNV_prediction), c("malignant", "normal", NA))
})

test_that("RunCNV supports gene-order files and store_matrix = FALSE", {
  srt <- make_cnv_seurat()
  gene_order <- data.frame(
    gene = rownames(srt),
    chr = "chr1",
    start = seq_len(nrow(srt)) * 100,
    end = seq_len(nrow(srt)) * 100 + 99
  )
  gene_file <- tempfile(fileext = ".tsv")
  utils::write.table(gene_order, gene_file, sep = "\t", quote = FALSE, row.names = FALSE)

  testthat::local_mocked_bindings(
    cnv_run_model = function(method, gene_order, ...) {
      expect_identical(method, "infercnv")
      expect_equal(colnames(gene_order), c("gene", "chr", "start", "end"))
      mock_cnv_backend()
    }
  )
  out <- RunCNV(
    srt,
    method = "infercnv",
    reference.by = "celltype",
    reference = "Normal",
    gene_order = gene_file,
    store_matrix = FALSE,
    verbose = FALSE
  )

  bundle <- out@tools$CNV$methods$infercnv
  expect_null(bundle$cnv_matrix)
  expect_equal(bundle$matrix_dim, c(3, 3))
  expect_equal(unname(out$CNV_infercnv_prediction), c("aneuploid", "diploid", NA))
  expect_error(CNVPlot(out, plot_type = "heatmap"), "store_matrix")
})

test_that("CNV standardization derives fallback predictions for every backend", {
  cells <- paste0("Cell", 1:3)
  mat <- matrix(
    c(
      0.3, 0.2, 0.0,
      0.4, 0.1, 0.0
    ),
    nrow = 2,
    dimnames = list(c("seg1", "seg2"), cells)
  )

  for (method in c("copykat", "fastCNV", "scevan", "infercnv", "numbat")) {
    bundle <- cnv_standardize_result(
      result = list(
        cnv_matrix = mat,
        cell_info = data.frame(
          cell = cells,
          prediction = NA_character_,
          stringsAsFactors = FALSE
        )
      ),
      method = method,
      cells = cells,
      features = rownames(mat),
      reference_cells = "Cell2"
    )

    expect_equal(
      unname(bundle$cell_info$prediction),
      c("aneuploid", "diploid", "diploid"),
      info = method
    )
  }
})

test_that("fastCNV extracts genomic scores even with DecontX UMAP present", {
  srt <- make_cnv_seurat()
  mat <- matrix(
    c(-0.2, 0.1, 0.3, 0.5, -0.4, 0.2),
    nrow = 2,
    dimnames = list(c("seg1", "seg2"), colnames(srt))
  )
  fast_obj <- srt
  suppressWarnings({
    fast_obj[["genomicScores"]] <- SeuratObject::CreateAssayObject(data = mat)
  })
  fast_obj@tools$decontX <- list(umap = matrix(
    c(8, 2, 3, 1, 4, 6),
    ncol = 2,
    dimnames = list(colnames(srt), c("DecontX_UMAP_1", "DecontX_UMAP_2"))
  ))
  fast_obj$CNV_fastCNV_score <- c(0.2, 0.1, 0.6)
  fast_obj$CNV_fastCNV_prediction <- NA_character_
  fast_obj$cnv_clusters <- c("A", "A", "B")

  extracted <- cnv_extract_fastcnv(
    fast_obj,
    cells = colnames(srt),
    reference.by = "celltype",
    reference = "Normal"
  )

  expect_equal(as.matrix(extracted$cnv_matrix), mat)
  expect_equal(extracted$cell_info$score, c(0.2, 0.1, 0.6))
  expect_equal(extracted$cell_info$prediction, c("aneuploid", "diploid", "aneuploid"))
  expect_equal(extracted$cell_info$cluster, c("A", "A", "B"))
})

test_that("fastCNV reports a missing genomicScores assay instead of using a tool matrix", {
  srt <- make_cnv_seurat()
  srt@tools$decontX <- list(umap = SeuratObject::Embeddings(srt[["umap"]]))

  expect_error(cnv_extract_fastcnv(srt, cells = colnames(srt)), "genomicScores")
})

test_that("fastCNV runner skips expensive upstream preparation by default", {
  srt <- make_cnv_seurat()
  mat <- matrix(
    c(-0.2, 0.1, 0.3, 0.5, -0.4, 0.2),
    nrow = 2,
    dimnames = list(c("seg1", "seg2"), colnames(srt))
  )
  fast_obj <- srt
  suppressWarnings({
    fast_obj[["genomicScores"]] <- SeuratObject::CreateAssayObject(data = mat)
  })
  fast_obj$cnv_fraction <- c(0.2, 0.1, 0.6)
  out_dir <- tempfile("fastcnv_")

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    cnv_attach_package_data = function(package) {
      expect_identical(package, "fastCNV")
      FALSE
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "fastCNV")
      expect_identical(name, "fastCNV")
      function(seuratObj,
               sampleName,
               referenceVar,
               referenceLabel,
               assay,
               prepareCounts,
               doPlot,
               savePath) {
        expect_s4_class(seuratObj, "Seurat")
        expect_identical(sampleName, "scop_cnv")
        expect_identical(referenceVar, "celltype")
        expect_identical(referenceLabel, "Normal")
        expect_identical(assay, "RNA")
        expect_false(prepareCounts)
        expect_false(doPlot)
        expect_identical(savePath, out_dir)
        fast_obj
      }
    }
  )

  extracted <- cnv_run_fastcnv(
    srt = srt,
    assay = "RNA",
    layer = "counts",
    reference.by = "celltype",
    reference = "Normal",
    genome = "hg38",
    output_dir = out_dir,
    verbose = FALSE
  )

  expect_equal(extracted$cell_info$prediction, c("aneuploid", "diploid", "aneuploid"))
})

test_that("RunCNV saves infercnv coordinates in matrix row order", {
  srt <- make_cnv_seurat()
  gene_order <- data.frame(
    gene = c("Gene4", "Gene2", "Gene3", "Gene1"),
    chr = c("chr2", "chr1", "chr2", "chr1"),
    start = c(400, 200, 300, 100),
    end = c(499, 299, 399, 199)
  )
  mat <- matrix(
    c(0.99, 1.2, 0.81, 0.96, 1.1, 0.85, 1.01, 1.3, 0.86),
    nrow = 3,
    dimnames = list(c("Gene3", "Gene1", "unmatched"), colnames(srt))
  )
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(package, name) {
      if (identical(name, "CreateInfercnvObject")) {
        return(function(...) list())
      }
      function(...) mat
    }
  )

  out <- RunCNV(
    srt,
    method = "infercnv",
    reference.by = "celltype",
    reference = "Normal",
    gene_order = gene_order,
    verbose = FALSE
  )
  bundle <- out@tools$CNV$methods$infercnv

  expect_identical(bundle$cnv_matrix, mat)
  expect_identical(bundle$bin_info$gene, rownames(mat))
  expect_identical(bundle$bin_info$chr, c("chr2", "chr1", NA_character_))
  expect_identical(bundle$bin_info$start, c(300, 100, NA_real_))
  expect_identical(bundle$bin_info$end, c(399, 199, NA_real_))
})

test_that("infercnv runner avoids subclusters by default and filters passthrough arguments", {
  srt <- make_cnv_seurat()
  counts <- cnv_get_counts(srt, assay = "RNA", layer = "counts")
  gene_order <- data.frame(
    gene = rownames(srt),
    chr = "chr1",
    start = seq_len(nrow(srt)) * 100,
    end = seq_len(nrow(srt)) * 100 + 99
  )
  out_dir <- tempfile("infercnv_")

  testthat::local_mocked_bindings(
    check_r = function(package, ...) {
      expect_identical(package, "broadinstitute/infercnv")
      TRUE
    },
    get_namespace_fun = function(package, name) {
      expect_identical(package, "infercnv")
      if (identical(name, "CreateInfercnvObject")) {
        return(function(raw_counts_matrix,
                        annotations_file,
                        delim,
                        gene_order_file,
                        ref_group_names) {
          expect_true(is.matrix(raw_counts_matrix))
          expect_equal(delim, "\t")
          expect_equal(ref_group_names, "reference")
          annotations <- utils::read.delim(annotations_file, header = FALSE)
          expect_equal(annotations[[2]], c("observation", "reference", "observation"))
          genes <- utils::read.delim(gene_order_file, header = FALSE)
          expect_equal(nrow(genes), nrow(gene_order))
          list(cells = colnames(raw_counts_matrix))
        })
      }
      if (identical(name, "run")) {
        return(function(infercnv_obj,
                        cutoff,
                        out_dir,
                        cluster_by_groups,
                        cluster_references,
                        denoise,
                        HMM,
                        analysis_mode,
                        tumor_subcluster_partition_method,
                        plot_steps,
                        resume_mode,
                        no_plot,
                        no_prelim_plot,
                        save_rds,
                        save_final_rds,
                        plot_probabilities,
                        write_expr_matrix,
                        write_phylo,
                        inspect_subclusters,
                        num_threads) {
          expect_equal(cutoff, 0.1)
          expect_true(dir.exists(out_dir))
          expect_false(cluster_by_groups)
          expect_false(cluster_references)
          expect_false(denoise)
          expect_false(HMM)
          expect_identical(analysis_mode, "samples")
          expect_identical(tumor_subcluster_partition_method, "qnorm")
          expect_false(plot_steps)
          expect_false(resume_mode)
          expect_true(no_plot)
          expect_true(no_prelim_plot)
          expect_false(save_rds)
          expect_false(save_final_rds)
          expect_false(plot_probabilities)
          expect_false(write_expr_matrix)
          expect_false(write_phylo)
          expect_false(inspect_subclusters)
          expect_true(is.numeric(num_threads))
          expect_true(num_threads >= 1)
          matrix(
            c(0.2, 0.1, 0, 0, 0.4, 0.3),
            nrow = 2,
            dimnames = list(c("seg1", "seg2"), infercnv_obj$cells)
          )
        })
      }
      stop(name)
    }
  )

  extracted <- cnv_run_infercnv(
    counts = counts,
    reference.by = "celltype",
    reference = "Normal",
    reference_cells = "Cell2",
    gene_order = gene_order,
    output_dir = out_dir,
    verbose = FALSE,
    layer = "data",
    assay = "RNA"
  )

  expect_equal(dim(extracted$cnv_matrix), c(2, 3))
  expect_identical(colnames(extracted$cnv_matrix), colnames(srt))
})

test_that("backend argument filtering drops unsupported named arguments", {
  out <- cnv_call_fun(
    function(a, b, ...) {
      list(a = a, b = b, passthrough = list(...))
    },
    list(a = 1, b = 2, layer = "counts", assay = "RNA")
  )

  expect_equal(out$a, 1)
  expect_equal(out$b, 2)
  expect_equal(out$passthrough, list())
})

test_that("SCEVAN extraction reads CNA RData output and cell assignments", {
  srt <- make_cnv_seurat()
  out_dir <- tempfile("scevan_")
  dir.create(file.path(out_dir, "output"), recursive = TRUE)
  CNA_mtx_relat <- matrix(
    c(-0.3, 0.2, 0.4, 0.1, -0.2, 0.5),
    nrow = 2,
    dimnames = list(c("seg1", "seg2"), colnames(srt))
  )
  save(CNA_mtx_relat, file = file.path(out_dir, "output", "scop_cnv_CNAmtx.RData"))
  count_mtx_annot <- data.frame(
    seqnames = c("chr1", "chr2"),
    start = c(1, 101),
    end = c(100, 200),
    gene = c("Gene1", "Gene2"),
    row.names = rownames(CNA_mtx_relat)
  )
  save(count_mtx_annot, file = file.path(out_dir, "output", "scop_cnv_count_mtx_annot.RData"))
  backend_result <- data.frame(
    cell.assignment = c("malignant", "normal", "malignant"),
    subclone = c("S1", NA, "S2"),
    row.names = colnames(srt)
  )

  extracted <- cnv_extract_scevan(
    result = backend_result,
    cells = colnames(srt),
    sample = "scop_cnv",
    output_dir = out_dir
  )

  expect_equal(dim(extracted$cnv_matrix), c(2, 3))
  expect_equal(extracted$bin_info$chr, c("chr1", "chr2"))
  expect_equal(extracted$cell_info$prediction, c("malignant", "normal", "malignant"))
  expect_equal(extracted$cell_info$cluster, c("S1", NA, "S2"))
})

test_that("Numbat extraction reads posterior and clone tables", {
  cells <- paste0("Cell", 1:3)
  out_dir <- tempfile("numbat_")
  dir.create(out_dir)
  joint_post <- data.frame(
    cell = rep(cells, each = 2),
    seg = rep(c("seg1", "seg2"), times = 3),
    cnv_state_post = c("amp", "neu", "del", "neu", "amp", "del"),
    p_cnv = c(0.8, 0, 0.6, 0, 0.4, 0.9),
    CHROM = rep(c("1", "2"), times = 3),
    seg_start = rep(c(1, 101), times = 3),
    seg_end = rep(c(100, 200), times = 3)
  )
  clone_post <- data.frame(
    cell = cells,
    compartment_opt = c("tumor", "normal", "tumor"),
    clone_opt = c("clone1", "normal", "clone2"),
    p_cnv = c(0.8, 0.1, 0.9)
  )
  utils::write.table(joint_post, file.path(out_dir, "joint_post_2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(clone_post, file.path(out_dir, "clone_post_2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  extracted <- cnv_extract_numbat(result = 0, cells = cells, output_dir = out_dir)

  expect_equal(dim(extracted$cnv_matrix), c(2, 3))
  expect_equal(extracted$cnv_matrix["seg1", "Cell1"], 0.8)
  expect_equal(extracted$cnv_matrix["seg1", "Cell2"], -0.6)
  expect_equal(extracted$bin_info$chr, c("1", "2"))
  expect_equal(extracted$cell_info$prediction, c("tumor", "normal", "tumor"))
  expect_equal(extracted$cell_info$cluster, c("clone1", "normal", "clone2"))
})

test_that("fastCNV rejects unsupported mouse genome before backend execution", {
  srt <- make_cnv_seurat()
  expect_error(
    cnv_run_fastcnv(
      srt = srt,
      assay = "RNA",
      layer = "counts",
      reference.by = "celltype",
      reference = "Normal",
      genome = "mm10"
    ),
    "human data only"
  )
})

test_that("CNVPlot returns bar, heatmap, and tree objects from stored results", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_model = function(...) mock_cnv_backend()
  )
  out <- RunCNV(srt, method = "copykat", verbose = FALSE)

  p <- CNVPlot(out, plot_type = "bar")
  expect_s3_class(p, "ggplot")

  p_dim <- CNVPlot(
    out,
    plot_type = "dim",
    reduction = "umap",
    value = "CNV_prediction",
    theme_use = "theme_blank"
  )
  expect_true(inherits(p_dim, "ggplot") || inherits(p_dim, "patchwork"))

  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  ht <- CNVPlot(out, plot_type = "heatmap", group.by = "celltype")
  expect_s4_class(ht, "Heatmap")

  tree <- CNVPlot(out, plot_type = "tree", cluster_tree_by = "cell")
  expect_s3_class(tree, "hclust")
})

test_that("CNVPlot centers infercnv colors on one and preserves the original matrix", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  mat <- matrix(
    c(0.6, 0.75, 1, 1.25, 1.5, 2), nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), colnames(srt))
  )
  srt@tools$CNV <- list(
    active_method = "infercnv",
    methods = list(infercnv = list(method = "infercnv", cnv_matrix = mat))
  )
  colors <- c("#2166AC", "#FFFFFF", "#B2182B")

  ht <- CNVPlot(srt, heatmap_palcolor = colors, cluster_columns = FALSE)
  expect_equal(ht@matrix, t(mat))
  mapped <- ht@matrix_color_mapping@col_fun(c(0.5, 0.75, 1, 1.25, 1.5))
  expect_equal(mapped[c(1, 2, 4, 5)], c("#2166ACFF", "#2166ACFF", "#B2182BFF", "#B2182BFF"))
  expect_true(all(grDevices::col2rgb(mapped[3]) >= 245))
  expect_equal(ht@matrix_legend_param$at, c(0.75, 1, 1.25))

  custom <- CNVPlot(
    srt, heatmap_limits = c(0.9, 1.4), heatmap_palcolor = colors,
    cluster_columns = FALSE
  )
  expect_equal(custom@matrix, t(mat))
  expect_equal(
    custom@matrix_color_mapping@col_fun(c(0.8, 1, 1.5)),
    c("#2166ACFF", mapped[3], "#B2182BFF")
  )
  expect_equal(custom@matrix_legend_param$at, c(0.9, 1, 1.4))

  two_colors <- CNVPlot(srt, heatmap_limits = c(0.9, 1.4), heatmap_palcolor = colors[c(1, 3)])
  two_colors_symmetric <- CNVPlot(srt, heatmap_palcolor = colors[c(1, 3)])
  expect_equal(
    two_colors@matrix_color_mapping@col_fun(1),
    two_colors_symmetric@matrix_color_mapping@col_fun(1)
  )
  expect_error(CNVPlot(srt, heatmap_limits = c(0.6, 0.8)), "heatmap_limits")
  expect_error(CNVPlot(srt, heatmap_limits = c(NA, 1.4)), "heatmap_limits")
})

test_that("CNVPlot retains zero-centered default colors for other methods", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  mat <- matrix(
    c(-2, -0.5, 0, 0.5, 1, 2), nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), colnames(srt))
  )
  for (method in c("copykat", "fastCNV", "scevan", "numbat")) {
    srt@tools$CNV <- list(
      active_method = method,
      methods = setNames(list(list(method = method, cnv_matrix = mat)), method)
    )
    ht <- CNVPlot(srt, heatmap_palcolor = c("#2166AC", "#FFFFFF", "#B2182B"))
    expect_equal(ht@matrix, t(mat))
    mapped <- ht@matrix_color_mapping@col_fun(c(-2, -1, 0, 1, 2))
    expect_equal(mapped[c(1, 2, 4, 5)], c("#2166ACFF", "#2166ACFF", "#B2182BFF", "#B2182BFF"))
    expect_true(all(grDevices::col2rgb(mapped[3]) >= 245))
    expect_equal(ht@matrix_legend_param$at, c(-1, 0, 1))
  }
})

test_that("CNVPlot orders chromosomes horizontally and cells by group", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  srt$group <- factor(c("Tumor", "Reference", "Tumor"), levels = c("Reference", "Tumor"))
  mat <- matrix(
    c(1.4, 0.7, 0.8, 1.3, 1, 0.65, 1.05, 1.2, 1.1), nrow = 3,
    dimnames = list(c("gene3", "gene2", "gene1"), colnames(srt))
  )
  bins <- data.frame(chr = c("chr2", "chr1", "chr1"), start = c(30, 20, 10))
  srt@tools$CNV <- list(
    active_method = "infercnv",
    methods = list(infercnv = list(method = "infercnv", cnv_matrix = mat, bin_info = bins))
  )

  paper <- CNVPlot(srt, group.by = "group", title = "", use_raster = FALSE)
  expect_equal(paper@matrix, t(mat[c(3, 2, 1), c(2, 1, 3)]))
  expect_s4_class(paper@top_annotation, "HeatmapAnnotation")
  expect_s4_class(paper@left_annotation, "HeatmapAnnotation")
  expect_equal(paper@matrix_legend_param$at, c(0.75, 1, 1.25))
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  rendered <- tryCatch(ComplexHeatmap::draw(paper), finally = grDevices::dev.off())
  expect_equal(names(ComplexHeatmap::column_order(rendered)), c("chr1", "chr2"))
  expect_null(rendered@ht_list[[1]]@column_title)
  expect_equal(unname(unlist(ComplexHeatmap::row_order(rendered))), 1:3)
  expect_equal(srt@tools$CNV$methods$infercnv$cnv_matrix, mat)
})

test_that("CNVPlot automatic display limits use the native midpoint and retain scores", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  mat <- matrix(
    c(0.88, 0.97, 1.04, 1.12, 1.2, 2), nrow = 2,
    dimnames = list(c("gene1", "gene2"), colnames(srt))
  )
  srt@tools$CNV <- list(
    active_method = "infercnv",
    methods = list(infercnv = list(method = "infercnv", cnv_matrix = mat))
  )

  span <- as.double(stats::quantile(abs(mat - 1), 0.95))
  auto <- CNVPlot(srt, heatmap_limits = "auto")
  expect_equal(auto@matrix_legend_param$at, 1 + c(-span, 0, span))
  expect_equal(auto@matrix, t(mat))
  manual <- CNVPlot(srt, heatmap_limits = c(0.7, 1.8))
  expect_equal(manual@matrix_legend_param$at, c(0.7, 1, 1.8))
  expect_equal(manual@matrix, auto@matrix)

  srt@tools$CNV$methods$infercnv$cnv_matrix[,] <- 1
  constant <- CNVPlot(srt, heatmap_limits = "auto")
  expect_equal(constant@matrix_legend_param$at, c(0.75, 1, 1.25))
  expect_error(CNVPlot(srt, heatmap_limits = "unknown"), "heatmap_limits")
})

test_that("CNVPlot labels native fastCNV windows with missing coordinates", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  mat <- matrix(
    seq_len(9) / 10, nrow = 3,
    dimnames = list(c("2.p1", "1.p1", "1.q1"), colnames(srt))
  )
  srt@tools$CNV <- list(
    active_method = "fastCNV",
    methods = list(fastCNV = list(
      method = "fastCNV", cnv_matrix = mat,
      bin_info = data.frame(chr = rep(NA_character_, nrow(mat)), start = NA_real_)
    ))
  )

  paper <- CNVPlot(srt)
  expect_equal(paper@matrix, t(mat[c(2, 3, 1), ]))
  expect_s4_class(paper@top_annotation, "HeatmapAnnotation")
  expect_true(all(is.na(srt@tools$CNV$methods$fastCNV$bin_info$chr)))
})

test_that("CNVPlot derives labels for custom legend ticks", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  mat <- matrix(
    c(0.8, 1, 1.1, 1.2, 1.3, 0.9), nrow = 2,
    dimnames = list(c("gene1", "gene2"), colnames(srt))
  )
  srt@tools$CNV <- list(
    active_method = "infercnv",
    methods = list(infercnv = list(method = "infercnv", cnv_matrix = mat))
  )
  at <- c(0.75, 1, 1.1, 1.25)
  ht <- CNVPlot(srt, heatmap_legend_param = list(at = at), use_raster = FALSE)
  expect_equal(as.numeric(ht@matrix_legend_param$labels), at)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  rendered <- tryCatch(ComplexHeatmap::draw(ht), finally = grDevices::dev.off())
  expect_s4_class(rendered, "HeatmapList")
  labels <- c("Loss", "Neutral", "Gain", "High")
  custom <- CNVPlot(srt, heatmap_legend_param = list(at = at, labels = labels))
  expect_identical(custom@matrix_legend_param$labels, labels)
})

test_that("CNVPlot keeps chromosome labels aligned when slices are clustered", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  srt <- make_cnv_seurat()
  mat <- outer(
    rep(c(0.1, 1.1, 0.11, 3.1), each = 2) + rep(c(0, 0.02), 4),
    c(0, 0.2, 0.4), "+"
  )
  dimnames(mat) <- list(paste0("gene", 1:8), colnames(srt))
  srt@tools$CNV <- list(
    active_method = "infercnv",
    methods = list(infercnv = list(
      method = "infercnv", cnv_matrix = mat,
      bin_info = data.frame(chr = rep(paste0("chr", 1:4), each = 2), start = rep(1:2, 4))
    ))
  )
  ht <- CNVPlot(srt, cluster_rows = TRUE, cluster_column_slices = TRUE, use_raster = FALSE)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  output <- tryCatch({
    rendered <- ComplexHeatmap::draw(ht)
    grob_names <- grid::grid.ls(print = FALSE)$name
    text <- unlist(lapply(grob_names[grepl("^GRID.text", grob_names)], function(name) {
      as.character(grid::grid.get(name)$label)
    }))
    list(
      order = names(ComplexHeatmap::column_order(rendered)),
      labels = text[text %in% as.character(1:4)]
    )
  }, finally = grDevices::dev.off())
  expect_false(identical(output$order, paste0("chr", 1:4)))
  expect_equal(output$labels, sub("^chr", "", output$order))
})
