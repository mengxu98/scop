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

test_that("RunCNV stores unified schema and aligned metadata", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_backend = function(method, counts, ...) {
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
    cnv_run_backend = function(...) {
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
    cnv_run_backend = function(...) stop("backend should not run")
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
    RunCNV(
      srt,
      method = "casper",
      reference.by = "celltype",
      reference = "Normal",
      gene_order = gene_order,
      verbose = FALSE
    ),
    "loh"
  )
  expect_error(
    RunCNV(
      srt,
      method = "casper",
      reference.by = "celltype",
      reference = "Normal",
      gene_order = gene_order,
      loh = list(S1 = matrix(0, nrow = 1, ncol = 1)),
      verbose = FALSE
    ),
    "cytoband"
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
    cnv_run_backend = function(method, ...) {
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
  RunCNV(
    srt,
    method = "CaSpER",
    reference.by = "celltype",
    reference = "Normal",
    gene_order = gene_order,
    loh = list(S1 = matrix(0, nrow = 1, ncol = 1)),
    cytoband = data.frame(chr = "chr1", start = 1, end = 100),
    verbose = FALSE
  )

  expect_equal(seen, c("fastCNV", "scevan", "numbat", "casper"))
})

test_that("RunCNV keeps method-specific metadata across multiple methods", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_backend = function(method, ...) {
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
    cnv_run_backend = function(method, gene_order, ...) {
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

  for (method in c("copykat", "fastCNV", "scevan", "infercnv", "numbat", "casper")) {
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

test_that("fastCNV extraction reads Seurat assays and derives cell predictions", {
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
  fast_obj$CNV_fastCNV_score <- c(0.2, 0.1, 0.6)
  fast_obj$CNV_fastCNV_prediction <- NA_character_
  fast_obj$cnv_clusters <- c("A", "A", "B")

  extracted <- cnv_extract_fastcnv(
    fast_obj,
    cells = colnames(srt),
    reference.by = "celltype",
    reference = "Normal"
  )

  expect_equal(dim(extracted$cnv_matrix), c(2, 3))
  expect_equal(extracted$cell_info$score, c(0.2, 0.1, 0.6))
  expect_equal(extracted$cell_info$prediction, c("aneuploid", "diploid", "aneuploid"))
  expect_equal(extracted$cell_info$cluster, c("A", "A", "B"))
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
    get_namespace_fun = function(package, name) {
      expect_identical(package, "fastCNV")
      expect_identical(name, "fastCNV")
      function(
        seuratObj,
        sampleName,
        referenceVar,
        referenceLabel,
        assay,
        prepareCounts,
        doPlot,
        savePath
      ) {
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

test_that("CaSpER extraction reads chromosome-arm event matrices", {
  cells <- paste0("Cell", 1:3)
  final_chr_mat <- matrix(
    c(1, 0, -1, 1, 0, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(cells, c("1p", "1q"))
  )

  extracted <- cnv_extract_casper(
    result = list(large_scale_events = final_chr_mat),
    cells = cells
  )

  expect_equal(dim(extracted$cnv_matrix), c(2, 3))
  expect_equal(extracted$cnv_matrix["1p", "Cell1"], 1)
  expect_equal(extracted$cnv_matrix["1p", "Cell2"], -1)
  expect_equal(extracted$bin_info$chr, c("chr1", "chr1"))
  expect_equal(extracted$cell_info$prediction, c("aneuploid", "aneuploid", "aneuploid"))
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

test_that("CNVPlot returns bar, heatmap, and tree objects from stored schema", {
  srt <- make_cnv_seurat()
  testthat::local_mocked_bindings(
    cnv_run_backend = function(...) mock_cnv_backend()
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
