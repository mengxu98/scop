test_that("LINGER metadata follows selected cell identities", {
  counts <- matrix(c(2L, 3L, 4L, 5L, 6L, 7L),
    nrow = 2L,
    dimnames = list(c("G1", "G2"), c("c1", "c2", "c3"))
  )
  object <- SeuratObject::CreateSeuratObject(methods::as(counts, "dgCMatrix"))
  object$sample <- c("s1", "s1", "s2")
  object$label <- c("B", "T", "Mono")
  input <- list(object = object, cells = c("c3", "c1"), group = "all")
  config <- list(method_args = list(sample.by = "sample", label.by = "label"))
  metadata <- grn_linger_metadata(input, config)
  expect_identical(metadata$cell_id, c("c3", "c1"))
  expect_identical(metadata$sample, c("s2", "s1"))
  expect_identical(metadata$label, c("Mono", "B"))
  config$method_args$sample.by <- "missing"
  expect_error(grn_linger_metadata(input, config), "existing metadata columns")
  config$method_args <- list()
  metadata <- grn_linger_metadata(input, config)
  expect_identical(metadata$sample, rep("sample1", 2L))
  expect_identical(metadata$label, rep("all", 2L))
  input$group <- "B"
  config$group.by <- c(c1 = "B", c2 = "T", c3 = "B")
  expect_identical(grn_linger_metadata(input, config)$label, c("B", "B"))
  config$method_args$label.by <- "label"
  expect_identical(grn_linger_metadata(input, config)$label, c("Mono", "B"))
  input$cells <- c("c1", "missing")
  expect_error(grn_linger_metadata(input, config), "missing from Seurat metadata")
})

test_that("LINGER checks counts at the modality boundary", {
  counts <- matrix(1:6,
    nrow = 2L,
    dimnames = list(c("01", "1"), c("c1", "c2", "c3"))
  )
  object <- SeuratObject::CreateSeuratObject(methods::as(counts, "dgCMatrix"))
  actual <- grn_linger_counts(object, "RNA", "counts", c("c3", "c1"))
  expect_identical(rownames(actual), c("01", "1"))
  expect_equal(as.matrix(actual), counts[, c("c3", "c1")])
  expect_error(grn_linger_counts(object, "RNA", "data", "c1"), "count layers")
  expect_error(grn_linger_counts(object, "ATAC", "counts", "c1"), "selected RNA and ATAC assays")
  expect_error(grn_linger_counts(object, "RNA", "counts", "missing"), "match the selected group")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(methods::as(counts, "dgCMatrix"))
  expect_equal(as.matrix(grn_linger_counts(object, "ATAC", "counts", c("c3", "c1"))), counts[, c("c3", "c1")])
  object <- SeuratObject::SetAssayData(object, assay = "RNA", layer = "counts", new.data = counts / 2)
  expect_error(grn_linger_counts(object, "RNA", "counts", "c1"), "integer RNA and ATAC counts")
  object[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = methods::as(counts, "dgCMatrix"))
  object[["RNA"]] <- split(object[["RNA"]], f = c("A", "B", "A"))
  expect_error(grn_linger_counts(object, "RNA", "counts", "c1"), "JoinLayers")
})

test_that("LINGER rejects unsupported inputs before starting a process", {
  expect_error(RunGRN(matrix(1, 2L, 2L), methods = "linger"), "Multiome inference requires a Seurat object")
  object <- SeuratObject::CreateSeuratObject(methods::as(matrix(1L, 2L, 3L,
    dimnames = list(c("G1", "G2"), c("c1", "c2", "c3"))
  ), "dgCMatrix"))
  input <- list(object = object, input_mode = "multiome", cells = colnames(object))
  expect_error(RunLinger(input, list(method_args = list(resume_dir = tempfile()))), "completed prepare checkpoint")
  expect_error(RunLinger(input, list(method_args = list(activation = "bad"))), "activation")
  expect_error(RunLinger(input, list(method_args = list(activation = NA_character_))), "activation")
  expect_error(RunLinger(input, list(method_args = list(patch_missing_keys = "yes"))), "TRUE or FALSE")
  expect_error(RunLinger(input, list(method_args = list(min_cells = -1))), "nonnegative integers")
  expect_error(RunLinger(input, list(method_args = list(), resources = list())), "bulk_reference")
})

test_that("LINGER preflight checks only the required reference files", {
  counts <- methods::as(matrix(1L, 2L, 3L,
    dimnames = list(c("G1", "G2"), c("c1", "c2", "c3"))
  ), "dgCMatrix")
  object <- SeuratObject::CreateSeuratObject(counts)
  peaks <- counts
  rownames(peaks) <- c("chr19-10-20", "chr19:30-40")
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(peaks)
  input <- list(object = object)
  directory <- tempfile("linger-reference-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  config <- list(atac_assay = "ATAC", method_args = list(genome = "hg38"), resources = list(bulk_reference = directory))
  selected <- grn_linger_reference_files(input, config)
  expect_identical(selected$files, c(
    "TFName.txt", "Match2.txt", "bulk_gene_all.txt", "all_hg19.txt", "motifWeight.txt",
    "RE_gene_corr_hg19.bed", "hg38_hg19_pair.bed", "MotifTarget_Matrix_chr19.txt",
    "MotifTarget_matrix_chr19.bed", "hg38_Peaks_chr19.bed", "chr19_gene.txt",
    "chr19_index.txt", "chr19_index_all.txt", "fisher_chr19.pt", "all_models_chr19.pt"
  ))
  expect_identical(selected$chromosomes, "chr19")
  expect_false("all_models_chr1.pt" %in% selected$files)
  for (file in selected$files) writeLines("A", file.path(directory, file))
  for (genome in list("mm10", c("hg19", "hg38"))) {
    expect_error(grn_linger_reference_files(input, utils::modifyList(
      config,
      list(method_args = list(genome = genome))
    )), "hg19.*hg38")
  }
  checked <- grn_linger_check_resources(input, config)
  expect_identical(checked$selection, selected)
  used <- file.path(directory, "all_models_chr19.pt")
  unlink(used)
  expect_error(grn_linger_check_resources(input, config), "all_models_chr19.pt")
  config$method_args$genome <- "hg19"
  expect_true("hg19_hg38_pair.bed" %in% grn_linger_reference_files(input, config)$files)
})
