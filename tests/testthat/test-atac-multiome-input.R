make_multiome_input_fixture <- function() {
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  counts <- Matrix::Matrix(matrix(c(
    1, 0, 2, 5, 0, 3, 1, 0, 2, 6, 7, 1,
    0, 3, 4, 2, 1, 8
  ), nrow = 6L), sparse = TRUE)
  rownames(counts) <- c("chr1:0-50", "chr1:100-150", "chr2:200-250", "ENSG1", "ENSG2", "ENSG3")
  colnames(counts) <- c("AA-1", "BB-1", "CC-1")
  list(counts = counts, type = rep(c("Peaks", "Gene Expression"), each = 3L))
}

test_that("paired 10x H5 ingestion retains both actual native count modalities", {
  skip_if_not_installed("hdf5r")
  x <- make_multiome_input_fixture()
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  matrix <- h5$create_group("matrix")
  features <- matrix$create_group("features")
  matrix[["data"]] <- x$counts@x
  matrix[["indices"]] <- x$counts@i
  matrix[["indptr"]] <- x$counts@p
  matrix[["shape"]] <- dim(x$counts)
  matrix[["barcodes"]] <- colnames(x$counts)
  features[["id"]] <- rownames(x$counts)
  features[["name"]] <- c(rownames(x$counts)[1:3], "DUP", "DUP", "GENE3")
  features[["feature_type"]] <- x$type
  h5$close_all()
  native <- suppressMessages(Seurat::Read10X_h5(path, use.names = FALSE, unique.features = FALSE))
  object <- suppressMessages(ReadMultiomeData(path,
    genome = "hg38", sample = "capture1",
    cell.names = paste0("capture1_", colnames(x$counts)), verbose = FALSE
  ))
  expect_s4_class(object[["peaks"]], "ChromatinAssay")
  expect_s4_class(object[["RNA"]], "Assay")
  expect_equal(unname(SeuratObject::GetAssayData(object, assay = "peaks", layer = "counts")), unname(native$Peaks))
  expect_equal(unname(SeuratObject::GetAssayData(object, assay = "RNA", layer = "counts")), unname(native[["Gene Expression"]]))
  expect_identical(rownames(object[["RNA"]]), c("ENSG1", "ENSG2", "ENSG3"))
  expect_identical(as.vector(colnames(object[["RNA"]])), as.vector(colnames(object[["peaks"]])))
  expect_identical(GenomicRanges::start(Signac::granges(object[["peaks"]])), c(1L, 101L, 201L))
  expect_identical(SeuratObject::DefaultAssay(object), "peaks")
  expect_identical(SeuratObject::Misc(object[["peaks"]], "scop_atac_input")$input_coordinate_system, "0-based half-open")
  receipt <- SeuratObject::Misc(object, "scop_multiome_input")
  expect_identical(receipt$expression_assay, "RNA")
  expect_identical(receipt$input_barcodes, colnames(x$counts))
  expect_identical(as.vector(object$sample), rep("capture1", 3))
})

test_that("paired compressed Matrix Market imports preserve feature IDs and barcodes", {
  x <- make_multiome_input_fixture()
  directory <- tempfile("multiome_mtx_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  Matrix::writeMM(x$counts, file.path(directory, "matrix.mtx"))
  writeLines(colnames(x$counts), file.path(directory, "barcodes.tsv"))
  features <- data.frame(id = rownames(x$counts), name = rownames(x$counts), type = x$type)
  utils::write.table(features, file.path(directory, "features.tsv"),
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  for (file in c("matrix.mtx", "features.tsv", "barcodes.tsv")) {
    input <- readBin(file.path(directory, file), "raw", n = file.info(file.path(directory, file))$size)
    con <- gzfile(file.path(directory, paste0(file, ".gz")), "wb")
    writeBin(input, con)
    close(con)
  }
  object <- suppressMessages(ReadMultiomeData(directory, genome = "hg38", verbose = FALSE))
  expect_equal(unname(SeuratObject::GetAssayData(object, assay = "RNA", layer = "counts")), unname(x$counts[4:6, ]))
  expect_equal(unname(SeuratObject::GetAssayData(object, assay = "peaks", layer = "counts")), unname(x$counts[1:3, ]))
  expect_identical(as.vector(colnames(object)), colnames(x$counts))
  expect_identical(SeuratObject::Misc(object, "scop_multiome_input")$source, "10x_arc_matrix_market")
  expect_error(ReadMultiomeData(directory, rna.assay = "peaks"), "distinct")
  con <- gzfile(file.path(directory, "features.tsv.gz"), "wt")
  features$id[5] <- features$id[4]
  utils::write.table(features, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(con)
  expect_error(ReadMultiomeData(directory, verbose = FALSE), "unique feature IDs")
})

test_that("multiome input rejects absent or unpaired modalities rather than inventing RNA", {
  x <- make_multiome_input_fixture()
  path <- tempfile()
  writeLines("contract-test placeholder", path)
  on.exit(unlink(path), add = TRUE)
  testthat::local_mocked_bindings(atac_backend_fun = function(package, fun) identity, .package = "scop")
  testthat::local_mocked_bindings(Read10X_h5 = function(filename, use.names, unique.features) {
    list(Peaks = x$counts[1:3, ])
  }, .package = "Seurat")
  expect_error(ReadMultiomeData(path, verbose = FALSE), "one Gene Expression and one Peaks")
  testthat::local_mocked_bindings(Read10X_h5 = function(filename, use.names, unique.features) {
    rna <- x$counts[4:6, ]
    colnames(rna)[3] <- "OTHER-1"
    list(Peaks = x$counts[1:3, ], "Gene Expression" = rna)
  }, .package = "Seurat")
  expect_error(ReadMultiomeData(path, verbose = FALSE), "barcodes must match exactly")
})


test_that("a real H5 with one RNA feature reports the native reader limitation", {
  skip_if_not_installed("hdf5r")
  x <- make_multiome_input_fixture()
  counts <- x$counts[1:4, , drop = FALSE]
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::H5File$new(path, mode = "w")
  matrix <- h5$create_group("matrix")
  features <- matrix$create_group("features")
  matrix[["data"]] <- counts@x
  matrix[["indices"]] <- counts@i
  matrix[["indptr"]] <- counts@p
  matrix[["shape"]] <- dim(counts)
  matrix[["barcodes"]] <- colnames(counts)
  features[["id"]] <- rownames(counts)
  features[["name"]] <- rownames(counts)
  features[["feature_type"]] <- x$type[1:4]
  h5$close_all()
  expect_error(
    suppressMessages(ReadMultiomeData(path, genome = "hg38", verbose = FALSE)),
    "native Seurat H5 reader dropped matrix dimensions.*two features per modality"
  )
})
