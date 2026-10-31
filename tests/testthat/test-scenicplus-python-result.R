write_scenicplus_tables <- function(directory, suffix = "tsv") {
  tables <- list(
    tf_gene = c("TF\ttarget\timportance", "TF1\tG1\t0.5"),
    region_gene = c("region\tgene\tscore", "r1\tG1\t0.5"),
    tf_region = c("TF\tregion\tscore", "TF1\tr1\t0.5"),
    triplets = c("TF\tregion\tgene\tscore", "TF1\tr1\tG1\t0.5"),
    auc = c("\tTF1_direct_+/+_(1g)", "c1\t0.1", "c2\t0.2")
  )
  for (name in names(tables)) {
    path <- file.path(directory, paste0(name, ".", suffix))
    con <- if (grepl("gz$", suffix)) gzfile(path, "w") else file(path, "w")
    writeLines(if (grepl("csv", suffix, ignore.case = TRUE)) gsub("\t", ",", tables[[name]]) else tables[[name]], con)
    close(con)
  }
}

test_that("SCENIC+ discovers objects while respecting explicit inputs", {
  root <- withr::local_tempdir()
  object <- file.path(root, "scplus.h5mu")
  file.create(object)
  resolved <- scenicplus_resolve_python_inputs(python_result_dir = root, work_dir = "export")
  expect_identical(resolved$scplus_object, object)
  expect_identical(resolved$python_result_dir, "export")
  expect_null(scenicplus_resolve_python_inputs(python_result_dir = root, auto_prepare = FALSE)$scplus_object)
  expect_identical(scenicplus_resolve_python_inputs(python_result_dir = object)$scplus_object, object)
  unsupported <- file.path(root, "result.rds")
  file.create(unsupported)
  expect_error(scenicplus_resolve_object(unsupported), "supports")
  file.create(file.path(root, "run2.h5mu"))
  expect_error(scenicplus_resolve_object(root), "multiple")
})

test_that("SCENIC+ reads supported tables without preparing Python", {
  counts <- matrix(1:6, nrow = 3, dimnames = list(c("G1", "G2", "G3"), c("c2", "c1")))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts))
  testthat::local_mocked_bindings(.package = "scop", PrepareEnv = function(...) stop("Python prepared"))
  for (suffix in c("tsv", "CSV.gz")) {
    root <- withr::local_tempdir()
    write_scenicplus_tables(root, suffix)
    out <- suppressWarnings(RunSCENICPlus(srt, backend = "python", python_result_dir = root, verbose = FALSE))
    expect_equal(out@tools$SCENICPlus$auc[[1L]], c(0.2, 0.1))
    expect_identical(colnames(out[["scenicplus"]]), c("c2", "c1"))
  }
})

test_that("SCENIC+ imports a single complete directory without mixing tables", {
  root <- withr::local_tempdir()
  first <- file.path(root, "first")
  second <- file.path(root, "second")
  dir.create(first)
  dir.create(second)
  write_scenicplus_tables(first)
  file.rename(file.path(first, "tf_gene.tsv"), file.path(second, "tf_gene.tsv"))
  expect_false(scenicplus_has_all_tables(root))
  expect_error(scenicplus_read_python_result(root, srt = NULL, cells = "c1"), "single.*directory")
  file.rename(file.path(second, "tf_gene.tsv"), file.path(first, "tf_gene.tsv"))
  expect_identical(scenicplus_table_dir(root), first)
  writeLines(c("regulon\ttarget", "TF1_direct_+/+_(1g)\tWRONG"), file.path(root, "eregulons.tsv"))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(matrix(1:3, nrow = 3, dimnames = list(c("G1", "G2", "G3"), "c1"))))
  expect_false("WRONG" %in% unlist(scenicplus_read_python_result(root, srt, "c1")$regulon_list))
  write_scenicplus_tables(second)
  expect_null(scenicplus_table_dir(root))
})

test_that("SCENIC+ AUC preserves numeric labels and rejects unrelated labels", {
  path <- withr::local_tempfile(fileext = ".tsv")
  for (header in c("", "Cell", "CELL_ID", "Barcode", "Index", "UNNAMED: 0")) {
    for (labels in list(c("2", "1"), c("002", "001"))) {
      writeLines(c(paste0(header, "\tTF1"), paste0(labels, "\t", c("0.2", "0.1"))), path)
      auc <- scenicplus_normalize_auc_table(scenicplus_read_delim(path), rev(labels))
      expect_identical(rownames(auc), rev(labels))
      expect_identical(colnames(auc), "TF1")
      expect_equal(auc$TF1, c(0.1, 0.2))
    }
  }
  expect_error(scenicplus_normalize_auc_table(data.frame(TF1 = 1:2, row.names = c("wrongA", "wrongB")), c("c1", "c2")), "cell row names")
})
