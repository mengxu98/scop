write_scenicplus_tables <- function(result_dir) {
  writeLines(
    c("TF\ttarget\timportance", "TF1\tG1\t0.5"),
    file.path(result_dir, "tf_gene.tsv")
  )
  writeLines(
    c("region\tgene\tscore", "r1\tG1\t0.5"),
    file.path(result_dir, "region_gene.tsv")
  )
  writeLines(
    c("TF\tregion\tscore", "TF1\tr1\t0.5"),
    file.path(result_dir, "tf_region.tsv")
  )
  writeLines(
    c("TF\tregion\tgene\tscore", "TF1\tr1\tG1\t0.5"),
    file.path(result_dir, "triplets.tsv")
  )
  writeLines(
    c("\tTF1_direct_+/+_(1g)", "c1\t0.1"),
    file.path(result_dir, "auc.tsv")
  )
}

test_that("SCENIC+ inputs are discovered from a result directory", {
  result_dir <- withr::local_tempdir()
  work_dir <- withr::local_tempdir()
  file.create(file.path(result_dir, "scplus_mdata.pkl"))

  resolved <- scenicplus_resolve_python_inputs(
    python_result_dir = result_dir,
    work_dir = work_dir
  )
  expect_identical(
    normalizePath(resolved$scplus_object),
    normalizePath(file.path(result_dir, "scplus_mdata.pkl"))
  )
  expect_identical(
    normalizePath(resolved$python_result_dir),
    normalizePath(work_dir)
  )
  expect_identical(
    normalizePath(resolved$work_dir),
    normalizePath(work_dir)
  )
})

test_that("existing SCENIC+ tables take precedence over a pickle object", {
  result_dir <- withr::local_tempdir()
  write_scenicplus_tables(result_dir)
  file.create(file.path(result_dir, "scplus_mdata.pkl"))

  resolved <- scenicplus_resolve_python_inputs(
    python_result_dir = result_dir
  )
  expect_null(resolved$scplus_object)
  expect_identical(
    normalizePath(resolved$python_result_dir),
    normalizePath(result_dir)
  )
})

test_that("a pickle path passed as python_result_dir becomes the object", {
  object <- tempfile(fileext = ".pkl")
  work_dir <- withr::local_tempdir()
  file.create(object)

  resolved <- scenicplus_resolve_python_inputs(
    python_result_dir = object,
    work_dir = work_dir
  )
  expect_identical(
    normalizePath(resolved$scplus_object),
    normalizePath(object)
  )
  expect_identical(
    normalizePath(resolved$python_result_dir),
    normalizePath(work_dir)
  )
})

test_that("SCENIC+ results are discovered in conventional directories", {
  root <- withr::local_tempdir()
  withr::local_dir(root)
  dir.create("scenic")
  file.create(file.path("scenic", "scplus_mdata.pkl"))

  resolved <- scenicplus_resolve_python_inputs()
  expect_identical(
    normalizePath(resolved$scplus_object),
    normalizePath(file.path("scenic", "scplus_mdata.pkl"))
  )
  expect_identical(
    normalizePath(resolved$python_result_dir, mustWork = FALSE),
    normalizePath("scenicplus_output/", mustWork = FALSE)
  )
})

test_that("automatic discovery can be disabled", {
  root <- withr::local_tempdir()
  withr::local_dir(root)
  dir.create("scenic")
  file.create(file.path("scenic", "scplus_mdata.pkl"))

  resolved <- scenicplus_resolve_python_inputs(auto_prepare = FALSE)
  expect_null(resolved$python_result_dir)
  expect_null(resolved$scplus_object)
})

test_that("discovered SCENIC+ objects prepare the Python environment", {
  result_dir <- withr::local_tempdir()
  file.create(file.path(result_dir, "scplus_mdata.pkl"))
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(modules, ...) stop("prepared")
  )

  expect_error(
    run_scenicplus_python(
      srt = NULL,
      envname = "scenicplus_env",
      conda = "auto",
      python_result_dir = result_dir,
      work_dir = result_dir
    ),
    "prepared"
  )
  expect_error(
    run_scenicplus_python(
      srt = NULL,
      envname = "scenicplus_env",
      conda = "auto",
      python_result_dir = result_dir,
      auto_prepare = FALSE
    ),
    "Missing official"
  )
})

test_that("h5mu objects are accepted as SCENIC+ input", {
  object <- tempfile(fileext = ".h5mu")
  file.create(object)

  expect_identical(
    normalizePath(scenicplus_resolve_object(object)),
    normalizePath(object)
  )
})

test_that("h5ad objects are rejected with export guidance", {
  object <- tempfile(fileext = ".h5ad")
  file.create(object)

  condition <- expect_error(scenicplus_resolve_object(object), "supports")
  expect_match(conditionMessage(condition), "python_result_dir", fixed = TRUE)
})

test_that("a directory with several SCENIC+ objects is rejected", {
  result_dir <- withr::local_tempdir()
  file.create(file.path(result_dir, "scplus_a.pkl"))
  file.create(file.path(result_dir, "scplus_b.pkl"))

  expect_error(scenicplus_resolve_object(result_dir), "multiple")
})

test_that("the Python backend explains its inputs when neither is supplied", {
  condition <- expect_error(
    run_scenicplus_python(
      srt = NULL,
      envname = "scenicplus_env",
      conda = "auto"
    ),
    "does not run"
  )
  message <- conditionMessage(condition)
  expect_match(message, "official SCENIC\\+ results", perl = TRUE)
  expect_match(message, "python_result_dir", fixed = TRUE)
  expect_match(message, "scplus_object", fixed = TRUE)
  expect_match(message, "allow_approximate", fixed = TRUE)
})

test_that("SCENIC+ result tables are discovered with supported suffixes", {
  result_dir <- withr::local_tempdir()
  writeLines(
    c("TF\ttarget\timportance", "TF1\tG1\t0.5"),
    file.path(result_dir, "tf_gene.tsv")
  )
  expect_identical(
    basename(scenicplus_table_path(result_dir, "tf_gene")),
    "tf_gene.tsv"
  )

  unlink(file.path(result_dir, "tf_gene.tsv"))
  con <- gzfile(file.path(result_dir, "tf_gene.csv.gz"), "w")
  writeLines(c("TF,target,importance", "TF1,G1,0.5"), con)
  close(con)

  path <- scenicplus_table_path(result_dir, "tf_gene")
  expect_identical(basename(path), "tf_gene.csv.gz")
  table <- scenicplus_read_required_table(
    result_dir,
    "tf_gene",
    c("TF", "target", "importance")
  )
  expect_identical(table$TF, "TF1")
  expect_identical(table$target, "G1")
  expect_equal(table$importance, 0.5)
})

test_that("SCENIC+ result tables are discovered case-insensitively", {
  result_dir <- withr::local_tempdir()
  writeLines(
    c("TF\ttarget\timportance", "TF1\tG1\t0.5"),
    file.path(result_dir, "TF_GENE.TSV")
  )
  expect_identical(
    basename(scenicplus_table_path(result_dir, "tf_gene")),
    "TF_GENE.TSV"
  )
})

test_that("a unique nested SCENIC+ table is used and ambiguous ones are rejected", {
  result_dir <- withr::local_tempdir()
  nested <- file.path(result_dir, "outs")
  dir.create(nested)
  writeLines(
    c("TF\ttarget\timportance", "TF1\tG1\t0.5"),
    file.path(nested, "tf_gene.tsv")
  )
  expect_identical(
    normalizePath(scenicplus_table_path(result_dir, "tf_gene")),
    normalizePath(file.path(nested, "tf_gene.tsv"))
  )

  second <- file.path(result_dir, "outs2")
  dir.create(second)
  writeLines(
    c("TF\ttarget\timportance", "TF2\tG2\t0.7"),
    file.path(second, "tf_gene.tsv")
  )
  expect_null(scenicplus_table_path(result_dir, "tf_gene"))
})

test_that("AUC tables with unnamed cell columns are normalized", {
  result_dir <- withr::local_tempdir()
  con <- gzfile(file.path(result_dir, "auc.tsv.gz"), "w")
  writeLines(
    c(
      "\tTF1_direct_+/+_(1g)",
      "c1\t0.1",
      "c2\t0.2",
      "c3\t0.3"
    ),
    con
  )
  close(con)

  auc <- scenicplus_read_required_table(result_dir, "auc", character())
  normalized <- scenicplus_normalize_auc_table(auc, cells = c("c1", "c2", "c3"))
  expect_identical(rownames(normalized), c("c1", "c2", "c3"))
  expect_identical(colnames(normalized), "TF1_direct_+/+_(1g)")
  expect_equal(unname(normalized[[1L]]), c(0.1, 0.2, 0.3))
})

test_that("AUC labels cannot be replaced by unrelated cells", {
  auc <- data.frame(TF = c(0.1, 0.9), row.names = c("wrong_A", "wrong_B"))
  expect_error(scenicplus_normalize_auc_table(auc, cells = c("cell_A", "cell_B")), "cell row names")
  rownames(auc) <- c("1", "2")
  expect_error(scenicplus_normalize_auc_table(auc, cells = c("cell_A", "cell_B")), "cell row names")
  expect_equal(scenicplus_normalize_auc_table(data.frame(TF = c(0.1, 0.9)), cells = c("cell_A", "cell_B")),
               data.frame(TF = c(0.1, 0.9), row.names = c("cell_A", "cell_B")))
})

test_that("numeric cell identifiers retain their labels and order", {
  path <- withr::local_tempfile(fileext = ".tsv")
  for (labels in list(c("2", "1"), c("002", "001"))) {
    writeLines(c("\tTF1", paste0(labels, "\t", c("0.2", "0.1"))), path)
    normalized <- scenicplus_normalize_auc_table(scenicplus_read_delim(path), cells = rev(labels))
    expect_identical(rownames(normalized), rev(labels))
    expect_identical(colnames(normalized), "TF1")
    expect_equal(normalized$TF1, c(0.1, 0.2))
  }
})

test_that("missing SCENIC+ tables report usable next steps", {
  result_dir <- withr::local_tempdir()
  writeLines("object", file.path(result_dir, "scplusmdata.pkl"))
  writeLines("metadata", file.path(result_dir, "eRegulon_direct.tsv"))

  condition <- expect_error(
    scenicplus_read_required_table(
      result_dir,
      "tf_gene",
      c("TF", "target", "importance")
    ),
    "tf_gene\\.tsv"
  )
  message <- conditionMessage(condition)
  expect_match(message, "scplusmdata\\.pkl", perl = TRUE)
  expect_match(message, "scplus_object", fixed = TRUE)
  expect_match(message, "eRegulon_direct\\.tsv", perl = TRUE)
  expect_match(message, "raw SCENIC\\+ output", perl = TRUE)
  expect_match(message, "TF", fixed = TRUE)
  expect_match(message, "importance", fixed = TRUE)
})

test_that("an empty SCENIC+ result directory is reported explicitly", {
  result_dir <- withr::local_tempdir()
  condition <- expect_error(
    scenicplus_read_required_table(result_dir, "auc", character()),
    "auc\\.tsv"
  )
  expect_match(conditionMessage(condition), "empty", fixed = TRUE)
})

test_that("ambiguous nested SCENIC+ tables suggest the containing directory", {
  result_dir <- withr::local_tempdir()
  for (sample in c("sample1", "sample2")) {
    dir.create(file.path(result_dir, sample))
    writeLines(
      c("TF\ttarget\timportance", "TF1\tG1\t0.5"),
      file.path(result_dir, sample, "tf_gene.tsv")
    )
  }
  condition <- expect_error(
    scenicplus_read_required_table(
      result_dir,
      "tf_gene",
      c("TF", "target", "importance")
    ),
    "Multiple"
  )
  message <- conditionMessage(condition)
  expect_match(message, "sample1", fixed = TRUE)
  expect_match(message, "sample2", fixed = TRUE)
})

test_that("a table with missing columns names the required schema", {
  result_dir <- withr::local_tempdir()
  writeLines(
    c("TF\tvalue", "TF1\t0.5"),
    file.path(result_dir, "tf_gene.tsv")
  )
  condition <- expect_error(
    scenicplus_read_required_table(
      result_dir,
      "tf_gene",
      c("TF", "target", "importance")
    ),
    "importance"
  )
  expect_match(conditionMessage(condition), "required columns", ignore.case = TRUE)
})

test_that("auto_prepare extracts tables from an h5mu object", {
  result_dir <- withr::local_tempdir()
  file.create(file.path(result_dir, "scplusmdata.h5mu"))
  counts <- matrix(
    1:9,
    nrow = 3,
    dimnames = list(c("G1", "G2", "G3"), c("c1", "c2", "c3"))
  )
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts))

  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(...) invisible(TRUE),
    resolve_conda = function(...) "conda",
    conda_python = function(...) "python",
    scenicplus_extract_python_tables = function(python_path, scplus_object, output_dir, verbose) {
      writeLines(
        c("TF\ttarget\timportance", "TF1\tG1\t0.9"),
        file.path(output_dir, "tf_gene.tsv")
      )
      writeLines(
        c("region\tgene\tscore", "r1\tG1\t0.8"),
        file.path(output_dir, "region_gene.tsv")
      )
      writeLines(
        c("TF\tregion\tscore", "TF1\tr1\t0.7"),
        file.path(output_dir, "tf_region.tsv")
      )
      writeLines(
        c("TF\tregion\tgene\tscore", "TF1\tr1\tG1\t0.6"),
        file.path(output_dir, "triplets.tsv")
      )
      writeLines(
        c(
          "\tTF1_direct_+/+_(1g)",
          "c1\t0.1",
          "c2\t0.2",
          "c3\t0.3"
        ),
        file.path(output_dir, "auc.tsv")
      )
      writeLines(
        c(
          "regulon\ttarget\ttarget_count",
          "TF1_direct_+/+_(1g)\tG1\t1"
        ),
        file.path(output_dir, "eregulons.tsv")
      )
    }
  )

  out <- suppressWarnings(
    run_scenicplus_python(
      srt = srt,
      envname = "scenicplus_env",
      conda = "auto",
      python_result_dir = result_dir,
      work_dir = result_dir
    )
  )
  expect_s4_class(out, "Seurat")
  expect_true("scenicplus" %in% names(out@assays))
  expect_identical(nrow(out[["scenicplus"]]), 1L)
  expect_identical(colnames(out[["scenicplus"]]), c("c1", "c2", "c3"))
  expect_identical(out@tools[["SCENICPlus"]][["tf_gene"]][["TF"]], "TF1")
})

test_that("RunSCENICPlus passes work_dir to the Python reader", {
  counts <- matrix(
    1:9,
    nrow = 3,
    dimnames = list(c("G1", "G2", "G3"), c("c1", "c2", "c3"))
  )
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts))
  seen <- NULL

  testthat::local_mocked_bindings(
    .package = "scop",
    run_scenicplus_python = function(..., work_dir) {
      seen <<- work_dir
      srt
    }
  )

  out <- RunSCENICPlus(
    srt,
    backend = "python",
    auto_prepare = FALSE,
    work_dir = "custom_output"
  )
  expect_identical(seen, "custom_output")
  expect_s4_class(out, "Seurat")
})

test_that("SCENIC+ tables from different result directories are never combined", {
  result_dir <- withr::local_tempdir()
  first <- file.path(result_dir, "sampleA")
  second <- file.path(result_dir, "sampleB")
  dir.create(first)
  dir.create(second)
  write_scenicplus_tables(second)
  file.rename(file.path(second, "tf_gene.tsv"), file.path(first, "tf_gene.tsv"))

  expect_false(scenicplus_has_all_tables(result_dir))
  expect_error(
    scenicplus_read_python_result(result_dir, srt = NULL, cells = "c1"),
    "single.*directory"
  )
})

test_that("SCENIC+ reading selects one complete directory including optional tables", {
  result_dir <- withr::local_tempdir()
  nested <- file.path(result_dir, "sampleA")
  dir.create(nested)
  write_scenicplus_tables(nested)
  writeLines(c("regulon\ttarget", "TF1_direct_+/+_(1g)\tWRONG"),
             file.path(result_dir, "eregulons.tsv"))
  counts <- matrix(1:3, nrow = 3, dimnames = list(c("G1", "G2", "G3"), "c1"))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts))

  expect_true(scenicplus_has_all_tables(result_dir))
  result <- scenicplus_read_python_result(result_dir, srt = srt, cells = "c1")
  expect_identical(result$tf_gene$TF, "TF1")
  expect_false("WRONG" %in% unlist(result$regulon_list))

  second <- file.path(result_dir, "sampleB")
  dir.create(second)
  write_scenicplus_tables(second)
  expect_false(scenicplus_has_all_tables(result_dir))
  expect_error(
    scenicplus_read_python_result(result_dir, srt = srt, cells = "c1"),
    "single.*directory"
  )
})

test_that("SCENIC+ optional tables cannot leak from a nested result", {
  result_dir <- withr::local_tempdir()
  write_scenicplus_tables(result_dir)
  dir.create(file.path(result_dir, "other"))
  writeLines(c("regulon\ttarget", "TF1_direct_+/+_(1g)\tWRONG"),
             file.path(result_dir, "other", "eregulons.tsv"))
  counts <- matrix(1:3, nrow = 3, dimnames = list(c("G1", "G2", "G3"), "c1"))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts))

  result <- scenicplus_read_python_result(result_dir, srt = srt, cells = "c1")
  expect_false("WRONG" %in% unlist(result$regulon_list))
})
