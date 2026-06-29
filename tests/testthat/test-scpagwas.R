pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)

make_scpagwas_seurat <- function() {
  counts <- matrix(
    c(
      4, 0, 1,
      0, 3, 2,
      1, 2, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:3), paste0("Cell", 1:3))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "A", "B")
  srt
}

make_scpagwas_gwas <- function() {
  data.frame(
    chrom = "1",
    pos = 123,
    rsid = "rs1",
    se = 0.1,
    beta = 0.2,
    maf = 0.3
  )
}

with_mock_scpagwas <- function(fun, code) {
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      expect_identical(packages, "sulab-wmu/scPagwas")
      invisible(TRUE)
    },
    scpagwas_namespace_available = function() {
      TRUE
    },
    get_namespace_fun = function(pkg, name) {
      expect_identical(pkg, "scPagwas")
      if (identical(name, "scPagwas_main")) {
        return(fun)
      }
      stop("not found")
    }
  )
  force(code)
}

test_that("RunscPagwas validates required GWAS columns early", {
  expect_error(
    RunscPagwas(single_data = "input.rds", gwas_data = data.frame(chrom = "1")),
    "missing required"
  )
})

test_that("RunscPagwas passes Seurat input and stores tools metadata", {
  srt <- make_scpagwas_seurat()
  gwas <- make_scpagwas_gwas()
  block <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chrom = "chr1", start = 1, end = 2, label = "Gene1"),
    block,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  runner <- function(Single_data, gwas_data, block_annotation, output.dirs) {
    expect_s4_class(Single_data, "Seurat")
    expect_identical(gwas_data, gwas)
    expect_s3_class(block_annotation, "data.frame")
    expect_true(grepl("^(/|[A-Za-z]:/)", output.dirs))
    Single_data
  }

  with_mock_scpagwas(runner, {
    out <- RunscPagwas(
      srt = srt,
      gwas_data = gwas,
      celltype_meta = "celltype",
      block_annotation = block,
      output.dirs = "relative-output",
      cleanup_soar = FALSE,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true("scPagwas" %in% names(out@tools))
  expect_equal(as.character(SeuratObject::Idents(out)), c("A", "A", "B"))
})

test_that("RunscPagwas supports RDS paths, custom block annotation, and list attr", {
  gwas <- make_scpagwas_gwas()
  single_data <- tempfile(fileext = ".rds")
  block <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chrom = "chr1", start = 1, end = 2, label = "Gene1"),
    block,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  runner <- function(Single_data, block_annotation, unused = NULL) {
    expect_identical(Single_data, normalizePath(single_data, mustWork = FALSE, winslash = "/"))
    expect_s3_class(block_annotation, "data.frame")
    list(score = 1)
  }

  with_mock_scpagwas(runner, {
    out <- RunscPaGWAS(
      single_data = single_data,
      gwas_data = gwas,
      block_annotation = block,
      return_seurat = FALSE,
      cleanup_soar = FALSE,
      unused = "kept only if formal exists",
      verbose = FALSE
    )
  })

  expect_equal(out$score, 1)
  expect_type(attr(out, "scPagwas"), "list")
})

test_that("RunscPagwas supplies upstream default package data explicitly", {
  srt <- make_scpagwas_seurat()
  gwas <- make_scpagwas_gwas()
  block <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chrom = "chr1", start = 1, end = 2, label = "Gene1"),
    block,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  runner <- function(Single_data, gwas_data, block_annotation, output.dirs, Pathway_list, chrom_ld) {
    expect_identical(Pathway_list, list(pathway = "genes"))
    expect_identical(chrom_ld, list(ld = "blocks"))
    list(ok = TRUE)
  }

  testthat::local_mocked_bindings(
    scpagwas_package_data_raw = function(name) {
      switch(name,
        Genes_by_pathway_kegg = list(pathway = "genes"),
        chrom_ld = list(ld = "blocks"),
        stop("unexpected data object")
      )
    }
  )
  with_mock_scpagwas(runner, {
    res <- RunscPagwas(
      srt = srt,
      gwas_data = gwas,
      block_annotation = block,
      output.dirs = "relative-output",
      return_seurat = FALSE,
      verbose = FALSE
    )
  })

  expect_true(isTRUE(res$ok))
})

test_that("RunscPagwas does not accept bare custom block selector", {
  expect_error(
    RunscPagwas(
      single_data = "input.rds",
      gwas_data = make_scpagwas_gwas(),
      block_annotation = "custom"
    ),
    "custom annotation path"
  )
})
