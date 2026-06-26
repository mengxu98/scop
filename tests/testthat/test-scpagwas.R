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
      expect_identical(packages, "scPagwas")
      invisible(TRUE)
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
  runner <- function(single_data, gwas_data, block_annotation, output.dirs) {
    expect_s4_class(single_data, "Seurat")
    expect_identical(gwas_data, gwas)
    expect_identical(block_annotation, "hg38")
    expect_true(grepl("^/", output.dirs))
    single_data
  }

  with_mock_scpagwas(runner, {
    out <- RunscPagwas(
      srt = srt,
      gwas_data = gwas,
      celltype_meta = "celltype",
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
  runner <- function(single_data, block_annotation, unused = NULL) {
    expect_identical(single_data, normalizePath(single_data, mustWork = FALSE))
    expect_identical(block_annotation, normalizePath(block, mustWork = FALSE))
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
