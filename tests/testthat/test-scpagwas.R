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

with_mock_scpagwas <- function(fun, code, backend_funs = list()) {
  testthat::local_mocked_bindings(
    .package = "scop",
    check_r = function(packages, ...) {
      expect_identical(packages, "sulab-wmu/scPagwas")
      invisible(TRUE)
    },
    get_namespace_fun = function(pkg, name) {
      expect_identical(pkg, "scPagwas")
      if (identical(name, "scPagwas_main2")) {
        return(fun)
      }
      if (name %in% names(backend_funs)) {
        return(backend_funs[[name]])
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
  runner <- function(
    Single_data,
    gwas_data,
    block_annotation,
    output.dirs,
    seurat_return,
    singlecell,
    celltype
  ) {
    expect_s4_class(Single_data, "Seurat")
    expect_identical(gwas_data, gwas)
    expect_s3_class(block_annotation, "data.frame")
    expect_true(grepl("^(/|[A-Za-z]:/)", output.dirs))
    expect_identical(seurat_return, TRUE)
    expect_identical(singlecell, FALSE)
    expect_identical(celltype, TRUE)
    Single_data
  }

  with_mock_scpagwas(runner, {
    out <- RunscPagwas(
      srt = srt,
      gwas_data = gwas,
      celltype_meta = "celltype",
      singlecell = FALSE,
      celltype = TRUE,
      block_annotation = block,
      output.dirs = "relative-output",
      cleanup_soar = FALSE,
      verbose = FALSE
    )
  })

  expect_s4_class(out, "Seurat")
  expect_true("scPagwas" %in% names(out@tools))
  expect_null(out@tools$scPagwas$result)
  expect_equal(as.character(SeuratObject::Idents(out)), c("A", "A", "B"))
})

test_that("RunscPagwas supports RDS paths, custom block annotation, and list attr", {
  gwas <- make_scpagwas_gwas()
  single_data <- tempfile(fileext = ".rds")
  saveRDS(make_scpagwas_seurat(), single_data)
  block <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chrom = "chr1", start = 1, end = 2, label = "Gene1"),
    block,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  runner <- function(Single_data, block_annotation, seurat_return, unused = NULL) {
    expect_s4_class(Single_data, "Seurat")
    expect_s3_class(block_annotation, "data.frame")
    expect_identical(seurat_return, FALSE)
    list(score = 1)
  }

  with_mock_scpagwas(runner, {
    out <- RunscPagwas(
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

test_that("RunscPagwas accepts and validates GWAS file paths", {
  gwas_file <- tempfile(fileext = ".tsv")
  write.table(
    make_scpagwas_gwas(),
    gwas_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  single_data <- tempfile(fileext = ".rds")
  saveRDS(make_scpagwas_seurat(), single_data)
  block <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chrom = "chr1", start = 1, end = 2, label = "Gene1"),
    block,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  runner <- function(gwas_data, ...) {
    expect_identical(
      gwas_data,
      normalizePath(gwas_file, mustWork = TRUE, winslash = "/")
    )
    list(ok = TRUE)
  }

  with_mock_scpagwas(runner, {
    out <- RunscPagwas(
      single_data = single_data,
      gwas_data = gwas_file,
      block_annotation = block,
      return_seurat = FALSE,
      verbose = FALSE
    )
  })

  expect_identical(out$ok, TRUE)
})

test_that("RunscPagwas rejects GWAS files with missing columns early", {
  gwas_file <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(chrom = "1", pos = 123),
    gwas_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  expect_error(
    RunscPagwas(single_data = "input.rds", gwas_data = gwas_file),
    "missing required"
  )
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
    .package = "scop",
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

test_that("scPagwas compatibility runner updates current Seurat calls", {
  runner <- function() {
    result <- Seurat::AddModuleScore(
      object,
      assay = assay,
      scPagwas_topgenes,
      name = "scPagwas.TRS.Score"
    )
    Seurat::FindAllMarkers(result, slot = "data")
  }
  single_input <- function() {
    GetAssayData(object, slot = "data", assay = assay)
  }
  correct_bg <- function() {
    GetAssayData(object, slot = "data", assay = assay)
  }

  testthat::local_mocked_bindings(
    .package = "scop",
    get_namespace_fun = function(pkg, name) {
      expect_identical(pkg, "scPagwas")
      switch(name,
        Single_data_input = single_input,
        Get_CorrectBg_p = correct_bg,
        stop("not found")
      )
    }
  )
  patched <- scpagwas_prepare_runner(runner)
  runner_code <- paste(deparse(body(patched)), collapse = " ")
  single_code <- paste(
    deparse(body(get("Single_data_input", envir = environment(patched)))),
    collapse = " "
  )

  expect_match(runner_code, "features = list\\(scPagwas_topgenes\\)")
  expect_match(runner_code, "layer = \\\"data\\\"")
  expect_match(single_code, "layer = \\\"data\\\"")
})

test_that("RunscPaGWAS is a deprecated compatibility alias", {
  testthat::local_mocked_bindings(
    .package = "scop",
    RunscPagwas = function(...) "ok"
  )

  expect_warning(out <- RunscPaGWAS(), "deprecated")
  expect_identical(out, "ok")
})

test_that("PlotScPagwas returns score and significance plots", {
  srt <- make_scpagwas_seurat()
  embeddings <- matrix(
    c(0, 0, 1, 0, 0, 1),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(colnames(srt), c("UMAP_1", "UMAP_2"))
  )
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "UMAP_",
    assay = "RNA"
  )
  srt$scPagwas.gPAS.score <- c(0.1, 0.5, 0.9)
  srt$scPagwas.TRS.Score1 <- c(-0.2, 0.3, 0.7)
  srt$Random_Correct_BG_adjp <- c(0.01, 0.2, 0.03)
  output_dir <- tempfile()

  plots <- PlotScPagwas(
    srt,
    reduction = "umap",
    output.dir = output_dir,
    do_plot = FALSE
  )

  expect_named(
    plots,
    c("scPagwas.gPAS.score", "scPagwas.TRS.Score1", "significant_cells")
  )
  expect_s3_class(plots[[1]], "ggplot")
  expect_length(list.files(output_dir, pattern = "\\.pdf$"), 3)
})
