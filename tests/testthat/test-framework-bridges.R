make_framework_seurat <- function(n = 4) {
  counts <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 1L, 2L, 3L),
    j = c(1L, 1L, 1L, 2L, 2L, 2L),
    x = c(1, 2, 3, 4, 5, 6),
    dims = c(3L, n),
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:n))
  )
  srt <- suppressMessages(Seurat::CreateSeuratObject(counts = counts))
  srt$col <- seq_len(n)
  srt$row <- seq_len(n)
  srt
}

test_that("giotto bridges request the converter matching the installed Seurat version", {
  seurat_major <- as.integer(strsplit(
    as.character(utils::packageVersion("Seurat")), "\\."
  )[[1L]][1L])
  expected <- if (seurat_major >= 5L) "V5" else "V4"
  seen <- character(0)
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "GiottoClass")
      seen <<- c(seen, name)
      function(...) NULL
    },
    .package = "scop"
  )
  suppressMessages(srt_to_giotto(make_framework_seurat()))
  suppressMessages(giotto_to_srt(structure(list(), class = "giotto")))
  expect_identical(
    seen,
    c(paste0("seuratToGiotto", expected), paste0("giottoToSeurat", expected))
  )
})

test_that("srt_to_giotto calls the official forward converter with the selected image", {
  srt <- make_framework_seurat()
  srt[["slice1"]] <- suppressMessages(SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 1), row.names = c("Spot1", "Spot2")),
    type = "centroids", assay = SeuratObject::DefaultAssay(srt), key = "s1_"
  ))
  seen <- NULL
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "GiottoClass")
      function(...) {
        seen <<- list(...)
        "converted"
      }
    },
    .package = "scop"
  )
  out <- suppressWarnings(suppressMessages(srt_to_giotto(srt, image = "slice1")))
  expect_identical(out, "converted")
  expect_identical(seen$sobject@project.name, srt@project.name)
  expect_equal(ncol(seen$sobject), 2)
  expect_false(is.null(seen$spatial_assay))
  expect_false(isTRUE(seen$verbose))
})

test_that("srt_to_giotto lets callers override converter defaults", {
  seen <- NULL
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "GiottoClass")
      function(sobject, spatial_assay = "Spatial", verbose = TRUE) {
        seen <<- list(
          sobject = sobject,
          spatial_assay = spatial_assay,
          verbose = verbose
        )
        "converted"
      }
    },
    .package = "scop"
  )

  out <- suppressMessages(srt_to_giotto(
    make_framework_seurat(),
    spatial_assay = "RNA",
    verbose = TRUE
  ))

  expect_identical(out, "converted")
  expect_identical(seen$spatial_assay, "RNA")
  expect_true(seen$verbose)
})

test_that("srt_to_giotto rejects multi-image objects without an explicit image", {
  suppressWarnings(srt <- make_framework_seurat())
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- suppressMessages(SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 0), row.names = c("Spot1", "Spot2")),
    type = "centroids", assay = assay, key = "s1_"
  ))
  srt[["slice2"]] <- suppressMessages(SeuratObject::CreateFOV(
    data.frame(x = c(2, 3), y = c(1, 1), row.names = c("Spot3", "Spot4")),
    type = "centroids", assay = assay, key = "s2_"
  ))
  expect_error(srt_to_giotto(srt), "Multiple spatial images")
})

test_that("giotto_to_srt calls the official reverse converter", {
  gobject <- structure(list(mock = TRUE), class = "giotto")
  seen <- NULL
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) {
      expect_identical(package, "GiottoClass")
      expect_match(name, "^giottoToSeurat")
      function(...) {
        seen <<- list(...)
        "seurat_object"
      }
    },
    .package = "scop"
  )
  out <- giotto_to_srt(gobject)
  expect_identical(out, "seurat_object")
  expect_identical(seen$gobject, gobject)
})

test_that("giotto bridges run check_r before resolving the converter", {
  testthat::local_mocked_bindings(
    check_r = function(...) {
      stop("installation declined", call. = FALSE)
    },
    .package = "scop"
  )
  expect_error(srt_to_giotto(make_framework_seurat()), "installation declined")
  expect_error(giotto_to_srt(structure(list(), class = "giotto")), "installation declined")
})

test_that("srt_to_giotto and giotto_to_srt round-trip with a real GiottoClass", {
  skip_if_not_installed("GiottoClass")
  skip_if_not_installed("callr")

  package_path <- getNamespaceInfo(asNamespace("scop"), "path")
  source_tree <- file.exists(file.path(package_path, ".Rbuildignore"))
  result <- callr::r(
    function(package_path, source_tree, libpath) {
      if (!source_tree) {
        libpath <- unique(c(dirname(package_path), libpath))
      }
      .libPaths(libpath)
      if (source_tree) {
        pkgload::load_all(package_path, quiet = TRUE)
      }

      data_env <- new.env(parent = emptyenv())
      utils::data(
        list = "visium_human_pancreas_sub",
        package = "scop",
        envir = data_env
      )
      srt <- Seurat::NormalizeData(
        data_env$visium_human_pancreas_sub,
        assay = "Spatial",
        verbose = FALSE
      )
      to_giotto <- getExportedValue("scop", "srt_to_giotto")
      to_seurat <- getExportedValue("scop", "giotto_to_srt")
      g <- suppressWarnings(suppressMessages(to_giotto(srt)))
      srt2 <- suppressWarnings(suppressMessages(to_seurat(g)))

      list(
        is_giotto = methods::is(g, "giotto"),
        is_seurat = methods::is(srt2, "Seurat"),
        cells = ncol(srt2),
        layers = SeuratObject::Layers(srt2, assay = "rna")
      )
    },
    args = list(package_path, source_tree, .libPaths())
  )

  expect_true(result$is_giotto)
  expect_true(result$is_seurat)
  expect_true(result$cells > 0)
  expect_true(all(c("counts", "data") %in% result$layers))
})
