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

test_that("installed GiottoClass is reused instead of reinstalling drieslab/Giotto", {
  seen <- list()
  testthat::local_mocked_bindings(
    check_r = function(packages, install = TRUE, ...) {
      seen <<- c(seen, list(list(packages = packages, install = install)))
      if (identical(packages, "GiottoClass") && !isTRUE(install)) {
        return(TRUE)
      }
      invisible(TRUE)
    },
    get_namespace_fun = function(...) {
      function(...) "converted"
    },
    .package = "scop"
  )
  expect_identical(srt_to_giotto(make_framework_seurat()), "converted")
  expect_identical(giotto_to_srt(structure(list(), class = "giotto")), "converted")
  expect_identical(
    vapply(seen, `[[`, character(1), "packages"),
    c("GiottoClass", "GiottoClass")
  )
  expect_false(any(vapply(seen, function(x) isTRUE(x$install), logical(1))))
})

test_that("missing GiottoClass requests drieslab/Giotto", {
  seen <- list()
  testthat::local_mocked_bindings(
    check_r = function(packages, install = TRUE, ...) {
      seen <<- c(seen, list(list(packages = packages, install = install)))
      if (identical(packages, "GiottoClass") && !isTRUE(install)) {
        return(FALSE)
      }
      invisible(TRUE)
    },
    get_namespace_fun = function(...) {
      function(...) "converted"
    },
    .package = "scop"
  )
  expect_identical(srt_to_giotto(make_framework_seurat()), "converted")
  expect_identical(
    vapply(seen, `[[`, character(1), "packages"),
    c("GiottoClass", "drieslab/Giotto")
  )
  expect_false(isTRUE(seen[[1]]$install))
  expect_true(isTRUE(seen[[2]]$install))
})

make_live_giotto_seurat <- function() {
  counts <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 3L),
    j = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L),
    x = c(4, 1, 2, 3, 5, 1, 2, 6),
    dims = c(3L, 3L),
    dimnames = list(paste0("Gene", 1:3), paste0("Spot", 1:3))
  )
  srt <- suppressMessages(Seurat::CreateSeuratObject(counts = counts, assay = "Spatial"))
  spots <- colnames(srt)
  # VisiumV1 matches the official GiottoClass converter's image/coordinate path.
  # A centroids FOV is not a reliable live fixture: seuratToGiottoV5 assigns
  # imagerow/imagecol names onto GetTissueCoordinates() output.
  srt[["slice1"]] <- methods::new(
    "VisiumV1",
    assay = "Spatial",
    key = "slice1_",
    image = array(0.5, dim = c(8L, 8L, 3L)),
    scale.factors = SeuratObject::scalefactors(
      spot = 1,
      fiducial = 1,
      hires = 1,
      lowres = 1
    ),
    coordinates = data.frame(
      tissue = 1L,
      row = c(0L, 0L, 1L),
      col = c(0L, 1L, 0L),
      imagerow = c(1, 1, 2),
      imagecol = c(1, 2, 1),
      row.names = spots
    ),
    spot.radius = 0.1
  )
  Seurat::NormalizeData(srt, assay = "Spatial", verbose = FALSE)
}

test_that("srt_to_giotto and giotto_to_srt round-trip with a real GiottoClass", {
  giotto_class <- thisutils::check_r("GiottoClass", install = FALSE, verbose = FALSE)
  if (!isTRUE(all(unlist(giotto_class, use.names = FALSE)))) {
    skip("GiottoClass is not installed")
  }

  # In-process tiny fixture: callr+load_all of the source tree exceeds 120s in
  # optional CI before conversion starts, and visium_human_pancreas_sub is far
  # larger than needed to exercise the official GiottoClass converters.
  old <- options(giotto.use_conda = FALSE, giotto.check_version = FALSE)
  on.exit(options(old), add = TRUE)

  srt <- make_live_giotto_seurat()
  g <- suppressWarnings(suppressMessages(srt_to_giotto(srt)))
  srt2 <- suppressWarnings(suppressMessages(giotto_to_srt(g)))

  expect_true(methods::is(g, "giotto"))
  expect_true(methods::is(srt2, "Seurat"))
  expect_true(ncol(srt2) > 0)
  layers <- tryCatch(
    SeuratObject::Layers(srt2, assay = "rna"),
    error = function(e) SeuratObject::Layers(srt2)
  )
  expect_true(all(c("counts", "data") %in% layers))
})
