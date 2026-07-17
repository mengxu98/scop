make_semla_spatial_seurat <- function(nspots = 80, nfeatures = 30) {
  data(visium_human_pancreas_sub, package = "scop")
  cells <- colnames(visium_human_pancreas_sub)[seq_len(nspots)]
  features <- rownames(visium_human_pancreas_sub)[seq_len(nfeatures)]
  srt <- suppressWarnings(visium_human_pancreas_sub[features, cells])
  srt <- Seurat::NormalizeData(srt, assay = SeuratObject::DefaultAssay(srt), verbose = FALSE)
  srt$semla_region <- rep(c("A", "B"), length.out = ncol(srt))
  srt
}

semla_installed_without_loading <- function() {
  requireNamespace("semla", quietly = TRUE)
}

skip_if_no_semla_backend <- function() {
  testthat::skip_on_os("mac")
  testthat::skip_if_not_installed("semla")
}

test_that("semla wrappers reject non-Seurat input before backend work", {
  expect_error(
    RunSemlaSpatialNetwork(matrix(1, nrow = 2, ncol = 2), verbose = FALSE),
    "Seurat"
  )
  expect_error(
    RunSemlaLocalG(matrix(1, nrow = 2, ncol = 2), features = "Gene1", verbose = FALSE),
    "Seurat"
  )
})

test_that("semla optional dependency error is clear", {
  testthat::skip_if(
    semla_installed_without_loading(),
    "semla is installed"
  )
  counts <- Matrix::sparseMatrix(
    i = c(1L, 2L),
    j = c(1L, 2L),
    x = c(1, 1),
    dims = c(2L, 2L)
  )
  srt <- Seurat::CreateSeuratObject(counts)
  expect_error(
    RunSemlaSpatialNetwork(srt, verbose = FALSE),
    "semla"
  )
})

test_that("semla_prepare_srt adds Staffli to spatial Seurat objects", {
  skip_if_no_semla_backend()
  srt <- make_semla_spatial_seurat()
  srt@tools[["Staffli"]] <- NULL
  out <- semla_prepare_srt(srt, verbose = FALSE)
  expect_s4_class(out, "Seurat")
  expect_false(is.null(out@tools[["Staffli"]]))
})

test_that("RunSemlaSpatialNetwork stores semla network results", {
  skip_if_no_semla_backend()
  srt <- make_semla_spatial_seurat()
  out <- RunSemlaSpatialNetwork(
    srt,
    nNeighbors = 3,
    tool_name = "SemlaSpatialNetwork",
    verbose = FALSE
  )
  expect_s4_class(out, "Seurat")
  expect_false(is.null(out@tools[["SemlaSpatialNetwork"]]))
  expect_type(out@tools[["SemlaSpatialNetwork"]][["network"]], "list")
  expect_equal(out@tools[["SemlaSpatialNetwork"]][["parameters"]][["nNeighbors"]], 3)
  expect_equal(out@tools[["SemlaSpatialNetwork"]][["schema_version"]], 1L)
  expect_equal(out@tools[["SemlaSpatialNetwork"]][["provenance"]][["backend_id"]], "semla")
})

test_that("RunSemlaLocalG writes local G metadata", {
  skip_if_no_semla_backend()
  srt <- make_semla_spatial_seurat()
  features <- rownames(srt)[1:2]
  out <- RunSemlaLocalG(
    srt,
    features = features,
    nNeighbors = 3,
    verbose = FALSE
  )
  expect_s4_class(out, "Seurat")
  expect_true(all(paste0("Gi[", features, "]") %in% colnames(out@meta.data)))
  expect_equal(out@tools[["SemlaLocalG"]][["schema_version"]], 1L)
  expect_equal(out@tools[["SemlaLocalG"]][["provenance"]][["producer"]], "RunSemlaLocalG")
})

test_that("Semla distance producers store schema-v1 provenance", {
  skip_if_no_semla_backend()
  srt <- make_semla_spatial_seurat()
  neighbors <- RunSemlaRegionNeighbors(
    srt, column_name = "semla_region", column_labels = "A",
    mode = "outer", column_key = "A_neighbor", verbose = FALSE
  )
  expect_equal(neighbors@tools[["SemlaRegionNeighbors"]][["schema_version"]], 1L)
  expect_equal(neighbors@tools[["SemlaRegionNeighbors"]][["provenance"]][["backend_id"]], "semla")

  radial <- RunSemlaRadialDistance(
    srt, column_name = "semla_region", selected_groups = "A",
    column_suffix = "A_distance", verbose = FALSE
  )
  expect_equal(radial@tools[["SemlaRadialDistance"]][["schema_version"]], 1L)
  expect_equal(radial@tools[["SemlaRadialDistance"]][["provenance"]][["producer"]], "RunSemlaRadialDistance")
})
