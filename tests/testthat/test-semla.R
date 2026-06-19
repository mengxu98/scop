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
  "semla" %in% rownames(utils::installed.packages())
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
  counts <- methods::as(Matrix::Matrix(1, nrow = 2, ncol = 2, sparse = TRUE), "dgCMatrix")
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
})
