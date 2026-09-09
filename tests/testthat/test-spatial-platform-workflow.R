make_platform_object <- function(n = 12L) {
  set.seed(421)
  counts <- matrix(rpois(40L * n, 3), 40L,
    dimnames = list(paste0("gene", 1:40), paste0("cell", seq_len(n))))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$x <- rep(seq_len(ceiling(n / 3)), length.out = n)
  object$y <- rep(1:3, length.out = n)
  object$label <- rep(c("A", "B"), length.out = n)
  object
}

test_that("spatial context is explicit and units are never guessed from metadata", {
  expect_false("SpatialDataInfo" %in% getNamespaceExports("scop"))
  expect_false("spatial_input_info" %in% getNamespaceExports("scop"))
  object <- make_platform_object()
  before <- object
  info <- spatial_input_info(object)
  expect_identical(info$data_type, "unknown")
  expect_identical(info$coordinate_units, "unknown")
  expect_equal(info$estimated_dense_bytes, 8 * 40 * 12)
  expect_identical(spatial_input_info(object, data_type = "cell", coordinate_units = "micron")$data_type, "cell")
  expect_identical(object, before)
  object[["fov1"]] <- SeuratObject::CreateFOV(object[[]][1:6, c("x", "y")], type = "centroids", assay = "RNA")
  object[["fov2"]] <- SeuratObject::CreateFOV(object[[]][7:12, c("x", "y")], type = "centroids", assay = "RNA", key = "second_")
  expect_error(spatial_input_info(object), "Multiple spatial images")
  expect_identical(spatial_input_info(object, image = "fov2")$cells, colnames(object)[7:12])
  object[["small"]] <- SeuratObject::CreateAssayObject(counts = GetAssayData5(object)[, 1:6])
  expect_error(spatial_input_info(object, assay = "small", image = "fov2"), "absent from assay")
})

test_that("Seurat plural segmentations are read as real cell polygons", {
  object <- make_platform_object()
  b <- data.frame(cell = rep(colnames(object), each = 4),
    x = rep(c(0, 2, 2, 0), ncol(object)), y = rep(c(0, 0, 2, 2), ncol(object)))
  fov <- SeuratObject::CreateFOV(object[[]][, c("x", "y")], type = "centroids", assay = "RNA")
  fov[["segmentations"]] <- SeuratObject::CreateSegmentation(b)
  object[["fov"]] <- fov
  expect_identical(spatial_input_info(object, image = "fov")$data_type, "cell")
  expect_s3_class(SpatialCellPlot(object, image = "fov"), "ggplot")
})

test_that("loaders validate provenance without installing or synthesizing data", {
  x <- make_platform_object()
  x[["fov"]] <- SeuratObject::CreateFOV(x[[]][, c("x", "y")], type = "centroids", assay = "RNA")
  testthat::local_mocked_bindings(.package = "Seurat", LoadXenium = function(data.dir, ...) x)
  out <- ReadSpatialData(tempdir(), "xenium", sample_id = "S1")
  expect_identical(out@misc$scop_spatial_input$fov$data_type, "cell")
  expect_identical(spatial_input_info(out, image = "fov")$coordinate_units, "micron")
  expect_error(spatial_input_info(out, image = "fov", coordinate_units = "pixel"), "conflicts")
  expect_identical(spatial_input_info(subset(out, cells = colnames(out)[1:4]), image = "fov")$data_type, "cell")
  expect_error(ReadSpatialData(tempdir(), "visium_hd", bin.size = c(8, 8)), "unique positive")
  testthat::local_mocked_bindings(.package = "Seurat", LoadXenium = function(...) stop("vendor read failed"))
  expect_error(ReadSpatialData(tempdir(), "xenium"), "vendor read failed")
  expect_null(x@misc$scop_spatial_input)
})

test_that("BANKSY and SmoothClust route to their own producers and reject stale labels", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  seen <- NULL
  producer <- function(srt, cluster_colname, tool_name = "BANKSY", ...) {
    seen <<- list(...)
    srt[[cluster_colname]] <- rep(c("1", "2"), length.out = ncol(srt))
    srt@tools[[tool_name]] <- list(status = "completed")
    srt
  }
  testthat::local_mocked_bindings(.package = "scop",
    RunStandardWorkflow = function(srt, ...) srt,
    RunBANKSY = producer,
    RunSmoothClust = function(srt, ..., tool_name = "SmoothClust") producer(srt, ..., tool_name = tool_name))
  for (method in c("BANKSY", "SmoothClust")) {
    out <- original(make_platform_object(), assay = "RNA", do_spot_qc = FALSE,
      do_spatial_variable_features = FALSE, do_deconvolution = FALSE,
      do_spatial_cluster = TRUE, spatial_cluster_method = method,
      spatial_q = if (method == "SmoothClust") 3L else NULL, verbose = FALSE)
    st <- out@tools$run_standard_spatial_workflow$stages
    expect_identical(st$actual_method[st$stage == "spatial_clustering"], paste0("Run", method))
    expect_identical(st$status[st$stage == "spatial_clustering"], "completed")
    if (method == "SmoothClust") expect_equal(seen$n_clusters, 3)
  }
  expect_error(original(make_platform_object(), assay = "RNA", do_spatial_cluster = TRUE,
    spatial_cluster_method = "BANKSY", spatial_q = 3, verbose = FALSE), "resolution")
  testthat::local_mocked_bindings(.package = "scop", RunBANKSY = function(srt, ...) srt)
  x <- make_platform_object(); x$BANKSY_cluster <- "old"; x@tools$BANKSY <- list(old = TRUE)
  expect_error(original(x, assay = "RNA", do_spot_qc = FALSE, do_spatial_variable_features = FALSE,
    do_spatial_cluster = TRUE, spatial_cluster_method = "BANKSY", verbose = FALSE))
  expect_true(all(x$BANKSY_cluster == "old"))
})

test_that("SpaNorm routes normalized data separately and QC partial is not success", {
  original <- getFromNamespace("run_standard_spatial_workflow", "scop")
  seen <- list()
  testthat::local_mocked_bindings(.package = "scop",
    RunSpaNorm = function(srt, assay, layer, ...) {
      seen$counts <<- GetAssayData5(srt, assay = assay, layer = layer)
      suppressWarnings(srt[["SpaNorm"]] <- SeuratObject::CreateAssayObject(data = log1p(seen$counts)))
      srt@tools$SpaNorm <- list(cells = colnames(srt)); srt
    },
    RunStandardWorkflow = function(srt, assay, do_normalization, HVF_method, ...) {
      seen$assay <<- assay; seen$normalize <<- do_normalization; seen$hvf <<- HVF_method; srt
    },
    RunSpatialVariableFeatures = function(srt, assay, ...) {
      seen$svf <<- assay; srt@tools$SpatialVariableFeatures <- list(result = data.frame(feature = rownames(srt))); srt
    },
    RunSpotSweeper = function(srt, ...) { srt@tools$SpotSweeper <- list(status = "partial"); srt })
  x <- make_platform_object()
  out <- original(x, assay = "RNA", do_spot_qc = FALSE, normalization_method = "SpaNorm", verbose = FALSE)
  expect_identical(seen$assay, "SpaNorm")
  expect_identical(seen$svf, "SpaNorm")
  expect_identical(seen$hvf, "mvp")
  expect_false(seen$normalize)
  expect_equal(GetAssayData5(out, assay = "RNA"), GetAssayData5(x))
  expect_equal(GetAssayData5(out, assay = "SpaNorm"), GetAssayData5(x))
  error <- tryCatch(original(x, assay = "RNA", do_spot_qc = FALSE,
    do_spatial_qc = TRUE, verbose = FALSE), error = identity)
  expect_s3_class(error, "error")
  st <- attr(error, "standard_spatial_stages")
  expect_identical(st$status[st$stage == "spatial_quality_control"], "failed")
  expect_error(original(x, assay = "RNA", spatial_data_type = "cell",
    do_spatial_cluster = TRUE, verbose = FALSE), "BayesSpace requires")
})

test_that("an explicit normalization refusal preserves externally normalized data", {
  x <- make_platform_object()
  x <- SeuratObject::SetAssayData(x, layer = "data", new.data = GetAssayData5(x) / 7)
  before <- GetAssayData5(x, layer = "data")
  testthat::local_mocked_bindings(.package = "scop",
    CheckDataType = function(...) "unknown",
    NormalizeData = function(...) stop("normalization must not run"))
  out <- CheckDataList(list(x), batch = "", assay = "RNA", do_normalization = FALSE,
    do_HVF_finding = FALSE, HVF = rownames(x)[1:5], nHVF = 5, verbose = FALSE)
  expect_equal(GetAssayData5(out$srt_list[[1]], layer = "data"), before)
})

test_that("SPARKX keeps sparse inputs and dense paths have a preflight bound", {
  x <- make_platform_object()
  testthat::local_mocked_bindings(.package = "scop",
    spatial_variable_run_sparkx = function(expr, ...) {
      expect_s4_class(expr, "dgCMatrix")
      data.frame(feature = rownames(expr), p_value = .2, q_value = .2, statistic = 1)
    })
  out <- RunSpatialVariableFeatures(x, method = "SPARKX", layer = "counts", verbose = FALSE)
  expect_true(length(out@tools$SpatialVariableFeatures$result$feature) > 0)
  expect_error(RunSpatialVariableFeatures(x, layer = "counts", backend = "r",
    max_dense_gb = 1e-12, verbose = FALSE), "max_dense_gb")
})

test_that("Assay5 scaling matches Seurat for sketch cells and repeated scaling", {
  object <- make_platform_object(30)
  selected <- colnames(object)[c(18:25, 1:8)]
  counts <- GetAssayData5(object)[, selected]
  object[["sketch"]] <- SeuratObject::CreateAssay5Object(counts = counts, data = log1p(counts))
  features <- rownames(counts)[1:12]
  original <- utils::getFromNamespace("ScaleData.Seurat", "Seurat")
  reference <- original(object, assay = "sketch", features = features, verbose = FALSE)
  out <- ScaleData(object, assay = "sketch", features = features, verbose = FALSE)
  expect_equal(GetAssayData5(out, assay = "sketch", layer = "scale.data"),
    GetAssayData5(reference, assay = "sketch", layer = "scale.data"), tolerance = 1e-8)
  expect_identical(colnames(GetAssayData5(out, assay = "sketch", layer = "scale.data")), colnames(object[["sketch"]]))
  out <- suppressWarnings(ScaleData(out, assay = "sketch", features = features[1:6], verbose = FALSE))
  expect_identical(rownames(GetAssayData5(out, assay = "sketch", layer = "scale.data")), features[1:6])
  expect_true(isTRUE(methods::validObject(out[["sketch"]])))
  object$batch <- rep(c("A", "B"), length.out = ncol(object))
  # Seurat's whole-object split.by path uses all object cells even for a sketch;
  # compare with its correctly subsetted reference context instead.
  reference <- original(subset(object, cells = colnames(object[["sketch"]])),
    assay = "sketch", features = features, split.by = "batch", verbose = FALSE)
  out <- ScaleData(object, assay = "sketch", features = features, split.by = "batch", verbose = FALSE)
  expect_equal(GetAssayData5(out, assay = "sketch", layer = "scale.data"),
    GetAssayData5(reference, assay = "sketch", layer = "scale.data"), tolerance = 1e-8)
})
