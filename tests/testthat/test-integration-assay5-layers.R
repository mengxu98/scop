test_that("legacy integration scaling handles split Assay5 data layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::rsparsematrix(80, 60, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 30)
  srt_list <- Seurat::SplitObject(srt, split.by = "batch")
  srt_list <- lapply(srt_list, function(x) {
    NormalizeData(x, verbose = FALSE)
  })
  merged <- Reduce(merge, srt_list)
  assay <- SeuratObject::DefaultAssay(merged)

  expect_false("data" %in% SeuratObject::Layers(merged[[assay]], search = NA))
  expect_gt(length(SeuratObject::Layers(merged[[assay]], search = "data")), 1)

  expect_no_error(
    suppressWarnings(
      RunIntegration(
        srt_merge = srt,
        batch = "batch",
        integration_methods = "Uncorrected",
        nHVF = 20,
        linear_reduction_dims = 5,
        linear_reduction_dims_use = 1:3,
        nonlinear_reduction = NULL,
        neighbor_k = 5,
        verbose = FALSE
      )
    )
  )
})

test_that("CheckDataMerge reuses log-normalized RNA merge without changing HVF output", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(2)
  counts <- Matrix::rsparsematrix(60, 40, density = 0.25)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 20)
  srt <- NormalizeData(srt, verbose = FALSE)

  checked <- CheckDataMerge(
    srt_merge = srt,
    batch = "batch",
    nHVF = 15,
    verbose = FALSE
  )

  expect_s4_class(checked$srt_merge, "Seurat")
  expect_equal(dim(checked$srt_merge), dim(srt))
  expect_equal(colnames(checked$srt_merge), colnames(srt))
  expect_length(checked$HVF, 15)
  expect_equal(SeuratObject::VariableFeatures(checked$srt_merge), checked$HVF)
})

test_that("CheckDataList preserves native SCT models and HVFs for split counts layers", {
  skip_if_not_installed("Seurat", "5.0.0")
  skip_if_not_installed("SeuratObject", "5.0.0")
  skip_if_not_installed("glmGamPoi")

  set.seed(31)
  counts <- Matrix::Matrix(matrix(rpois(120 * 100, lambda = 5), 120, 100), sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt <- Seurat::CreateSeuratObject(counts)
  srt[["RNA"]] <- split(srt[["RNA"]], f = rep(c("a", "b"), each = 50))
  layers <- SeuratObject::Layers(srt[["RNA"]], search = "counts")
  expected <- Seurat::SCTransform(srt, assay = "RNA", variable.features.n = 30, verbose = FALSE)
  expected <- Seurat::PrepSCTIntegration(
    object.list = list(expected),
    anchor.features = SeuratObject::VariableFeatures(expected[["SCT"]]),
    assay = "SCT",
    verbose = FALSE
  )[[1L]]

  checked <- CheckDataList(
    srt_list = list(srt),
    batch = "sample",
    assay = "RNA",
    normalization_method = "SCT",
    nHVF = 30,
    verbose = FALSE
  )
  actual <- checked$srt_list[[1L]]

  expect_identical(SeuratObject::Layers(actual[["RNA"]], search = "counts"), layers)
  expect_length(actual[["SCT"]]@SCTModel.list, 2L)
  expect_identical(SeuratObject::VariableFeatures(actual[["SCT"]]), SeuratObject::VariableFeatures(expected[["SCT"]]))
  expect_equal(Seurat::SCTResults(actual[["SCT"]], slot = "feature.attributes"), Seurat::SCTResults(expected[["SCT"]], slot = "feature.attributes"))
  expect_equal(Seurat::SCTResults(actual[["SCT"]], slot = "cell.attributes"), Seurat::SCTResults(expected[["SCT"]], slot = "cell.attributes"))
  expect_equal(SeuratObject::GetAssayData(actual, assay = "SCT", layer = "scale.data"), SeuratObject::GetAssayData(expected, assay = "SCT", layer = "scale.data"))
  expect_setequal(checked$HVF, SeuratObject::VariableFeatures(expected[["SCT"]]))
})

test_that("CheckDataList preserves native SCT models when merging the same batch", {
  skip_if_not_installed("Seurat", "5.0.0")
  skip_if_not_installed("SeuratObject", "5.0.0")
  skip_if_not_installed("glmGamPoi")

  set.seed(31)
  counts <- Matrix::Matrix(matrix(rpois(120 * 100, lambda = 5), 120, 100), sparse = TRUE)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  srt_list <- list(
    Seurat::CreateSeuratObject(counts[, 1:50]),
    Seurat::CreateSeuratObject(counts[, 51:100])
  )
  merged <- Reduce(merge, srt_list)
  expected <- Seurat::SCTransform(merged, assay = "RNA", variable.features.n = 30, verbose = FALSE)
  expected <- Seurat::PrepSCTIntegration(
    object.list = list(expected),
    anchor.features = SeuratObject::VariableFeatures(expected[["SCT"]]),
    assay = "SCT",
    verbose = FALSE
  )[[1L]]

  checked <- CheckDataList(
    srt_list = srt_list,
    batch = c("same", "same"),
    assay = "RNA",
    normalization_method = "SCT",
    nHVF = 30,
    verbose = FALSE
  )
  actual <- checked$srt_list[[1L]]

  expect_length(checked$srt_list, 1L)
  expect_identical(colnames(actual), colnames(merged))
  expect_identical(SeuratObject::Layers(actual[["RNA"]], search = "counts"), SeuratObject::Layers(merged[["RNA"]], search = "counts"))
  expect_length(actual[["SCT"]]@SCTModel.list, 2L)
  expect_identical(SeuratObject::VariableFeatures(actual[["SCT"]]), SeuratObject::VariableFeatures(expected[["SCT"]]))
  expect_equal(Seurat::SCTResults(actual[["SCT"]], slot = "feature.attributes"), Seurat::SCTResults(expected[["SCT"]], slot = "feature.attributes"))
  expect_equal(Seurat::SCTResults(actual[["SCT"]], slot = "cell.attributes"), Seurat::SCTResults(expected[["SCT"]], slot = "cell.attributes"))
  expect_equal(SeuratObject::GetAssayData(actual, assay = "SCT", layer = "scale.data"), SeuratObject::GetAssayData(expected, assay = "SCT", layer = "scale.data"))
  expect_setequal(checked$HVF, SeuratObject::VariableFeatures(expected[["SCT"]]))
})

test_that("NormalizeData handles split Assay5 counts layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::rsparsematrix(80, 60, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 30)
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)

  expect_setequal(
    SeuratObject::Layers(srt[[assay]], search = "counts"),
    c("counts.a", "counts.b")
  )
  expect_no_error(srt <- NormalizeData(srt, verbose = FALSE))
  expect_setequal(
    SeuratObject::Layers(srt[[assay]], search = "data"),
    c("data.a", "data.b")
  )
  expect_true(all(Matrix::colSums(SeuratObject::LayerData(srt, layer = "data.a")) > 0))
})

test_that("FindVariableFeatures handles split Assay5 counts layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::rsparsematrix(80, 60, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 30)
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)

  expect_no_error(
    srt <- FindVariableFeatures(srt, nfeatures = 20, verbose = FALSE)
  )
  expect_length(SeuratObject::VariableFeatures(srt), 20)
})

test_that("split-layer variable feature statistics match joined sparse statistics", {
  skip_if_not_installed("Matrix")

  set.seed(10)
  counts <- Matrix::rsparsematrix(100, 48, density = 0.18)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  layers <- list(
    counts[, 1:13, drop = FALSE],
    counts[, 14:31, drop = FALSE],
    counts[, 32:48, drop = FALSE]
  )
  joined <- do.call(cbind, layers)
  split_vst <- get("variable_features_vst_sparse_layers", asNamespace("scop"))
  joined_vst <- get("variable_features_vst_sparse", asNamespace("scop"))

  split_info <- split_vst(layers, nselect = 25, verbose = FALSE)
  joined_info <- joined_vst(joined, nselect = 25, verbose = FALSE)

  expect_equal(split_info$mean, joined_info$mean, tolerance = 1e-12)
  expect_equal(split_info$variance, joined_info$variance, tolerance = 1e-12)
  expect_equal(split_info$variance.expected, joined_info$variance.expected, tolerance = 1e-12)
  expect_equal(split_info$variance.standardized, joined_info$variance.standardized, tolerance = 1e-12)
  expect_equal(split_info$rank, joined_info$rank)
})

test_that("ScaleData handles split Assay5 data layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  set.seed(1)
  counts <- Matrix::rsparsematrix(80, 60, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 30)
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)
  srt <- NormalizeData(srt, verbose = FALSE)
  seurat_layered <- Seurat::FindVariableFeatures(
    srt,
    nfeatures = 20,
    verbose = FALSE
  )
  seurat_layered <- Seurat::ScaleData(
    seurat_layered,
    features = SeuratObject::VariableFeatures(seurat_layered),
    verbose = FALSE
  )
  srt <- FindVariableFeatures(srt, nfeatures = 20, verbose = FALSE)

  expect_no_error(
    srt <- ScaleData(
      srt,
      features = SeuratObject::VariableFeatures(srt),
      verbose = FALSE
    )
  )
  expect_true("scale.data" %in% SeuratObject::Layers(srt[[assay]], search = NA))
  expect_equal(
    nrow(SeuratObject::LayerData(srt, layer = "scale.data")),
    20L
  )
  expect_identical(
    SeuratObject::VariableFeatures(srt),
    SeuratObject::VariableFeatures(seurat_layered)
  )
  scale_data <- SeuratObject::LayerData(srt, layer = "scale.data")
  expect_equal(
    as.matrix(scale_data),
    as.matrix(SeuratObject::LayerData(seurat_layered, layer = "scale.data")),
    tolerance = 1e-12
  )
})

test_that("srt_to_adata stacks split Assay5 layers into sparse X", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("callr")

  set.seed(1)
  counts <- Matrix::rsparsematrix(20, 12, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)
  srt$batch <- rep(c("a", "b"), each = 6)
  assay <- SeuratObject::DefaultAssay(srt)
  srt[[assay]] <- split(srt[[assay]], f = srt$batch)

  layers <- SeuratObject::Layers(srt[[assay]], search = "counts")
  expect_setequal(layers, c("counts.a", "counts.b"))

  package_path <- getNamespaceInfo(asNamespace("scop"), "path")
  source_tree <- file.exists(file.path(package_path, ".Rbuildignore"))
  result <- callr::r(
    function(srt, assay, package_path, source_tree, libpath) {
      if (!source_tree) {
        libpath <- unique(c(dirname(package_path), libpath))
      }
      .libPaths(libpath)
      if (source_tree) {
        pkgload::load_all(package_path, quiet = TRUE)
      }
      if (!reticulate::py_module_available("anndata") ||
        !reticulate::py_module_available("scipy")) {
        return(list(available = FALSE))
      }

      options(scop_skip_python_prepare = TRUE)
      to_adata <- getExportedValue("scop", "srt_to_adata")
      adata <- to_adata(
        srt,
        assay_x = assay,
        layer_x = "counts",
        verbose = FALSE
      )
      list(
        available = TRUE,
        format = reticulate::py_to_r(adata$X$format),
        shape = unlist(reticulate::py_to_r(adata$X$shape)),
        obs_names = reticulate::py_to_r(adata$obs_names$to_list()),
        obs_rows = nrow(reticulate::py_to_r(adata$obs))
      )
    },
    args = list(srt, assay, package_path, source_tree, .libPaths())
  )
  skip_if_not(isTRUE(result$available), "Python modules anndata and scipy are unavailable")

  expect_identical(result$format, "csr")
  expect_equal(
    result$shape,
    c(ncol(srt), nrow(srt[[assay]]))
  )
  expect_identical(
    result$obs_names,
    colnames(srt)
  )
  expect_equal(result$obs_rows, ncol(srt))
})


test_that("SCT integration features are shared by every sparse layer model", {
  skip_if_not_installed("Seurat", "5.1.0")
  skip_if_not_installed("SeuratObject", "5.0.0")
  skip_if_not_installed("glmGamPoi")

  set.seed(491)
  counts <- rbind(matrix(rpois(100 * 120, lambda = 5), 100, 120), matrix(0, 40, 120))
  counts[101:120, 1:60] <- rnbinom(20 * 60, mu = 5, size = 0.2)
  counts[121:140, 61:120] <- rnbinom(20 * 60, mu = 5, size = 0.2)
  rownames(counts) <- c(paste0("shared", 1:100), paste0("only", 1:40))
  colnames(counts) <- paste0("cell", 1:120)
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object[["RNA"]] <- split(object[["RNA"]], f = rep(c("a", "b"), each = 60))
  native <- Seurat::SCTransform(object, variable.features.n = 30, verbose = FALSE)
  models <- native[["SCT"]]@SCTModel.list
  shared <- Reduce(intersect, lapply(models, function(model) rownames(model@feature.attributes)))
  layers <- SeuratObject::Layers(object[["RNA"]], search = "counts")
  expect_gt(length(setdiff(SeuratObject::VariableFeatures(native[["SCT"]]), shared)), 0)

  for (params in list(list(), list(HVF_source = "global"), list(HVF = c("only1", paste0("shared", 1:20))))) {
    checked <- suppressWarnings(do.call(CheckDataList, c(list(
      srt_list = list(native), batch = "sample", assay = "RNA", normalization_method = "SCT",
      do_normalization = FALSE, do_HVF_finding = FALSE, nHVF = 30, verbose = FALSE
    ), params)))
    actual <- checked$srt_list[[1L]]
    expected <- native
    suppressWarnings(expected[["RNA"]] <- SeuratObject::CreateAssayObject(counts = GetAssayData5(expected, assay = "RNA", layer = "counts")))
    expected <- Seurat::PrepSCTIntegration(list(expected), anchor.features = checked$HVF, assay = "SCT", verbose = FALSE)[[1L]]
    expect_true(all(checked$HVF %in% shared))
    expect_true(all(SeuratObject::VariableFeatures(actual[["SCT"]]) %in% shared))
    expect_identical(SeuratObject::Layers(actual[["RNA"]], search = "counts"), layers)
    expect_equal(actual[["SCT"]]@SCTModel.list, models)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "SCT", layer = "scale.data"), SeuratObject::GetAssayData(expected, assay = "SCT", layer = "scale.data"))
  }
  expect_error(suppressWarnings(CheckDataList(
    srt_list = list(native), batch = "sample", assay = "RNA", normalization_method = "SCT",
    do_normalization = FALSE, do_HVF_finding = FALSE, HVF = "only1", verbose = FALSE
  )), "No HVF available in every SCT model")
  object$batch <- rep(c("a", "b"), each = 60)
  integrated <- suppressWarnings(RPCA_integrate(
    srt_merge = object, batch = "batch", assay = "RNA", append = FALSE,
    normalization_method = "SCT", HVF = c("only1", paste0("shared", 1:20)), nHVF = 30,
    linear_reduction_dims = 5, linear_reduction_dims_use = 1:5, nonlinear_reduction = NULL,
    neighbor_k = 10L, IntegrateLayers_params = list(dims = 1:5, k.anchor = 5, k.weight = 10), verbose = FALSE
  ))
  expect_true(all(SeuratObject::VariableFeatures(integrated[["SCT"]]) %in% shared))
  expect_false("only1" %in% integrated@misc$RPCA_HVF)
  expect_identical(SeuratObject::Layers(integrated[["RNA"]], search = "counts"), layers)
  expect_equal(Seurat::SCTResults(integrated[["SCT"]], slot = "feature.attributes"), Seurat::SCTResults(native[["SCT"]], slot = "feature.attributes"))
})

test_that("legacy SCT integration preserves native layer models and reference residuals", {
  skip_if_not_installed("Seurat", "5.1.0")
  skip_if_not_installed("SeuratObject", "5.0.0")
  skip_if_not_installed("glmGamPoi")

  set.seed(490)
  counts <- matrix(rpois(120 * 230, lambda = 5), 120, 230,
    dimnames = list(paste0("g", 1:120), paste0("c", 1:230))
  )
  objects <- lapply(list(1:120, 121:230), function(cells) {
    object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts[, cells], sparse = TRUE))
    object[["RNA"]] <- split(object[["RNA"]], f = c(rep("a", length(cells) / 2 - 5), rep("b", length(cells) / 2 + 5)))
    Seurat::SCTransform(object, variable.features.n = 60, verbose = FALSE)
  })
  names(objects) <- c("one", "two")
  features <- paste0("g", 1:60)

  for (config in list(
    list(reference = 1L, tree = NULL, preserve = FALSE, features = features),
    list(reference = NULL, tree = matrix(c(-2L, -1L), nrow = 1), preserve = TRUE, features = features),
    list(reference = NULL, tree = matrix(c(-2L, -1L), nrow = 1), preserve = FALSE, features = paste0("g", 1:80)),
    list(reference = NULL, tree = NULL, preserve = FALSE, features = NULL, weight_features = paste0("g", 1:80))
  )) {
    set.seed(11)
    prepared <- lapply(objects, function(object) {
      suppressWarnings(object[["RNA"]] <- SeuratObject::CreateAssayObject(counts = GetAssayData5(object, assay = "RNA", layer = "counts")))
      object
    })
    prepared <- Seurat::PrepSCTIntegration(prepared, anchor.features = features, assay = "SCT", verbose = FALSE)
    anchors <- suppressWarnings(Seurat::FindIntegrationAnchors(prepared,
      normalization.method = "SCT", anchor.features = features, dims = 1:5,
      k.anchor = 5, k.filter = NA, reference = config$reference, verbose = FALSE
    ))
    residual_features <- if (is.null(config$features)) config$weight_features else config$features
    for (i in seq_along(anchors@object.list)) {
      object <- Seurat::GetResidual(anchors@object.list[[i]], features = residual_features, assay = "SCT", verbose = FALSE)
      models <- object[["SCT"]]@SCTModel.list
      sizes <- vapply(models, function(model) nrow(model@cell.attributes), integer(1))
      object[["SCT"]]@SCTModel.list <- models[which.max(sizes)]
      anchors@object.list[[i]] <- object
    }
    expected <- suppressWarnings(Seurat::IntegrateData(anchors,
      normalization.method = "SCT", new.assay.name = "Seuratcorrected", features = config$weight_features, features.to.integrate = config$features,
      dims = 1:5, k.weight = 10, sample.tree = config$tree, preserve.order = config$preserve, verbose = FALSE
    ))
    actual <- suppressWarnings(Seurat_integrate(
      srt_list = objects, batch = names(objects), assay = "RNA", normalization_method = "SCT",
      do_normalization = FALSE, do_HVF_finding = FALSE, HVF = features, nHVF = 60,
      do_scaling = FALSE, linear_reduction_dims = 5, linear_reduction_dims_use = 1:5,
      nonlinear_reduction = NULL, neighbor_k = 10, verbose = FALSE,
      FindIntegrationAnchors_params = list(dims = 1:5, k.anchor = 5, k.filter = NA, reference = config$reference),
      IntegrateData_params = list(dims = 1:5, k.weight = 10, sample.tree = config$tree, preserve.order = config$preserve, features = config$weight_features, features.to.integrate = config$features)
    ))
    expected_sct <- merge(objects[[1L]][["SCT"]], y = objects[[2L]][["SCT"]])
    expect_identical(colnames(actual), colnames(expected))
    expect_length(actual[["SCT"]]@SCTModel.list, 4L)
    expect_equal(actual[["SCT"]]@SCTModel.list, expected_sct@SCTModel.list)
    expect_equal(actual[["Seuratcorrected"]]@SCTModel.list, expected[["Seuratcorrected"]]@SCTModel.list)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "Seuratcorrected", layer = "scale.data"), SeuratObject::GetAssayData(expected, assay = "Seuratcorrected", layer = "scale.data"), tolerance = 1e-10)
  }
})


test_that("SCT residual completion uses the recorded source assay", {
  skip_if_not_installed("Seurat", "5.1.0")
  skip_if_not_installed("SeuratObject", "5.0.0")
  skip_if_not_installed("glmGamPoi")

  set.seed(495)
  counts <- Matrix::Matrix(matrix(rpois(100 * 100, lambda = 5), 100, 100), sparse = TRUE)
  rownames(counts) <- paste0("g", 1:100)
  colnames(counts) <- paste0("c", 1:100)
  object <- Seurat::CreateSeuratObject(counts)
  object <- suppressWarnings(SeuratObject::RenameAssays(object, RNA = "GeneExpression"))
  object[["GeneExpression"]] <- split(object[["GeneExpression"]], f = rep(c("a", "b"), each = 50))
  native <- Seurat::SCTransform(object, assay = "GeneExpression", variable.features.n = 30, verbose = FALSE)
  for (i in seq_along(native[["SCT"]]@SCTModel.list)) native[["SCT"]]@SCTModel.list[[i]]@umi.assay <- "GeneExpression"
  features <- head(setdiff(rownames(native[["SCT"]]), rownames(SeuratObject::GetAssayData(native, assay = "SCT", layer = "scale.data"))), 10)
  expect_length(features, 10L)
  for (decoy in c(FALSE, TRUE)) {
    input <- native
    if (decoy) input[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = counts * 10)
    expected <- input
    suppressWarnings(expected[["GeneExpression"]] <- SeuratObject::CreateAssayObject(counts = GetAssayData5(expected, assay = "GeneExpression", layer = "counts")))
    expected <- Seurat::GetResidual(expected, assay = "SCT", umi.assay = "GeneExpression", features = features, verbose = FALSE)
    expected <- Seurat::PrepSCTIntegration(list(expected), anchor.features = features, assay = "SCT", verbose = FALSE)[[1L]]
    checked <- suppressWarnings(CheckDataList(list(input),
      batch = "sample", assay = "GeneExpression",
      normalization_method = "SCT", do_normalization = FALSE, do_HVF_finding = FALSE,
      HVF = features, nHVF = 30, verbose = FALSE
    ))
    actual <- checked$srt_list[[1L]]
    expect_identical(SeuratObject::Layers(actual[["GeneExpression"]], search = "counts"), SeuratObject::Layers(input[["GeneExpression"]], search = "counts"))
    expect_equal(actual[["SCT"]]@SCTModel.list, native[["SCT"]]@SCTModel.list)
    expect_equal(SeuratObject::GetAssayData(actual, assay = "SCT", layer = "scale.data"), SeuratObject::GetAssayData(expected, assay = "SCT", layer = "scale.data"))
  }
  native[["SCT"]]@SCTModel.list[[2L]]@umi.assay <- "RNA"
  expect_error(suppressWarnings(CheckDataList(list(native),
    batch = "sample", assay = "GeneExpression",
    normalization_method = "SCT", do_normalization = FALSE, do_HVF_finding = FALSE,
    HVF = features, verbose = FALSE
  )), "must share one UMI assay")
})
