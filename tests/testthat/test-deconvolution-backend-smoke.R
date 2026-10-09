# Optional genuine-backend checks. Package checks are read-only: these tests
# never install dependencies or download spacexr's likelihood matrices.
make_real_deconv_smoke_pair <- function(reverse_order = FALSE, full_depth = FALSE) {
  withr::local_seed(4815)
  genes <- paste0("g", seq_len(80))
  ref_ids <- paste0("ref", seq_len(40))
  labels <- stats::setNames(rep(c("A", "B"), each = 20), ref_ids)
  ref <- matrix(stats::rpois(80 * 40, 2), 80,
    dimnames = list(genes, ref_ids))
  ref[1:25, 1:20] <- ref[1:25, 1:20] + matrix(stats::rpois(25 * 20, 60), 25)
  ref[26:50, 21:40] <- ref[26:50, 21:40] + matrix(stats::rpois(25 * 20, 60), 25)
  st <- round(3 * cbind(
    rowMeans(ref[, 1:20]), rowMeans(ref[, 21:40]), rowMeans(ref),
    rowMeans(ref[, c(1:15, 21:25)]), rowMeans(ref[, c(1:5, 21:35)]),
    rep(0, 80), rowMeans(ref)
  ))
  dimnames(st) <- list(genes, paste0("spot", seq_len(7)))
  if (full_depth) {
    ref[61:80, ] <- ref[61:80, ] + 1000
    st[61:80, ] <- st[61:80, ] + 3000
    ref[1:60, "ref3"] <- 0
  }
  if (reverse_order) {
    # Construct in reverse order: Seurat's [, rev(cells)] preserves its old order.
    st <- st[, rev(colnames(st)), drop = FALSE]
    ref <- ref[, rev(colnames(ref)), drop = FALSE]
  }
  spatial <- SeuratObject::CreateSeuratObject(Matrix::Matrix(st, sparse = TRUE))
  reference <- SeuratObject::CreateSeuratObject(Matrix::Matrix(ref, sparse = TRUE))
  reference$celltype <- labels[colnames(reference)]
  spatial$x <- as.integer(sub("spot", "", colnames(spatial)))
  spatial$y <- spatial$x
  ref_subset <- c("ref39", "ref3", "ref28", "ref9", "ref25", "ref5",
    "ref32", "ref18", "ref37", "ref12", "ref23", "ref2")
  spot_subset <- c("spot5", "spot2", "spot6", "spot1", "spot4")
  spatial[["ALT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(st[, spot_subset, drop = FALSE], sparse = TRUE))
  reference[["ALT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(ref[, ref_subset, drop = FALSE], sparse = TRUE))
  list(spatial = spatial, reference = reference, features = genes[1:60])
}

deconv_smoke_check_installed <- function(packages, ...) {
  package_names <- sub("^.*/", "", packages)
  stats::setNames(lapply(package_names, requireNamespace, quietly = TRUE), package_names)
}

test_that("real SPOTlight matches explicitly aligned inputs for subset assays", {
  skip_on_cran()
  skip_if_not_installed("SPOTlight")
  skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
  testthat::local_mocked_bindings(
    check_r = deconv_smoke_check_installed,
    .package = "scop"
  )

  for (assay in c("RNA", "ALT")) {
    pair <- make_real_deconv_smoke_pair(reverse_order = identical(assay, "ALT"))
    spatial <- pair$spatial
    reference <- pair$reference
    if (identical(assay, "ALT")) {
      expect_identical(colnames(spatial), paste0("spot", 7:1))
      expect_identical(colnames(reference), paste0("ref", 40:1))
    }
    set.seed(813)
    wrapped <- RunSPOTlight(spatial, reference = reference,
      reference_label = "celltype", assay = assay, reference_assay = assay,
      features = pair$features, marker_top_n = 20, min_prop = 0,
      threads = 1, maxit = 100, verbose = FALSE)

    ref <- GetAssayData5(reference, assay = assay, layer = "counts")[pair$features, , drop = FALSE]
    st <- GetAssayData5(spatial, assay = assay, layer = "counts")[pair$features, , drop = FALSE]
    ref <- ref[, Matrix::colSums(ref) > 0, drop = FALSE]
    st <- st[, Matrix::colSums(st) > 0, drop = FALSE]
    keep_genes <- Matrix::rowSums(ref) > 0 & Matrix::rowSums(st) > 0
    ref <- ref[keep_genes, , drop = FALSE]
    st <- st[keep_genes, , drop = FALSE]
    labels <- as.character(reference[[]][colnames(ref), "celltype"])

    set.seed(813)
    direct <- SPOTlight::SPOTlight(x = ref, y = st, groups = labels,
      mgs = wrapped@tools$SPOTlight$marker_genes,
      gene_id = "gene", group_id = "cluster", weight_id = "weight",
      min_prop = 0, scale = TRUE, threads = 1, maxit = 100,
      verbose = FALSE)$mat
    direct <- direct / rowSums(direct)
    weights <- wrapped@tools$SPOTlight$weights
    expect_identical(rownames(weights), colnames(st))
    expect_equal(weights, direct[rownames(weights), colnames(weights), drop = FALSE],
      tolerance = 1e-8)
    expect_true(all(is.finite(weights)))
    expect_true(all(weights >= 0))
    expect_equal(unname(rowSums(weights)), rep(1, nrow(weights)), tolerance = 1e-8)
    omitted <- setdiff(colnames(spatial), rownames(weights))
    expect_true(all(is.na(wrapped@tools$SPOTlight$proportions[omitted, , drop = FALSE])))
    expect_gt(weights["spot1", "A"], weights["spot1", "B"])
    expect_gt(weights["spot2", "B"], weights["spot2", "A"])
    expect_identical(nrow(weights), if (identical(assay, "ALT")) 4L else 6L)
  }
})

test_that("real spacexr preprocessing preserves full-assay depths and identity", {
  skip_on_cran()
  skip_if_not_installed("spacexr")
  skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not("createRctd" %in% getNamespaceExports("spacexr"),
    "Requires the SummarizedExperiment spacexr API")
  pair <- make_real_deconv_smoke_pair(reverse_order = TRUE, full_depth = TRUE)
  spatial <- pair$spatial
  reference <- pair$reference
  spatial_depth <- Matrix::colSums(GetAssayData5(spatial, assay = "ALT", layer = "counts"))
  reference_depth <- Matrix::colSums(GetAssayData5(reference, assay = "ALT", layer = "counts"))
  expect_true(all(spatial_depth > 10000))
  expect_true(all(reference_depth > 5000))
  namespace_fun <- get_namespace_fun
  preprocessing_verified <- FALSE

  testthat::local_mocked_bindings(
    # Accept the installed Bioconductor API without requesting a GitHub-source
    # replacement. The constructor and profile calculations below are genuine.
    check_r = deconv_smoke_check_installed,
    get_namespace_fun = function(pkg, fun) {
      original <- namespace_fun(pkg, fun)
      if (!identical(pkg, "spacexr") || !identical(fun, "createRctd")) return(original)
      function(spatial_experiment, reference_experiment, ...) {
        st <- SummarizedExperiment::assay(spatial_experiment, "counts")
        ref <- SummarizedExperiment::assay(reference_experiment, "counts")
        st_numi <- stats::setNames(as.numeric(SummarizedExperiment::colData(spatial_experiment)$nUMI), colnames(st))
        ref_numi <- stats::setNames(as.numeric(SummarizedExperiment::colData(reference_experiment)$nUMI), colnames(ref))
        labels <- SummarizedExperiment::colData(reference_experiment)$cell_type
        expect_identical(st_numi, spatial_depth[colnames(st)])
        expect_identical(ref_numi, reference_depth[colnames(ref)])
        expect_identical(as.character(labels), as.character(reference[[]][colnames(ref), "celltype"]))
        expect_identical(rownames(SpatialExperiment::spatialCoords(spatial_experiment)), colnames(st))
        expect_false("ref3" %in% colnames(ref))
        expect_false("spot6" %in% colnames(st))
        expect_true(all(Matrix::colSums(st) < 10000))
        expect_true(all(Matrix::colSums(ref) < 5000))

        result <- original(spatial_experiment, reference_experiment, ...)
        result_numi <- stats::setNames(
          as.numeric(SummarizedExperiment::colData(result$spatial_experiment)$nUMI),
          colnames(result$spatial_experiment))
        expect_identical(result_numi, spatial_depth[colnames(result$spatial_experiment)])
        profiles <- result$cell_type_info$info[[1]]
        normalized <- sweep(as.matrix(ref), 2, ref_numi, "/")
        for (cell_type in levels(labels)) {
          expected <- rowMeans(normalized[, labels == cell_type, drop = FALSE])
          expect_equal(unname(profiles[, cell_type]), unname(expected), tolerance = 1e-12)
        }
        preprocessing_verified <<- TRUE
        # Deliberately stop before runRctd: this is a preprocessing check, not
        # a fitted-weight test, and must not download external Q matrices.
        stop("Completed genuine spacexr preprocessing smoke", call. = FALSE)
      }
    },
    .package = "scop"
  )
  expect_error(RunRCTD(spatial, reference = reference,
    reference_label = "celltype", assay = "ALT", reference_assay = "ALT",
    features = pair$features, min_cells = 2, max_cores = 1,
    rctd_mode = "full", verbose = FALSE,
    create_rctd_params = list(ref_n_cells_min = 2, ref_UMI_min = 5000, UMI_min = 10000)),
    "Completed genuine spacexr preprocessing smoke")
  expect_true(preprocessing_verified)
})
