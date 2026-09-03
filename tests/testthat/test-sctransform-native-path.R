test_that("SCTransform gates the validated native and fallback branches", {
  gate <- get("sct_fast_path_supported", asNamespace("scop"))
  baseline <- list(
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    residual.features = NULL,
    conserve.memory = FALSE,
    vst.flavor = "v2",
    do.scale = FALSE,
    do.center = TRUE,
    return.only.var.genes = TRUE,
    extra_args = list(),
    ncells = 5000,
    variable.features.n = 3000,
    variable.features.rv.th = 1.3,
    clip.range = c(-2, 2),
    regression_ok = TRUE
  )
  expect_true(do.call(gate, baseline))

  supported <- list(
    do.correct.umi = FALSE,
    do.scale = TRUE,
    do.center = FALSE,
    return.only.var.genes = FALSE
  )
  for (name in names(supported)) {
    candidate <- baseline
    candidate[[name]] <- supported[[name]]
    expect_true(do.call(gate, candidate), info = name)
  }
  candidate <- baseline
  candidate["variable.features.n"] <- list(NULL)
  expect_true(do.call(gate, candidate), info = "variable.features.n=NULL")

  unsupported <- list(
    reference.SCT.model = structure(list(), class = "SCTModel"),
    residual.features = "G1",
    conserve.memory = TRUE,
    vst.flavor = "v1",
    extra_args = list(min_cells = 1),
    ncells = 0,
    variable.features.n = NA_real_,
    variable.features.rv.th = -1,
    clip.range = c(1, -1),
    regression_ok = FALSE
  )
  for (name in names(unsupported)) {
    candidate <- baseline
    candidate[[name]] <- unsupported[[name]]
    expect_false(do.call(gate, candidate), info = name)
  }
})

test_that("SCTransform methods track Seurat parameters", {
  for (method in c("SCTransform.default", "SCTransform.Seurat")) {
    reference <- names(formals(get(method, asNamespace("Seurat"))))
    candidate <- names(formals(get(method, asNamespace("scop"))))
    expect_identical(sort(candidate), sort(c(reference, "cores")))
  }
})

test_that("SCTransform validates native regression designs", {
  validate <- get("sct_regression_supported", asNamespace("scop"))
  cells <- paste0("C", seq_len(12L))
  cell_attr <- data.frame(
    numeric_covariate = seq_along(cells),
    batch = factor(rep(c("a", "b"), length.out = length(cells))),
    row.names = cells
  )
  latent <- data.frame(
    latent_covariate = stats::runif(length(cells)),
    row.names = rev(cells)
  )
  expect_true(validate(NULL, NULL, cell_attr, cells))
  expect_true(validate("numeric_covariate", NULL, cell_attr, cells))
  expect_true(validate(c("numeric_covariate", "batch"), latent, cell_attr, cells))
  expect_false(validate("missing", NULL, cell_attr, cells))
  cell_attr$numeric_covariate[[1L]] <- NA_real_
  expect_false(validate("numeric_covariate", NULL, cell_attr, cells))
  expect_false(validate(NULL, latent[-1L, , drop = FALSE], cell_attr, cells))
})

test_that("SCTransform native path preserves Seurat output semantics", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glmGamPoi")

  set.seed(3)
  genes <- paste0("G", seq_len(250L))
  cells <- paste0("C", seq_len(150L))
  umi <- Matrix::rsparsematrix(
    nrow = length(genes),
    ncol = length(cells),
    density = 0.01,
    rand.x = function(n) stats::rpois(n, lambda = 4) + 1
  )
  dimnames(umi) <- list(genes, cells)
  umi[1L, ] <- 1
  reference <- Seurat::CreateSeuratObject(counts = umi)
  candidate <- Seurat::CreateSeuratObject(counts = umi)
  seurat_sct <- get("SCTransform.Seurat", asNamespace("Seurat"))
  reference <- seurat_sct(
    object = reference,
    vst.flavor = "v2",
    variable.features.n = 100,
    ncells = 150,
    seed.use = 103,
    verbose = FALSE
  )
  candidate <- SCTransform(
    candidate,
    vst.flavor = "v2",
    variable.features.n = 100,
    ncells = 150,
    seed.use = 103,
    verbose = FALSE
  )

  expect_identical(
    SeuratObject::VariableFeatures(candidate),
    SeuratObject::VariableFeatures(reference)
  )
  reference_counts <- SeuratObject::LayerData(reference[["SCT"]], "counts")
  candidate_counts <- SeuratObject::LayerData(candidate[["SCT"]], "counts")
  expect_identical(candidate_counts, reference_counts)
  reference_scale <- SeuratObject::LayerData(reference[["SCT"]], "scale.data")
  candidate_scale <- SeuratObject::LayerData(candidate[["SCT"]], "scale.data")
  expect_lt(max(abs(candidate_scale - reference_scale)), 1e-10)

  branch_args <- list(
    list(do.correct.umi = FALSE),
    list(do.scale = TRUE),
    list(do.center = FALSE),
    list(return.only.var.genes = FALSE),
    list(variable.features.n = NULL, variable.features.rv.th = 0.5)
  )
  for (branch in branch_args) {
    common <- list(
      object = Seurat::CreateSeuratObject(counts = umi),
      vst.flavor = "v2",
      variable.features.n = 100,
      ncells = 150,
      seed.use = 7,
      verbose = FALSE
    )
    for (key in names(branch)) {
      common[key] <- branch[key]
    }
    expected <- suppressWarnings(do.call(seurat_sct, common))
    actual <- suppressWarnings(do.call(SCTransform, common))
    expect_identical(
      SeuratObject::VariableFeatures(actual),
      SeuratObject::VariableFeatures(expected)
    )
    for (layer in c("counts", "data", "scale.data")) {
      expect_equal(
        as.matrix(SeuratObject::LayerData(actual[["SCT"]], layer = layer)),
        as.matrix(SeuratObject::LayerData(expected[["SCT"]], layer = layer)),
        tolerance = 1e-10,
        info = paste(names(branch), collapse = ",")
      )
    }
  }

  regression_specs <- list(
    "numeric_covariate",
    "batch",
    c("numeric_covariate", "batch")
  )
  for (vars in regression_specs) {
    expected_object <- Seurat::CreateSeuratObject(counts = umi)
    actual_object <- Seurat::CreateSeuratObject(counts = umi)
    covariate <- seq_len(ncol(umi)) / ncol(umi)
    batch <- factor(rep(c("a", "b", "c"), length.out = ncol(umi)))
    expected_object$numeric_covariate <- covariate
    actual_object$numeric_covariate <- covariate
    expected_object$batch <- batch
    actual_object$batch <- batch
    args <- list(
      vst.flavor = "v2",
      variable.features.n = 100,
      ncells = 150,
      vars.to.regress = vars,
      seed.use = 31,
      verbose = FALSE
    )
    expected <- suppressWarnings(do.call(
      seurat_sct,
      c(list(object = expected_object), args)
    ))
    actual <- suppressWarnings(do.call(
      SCTransform,
      c(list(object = actual_object), args)
    ))
    expect_identical(
      SeuratObject::VariableFeatures(actual),
      SeuratObject::VariableFeatures(expected),
      info = paste(vars, collapse = ",")
    )
    expect_equal(
      SeuratObject::LayerData(actual[["SCT"]], layer = "scale.data"),
      SeuratObject::LayerData(expected[["SCT"]], layer = "scale.data"),
      tolerance = 1e-9,
      info = paste(vars, collapse = ",")
    )
  }
})

test_that("SCTransform native cores keep residualization identical", {
  skip_if_not_installed("glmGamPoi")
  set.seed(20260903)
  umi <- Matrix::rsparsematrix(
    nrow = 200L,
    ncol = 140L,
    density = 0.04,
    rand.x = function(n) stats::rpois(n, lambda = 4) + 1
  )
  dimnames(umi) <- list(
    paste0("G", seq_len(nrow(umi))),
    paste0("C", seq_len(ncol(umi)))
  )
  umi[1L, ] <- 1
  cell_attr <- data.frame(
    log_umi = log10(Matrix::colSums(umi)),
    row.names = colnames(umi)
  )
  args <- list(
    object = umi,
    cell.attr = cell_attr,
    vst.flavor = "v2",
    variable.features.n = 80,
    ncells = 140,
    seed.use = 11,
    verbose = FALSE
  )
  single <- do.call(
    getS3method("SCTransform", "default"),
    c(args, list(cores = 1L))
  )
  multi <- do.call(
    getS3method("SCTransform", "default"),
    c(args, list(cores = 4L))
  )
  expect_identical(multi$variable_features, single$variable_features)
  expect_equal(multi$y, single$y, tolerance = 1e-12)
  expect_equal(multi$umi_corrected, single$umi_corrected, tolerance = 0)
})

test_that("SCTransform default native path supports latent.data", {
  skip_if_not_installed("glmGamPoi")
  set.seed(20260901)
  umi <- Matrix::rsparsematrix(
    nrow = 180L,
    ncol = 120L,
    density = 0.03,
    rand.x = function(n) stats::rpois(n, lambda = 3) + 1
  )
  dimnames(umi) <- list(paste0("G", seq_len(nrow(umi))), paste0("C", seq_len(ncol(umi))))
  umi[1L, ] <- 1
  cell_attr <- data.frame(
    log_umi = log10(Matrix::colSums(umi)),
    batch = factor(rep(c("a", "b"), length.out = ncol(umi))),
    numeric_covariate = seq_len(ncol(umi)) / ncol(umi),
    row.names = colnames(umi)
  )
  latent <- data.frame(
    batch = cell_attr$batch,
    numeric_covariate = cell_attr$numeric_covariate,
    row.names = colnames(umi)
  )
  args <- list(
    object = umi,
    cell.attr = cell_attr,
    vars.to.regress = c("batch", "numeric_covariate"),
    latent.data = latent,
    vst.flavor = "v2",
    variable.features.n = 80,
    ncells = 120,
    seed.use = 91,
    verbose = FALSE
  )
  expected <- suppressWarnings(do.call(
    get("SCTransform.default", asNamespace("Seurat")),
    args
  ))
  actual <- suppressWarnings(do.call(getS3method("SCTransform", "default"), args))
  expect_identical(actual$variable_features, expected$variable_features)
  expect_equal(actual$y, expected$y, tolerance = 1e-9)
  expect_equal(actual$umi_corrected, expected$umi_corrected, tolerance = 0)
})
