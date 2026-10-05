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
    regression_ok = TRUE,
    defer.residual.matrix = FALSE
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
    regression_ok = FALSE,
    defer.residual.matrix = TRUE
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

test_that("SCTransform adapts its default API across Seurat upgrades and downgrades", {
  adapt <- get("sct_adapt_default_formals", asNamespace("scop"))
  candidate <- get("SCTransform.default", asNamespace("scop"))
  reference <- get("SCTransform.default", asNamespace("Seurat"))
  old_reference <- reference
  formals(old_reference)$defer.residual.matrix <- NULL
  new_reference <- old_reference
  formals(new_reference) <- as.pairlist(append(
    as.list(formals(old_reference)),
    list(defer.residual.matrix = FALSE),
    after = match("return.only.var.genes", names(formals(old_reference)))
  ))

  # Exercise an upgrade, downgrade, and another upgrade of the same method.
  for (ref in list(old_reference, new_reference, old_reference, new_reference)) {
    candidate <- adapt(candidate, ref)
    expect_identical(
      names(formals(candidate)),
      append(names(formals(ref)), "cores", after = length(formals(ref)) - 1L)
    )
    expect_identical(body(candidate), body(get("SCTransform.default", asNamespace("scop"))))
    expect_identical(environment(candidate), asNamespace("scop"))
    if ("defer.residual.matrix" %in% names(formals(ref))) {
      expect_identical(formals(candidate)$defer.residual.matrix, FALSE)
    }
  }
  expect_identical(
    getS3method("SCTransform", "default"),
    get("SCTransform.default", asNamespace("scop"))
  )
})

test_that("SCTransform forwards deferred residuals without leaking cores to Seurat", {
  adapt <- get("sct_adapt_default_formals", asNamespace("scop"))
  upstream <- function(object, cell.attr, defer.residual.matrix = FALSE, ...) {
    list(defer = defer.residual.matrix, dots = list(...))
  }
  candidate <- adapt(get("SCTransform.default", asNamespace("scop")), upstream)
  testthat::local_mocked_bindings(SCTransform.default = upstream, .package = "Seurat")
  umi <- matrix(1, nrow = 2L, ncol = 3L, dimnames = list(c("G1", "G2"), c("C1", "C2", "C3")))
  cell_attr <- data.frame(row.names = colnames(umi))

  # TRUE must itself select fallback, even with no other non-native arguments.
  expect_message(
    deferred <- candidate(umi, cell_attr, defer.residual.matrix = TRUE, cores = 2L),
    "delegating"
  )
  expect_identical(deferred$defer, TRUE)
  expect_false("cores" %in% names(deferred$dots))
  expect_identical(deferred$dots$return.only.var.genes, TRUE)

  # FALSE must be consumed by this wrapper even when another argument delegates.
  expect_message(
    eager <- candidate(
      umi, cell_attr, defer.residual.matrix = FALSE, conserve.memory = TRUE,
      min_cells = 1L, cores = 2L
    ),
    "delegating"
  )
  expect_identical(eager$defer, FALSE)
  expect_identical(eager$dots$conserve.memory, TRUE)
  expect_identical(eager$dots$min_cells, 1L)
  expect_false("cores" %in% names(eager$dots))
})

test_that("SCTransform keeps default and explicit eager residuals on the native path", {
  adapt <- get("sct_adapt_default_formals", asNamespace("scop"))
  upstream <- function(object, cell.attr, defer.residual.matrix = FALSE, ...) {
    stop("unexpected Seurat fallback")
  }
  candidate <- adapt(get("SCTransform.default", asNamespace("scop")), upstream)
  testthat::local_mocked_bindings(SCTransform.default = upstream, .package = "Seurat")
  observed <- NULL
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(package, fun) {
      expect_identical(package, "sctransform")
      expect_identical(fun, "vst")
      function(...) {
        observed <<- list(...)
        stop("native VST reached")
      }
    },
    .package = "scop"
  )
  umi <- matrix(1, nrow = 2L, ncol = 3L, dimnames = list(c("G1", "G2"), c("C1", "C2", "C3")))
  cell_attr <- data.frame(row.names = colnames(umi))
  expect_error(candidate(umi, cell_attr), "native VST reached")
  expect_error(candidate(umi, cell_attr, defer.residual.matrix = FALSE), "native VST reached")
  expect_identical(observed$residual_type, "none")
  expect_false("defer.residual.matrix" %in% names(observed))
  expect_false("cores" %in% names(observed))
})

test_that("SCTransform does not pass the new formal into older Seurat VST arguments", {
  upstream <- function(object, cell.attr, ...) list(...)
  adapt <- get("sct_adapt_default_formals", asNamespace("scop"))
  candidate <- adapt(get("SCTransform.default", asNamespace("scop")), upstream)
  testthat::local_mocked_bindings(SCTransform.default = upstream, .package = "Seurat")
  umi <- matrix(1, nrow = 2L, ncol = 3L, dimnames = list(c("G1", "G2"), c("C1", "C2", "C3")))
  cell_attr <- data.frame(row.names = colnames(umi))
  expect_message(
    out <- candidate(umi, cell_attr, conserve.memory = TRUE, min_cells = 1L),
    "delegating"
  )
  expect_false("defer.residual.matrix" %in% names(out))
  expect_false("cores" %in% names(out))
  expect_identical(out$min_cells, 1L)
})

test_that("SCTransform retains upstream deferred-output and fallback semantics", {
  reference <- get("SCTransform.default", asNamespace("Seurat"))
  candidate <- get("SCTransform.default", asNamespace("scop"))
  if ("defer.residual.matrix" %in% names(formals(reference))) {
    set.seed(107)
    umi <- Matrix::Matrix(
      matrix(stats::rpois(180L * 100L, lambda = 3), nrow = 180L),
      sparse = TRUE
    )
    dimnames(umi) <- list(paste0("G", seq_len(nrow(umi))), paste0("C", seq_len(ncol(umi))))
    cell_attr <- data.frame(row.names = colnames(umi))
    args <- list(
      object = umi, cell.attr = cell_attr, defer.residual.matrix = TRUE,
      variable.features.n = 30L, seed.use = 17L, verbose = FALSE
    )
    for (conserve in c(FALSE, TRUE)) {
      args$conserve.memory <- conserve
      expected <- suppressWarnings(suppressMessages(do.call(reference, args)))
      expect_message(
        actual <- suppressWarnings(do.call(candidate, args)),
        "delegating"
      )
      expect_identical(actual$variable_features, expected$variable_features)
      expect_identical(actual$umi_corrected, expected$umi_corrected)
      expect_equal(actual$y, expected$y, tolerance = 1e-12)
      expect_identical(ncol(actual$y), ncol(umi))
      expect_identical(nrow(actual$y), if (conserve) 30L else 0L)
      expect_identical(colnames(actual$y), colnames(umi))
      expect_equal(actual$model_pars_fit, expected$model_pars_fit, tolerance = 0)
    }
  } else {
    expect_false("defer.residual.matrix" %in% names(formals(candidate)))
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
