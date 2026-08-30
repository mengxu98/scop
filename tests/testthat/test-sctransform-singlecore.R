test_that("single-core glmGamPoi offset kernel preserves fitted values", {
  skip_if_not_installed("glmGamPoi")
  available <- get("sct_single_offset_available", asNamespace("scop"))
  skip_if(!available())
  set.seed(20260829)
  genes <- paste0("G", seq_len(80L))
  cells <- paste0("C", seq_len(60L))
  umi <- matrix(
    stats::rpois(length(genes) * length(cells), lambda = 1.5),
    nrow = length(genes),
    dimnames = list(genes, cells)
  )
  cell_attr <- data.frame(
    log_umi = log10(colSums(umi)),
    row.names = cells
  )
  reference <- get("fit_glmGamPoi_offset", asNamespace("sctransform"))(
    umi = umi,
    model_str = "y ~ log_umi",
    data = cell_attr,
    allow_inf_theta = TRUE
  )
  candidate <- get(
    "sct_fit_glmgampoi_offset_single",
    asNamespace("scop")
  )(
    umi = umi,
    data = cell_attr
  )
  expect_identical(dimnames(candidate), dimnames(reference))
  expect_equal(
    1 / candidate[, "theta"],
    1 / reference[, "theta"],
    tolerance = 2e-6
  )
  expect_equal(
    candidate[, c("(Intercept)", "log_umi")],
    reference[, c("(Intercept)", "log_umi")],
    tolerance = 1e-7
  )

  offset <- matrix(
    rep(log(10^cell_attr$log_umi), each = nrow(umi)),
    nrow = nrow(umi)
  )
  mu <- exp(candidate[, "(Intercept)"]) * exp(offset)
  reference_residual <- (umi - mu) /
    sqrt(mu + mu^2 / reference[, "theta"])
  candidate_residual <- (umi - mu) /
    sqrt(mu + mu^2 / candidate[, "theta"])
  expect_lt(max(abs(candidate_residual - reference_residual)), 1e-5)
})

test_that("single-core offset kernel is guarded by validated dependencies", {
  skip_if_not_installed("glmGamPoi")
  skip_if_not_installed("sctransform")
  available <- get("sct_single_offset_available", asNamespace("scop"))
  validated_versions <-
    as.character(packageVersion("glmGamPoi")) %in% c("1.18.0", "1.22.0") &&
      packageVersion("sctransform") == "0.4.3"
  expect_identical(available(), validated_versions)
})

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

test_that("SCTransform native path preserves Seurat output semantics", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glmGamPoi")
  available <- get("sct_single_offset_available", asNamespace("scop"))
  skip_if(!available())

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
})
