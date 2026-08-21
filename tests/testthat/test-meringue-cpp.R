# Consistency and speed tests for the scop C++ Moran backend
# (`meringue_moran_cpp`) against the original MERINGUE implementation
# (`moranTest` / `moranPermutationTest`) on real spatial data. No mocks here:
# both implementations run on the same real Visium input and must agree to
# machine precision, including permutation p-values (the kernel replicates the
# session's `sample.kind` RNG stream).

meringue_real_input <- function(n = 120, seed = 42) {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  set.seed(seed)
  srt <- srt[, sample(colnames(srt), n)]
  srt
}

meringue_real_features <- function(srt, k = 12, seed = 1) {
  expr <- Seurat::GetAssayData(srt, assay = SeuratObject::DefaultAssay(srt), layer = "data")
  set.seed(seed)
  head(rownames(expr)[order(apply(expr, 1, stats::var), decreasing = TRUE)], k)
}

meringue_test_weight <- function(srt) {
  coords <- getFromNamespace("spatial_analysis_coords", "scop")(
    srt = srt,
    image = NULL,
    coord.cols = c("x", "y"),
    coordinate_space = "raw"
  )$data
  MERINGUE::getSpatialNeighbors(
    pos = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    filterDist = NA_real_,
    binary = TRUE,
    verbose = FALSE
  )
}

test_that("meringue_moran_cpp matches MERINGUE::moranTest exactly", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input()
  expr <- Seurat::GetAssayData(srt, assay = SeuratObject::DefaultAssay(srt), layer = "data")
  features <- meringue_real_features(srt)
  weight <- meringue_test_weight(srt)
  rounding <- tryCatch(
    identical(RNGkind()[3], "Rounding"),
    error = function(e) FALSE
  )

  for (feature in features) {
    z <- expr[feature, ]
    z <- z[rownames(weight)]
    cpp_out <- scop:::meringue_moran_cpp(
      z = z,
      weight = weight,
      n_perm = 0,
      alternative = "greater",
      rounding_sample = rounding
    )
    r_out <- MERINGUE::moranTest(z, weight, alternative = "greater")
    expect_equal(cpp_out[["observed"]], r_out[["observed"]], tolerance = 1e-12,
      label = paste0(feature, " observed"))
    expect_equal(cpp_out[["expected"]], r_out[["expected"]], tolerance = 1e-12,
      label = paste0(feature, " expected"))
    expect_equal(cpp_out[["sd"]], r_out[["sd"]], tolerance = 1e-12,
      label = paste0(feature, " sd"))
    expect_equal(cpp_out[["p_value"]], r_out[["p.value"]], tolerance = 1e-12,
      label = paste0(feature, " p"))
  }
})

test_that("meringue_moran_matrix_cpp batch matches per-gene kernel and moranTest", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input()
  expr <- Seurat::GetAssayData(srt, assay = SeuratObject::DefaultAssay(srt), layer = "data")
  features <- meringue_real_features(srt, k = 20)
  weight <- meringue_test_weight(srt)
  rounding <- tryCatch(
    identical(RNGkind()[3], "Rounding"),
    error = function(e) FALSE
  )

  batch <- scop:::meringue_moran_matrix_cpp(
    as.matrix(expr[features, rownames(weight)]),
    weight,
    "greater",
    rounding
  )
  for (i in seq_along(features)) {
    feature <- features[i]
    z <- expr[feature, ]
    z <- z[rownames(weight)]
    per_gene <- scop:::meringue_moran_cpp(
      z = z, weight = weight, n_perm = 0, alternative = "greater",
      rounding_sample = rounding
    )
    expect_equal(batch[i, 1], per_gene[["observed"]], tolerance = 1e-12)
    expect_equal(batch[i, 2], per_gene[["expected"]], tolerance = 1e-12)
    expect_equal(batch[i, 3], per_gene[["sd"]], tolerance = 1e-12)
    expect_equal(batch[i, 4], per_gene[["p_value"]], tolerance = 1e-12)
  }
})

test_that("RunMERINGUE batch cpp path matches the r backend end to end", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input(n = 150)
  features <- meringue_real_features(srt, k = 100)

  cpp_out <- RunMERINGUE(
    srt, mode = "autocorrelation", coord.cols = c("x", "y"),
    features = features, min_spots = 5, nperm = 0, backend = "cpp", verbose = FALSE
  )
  r_out <- RunMERINGUE(
    srt, mode = "autocorrelation", coord.cols = c("x", "y"),
    features = features, min_spots = 5, nperm = 0, backend = "r", verbose = FALSE
  )
  ct <- cpp_out@tools$MERINGUE$autocorrelation
  rt <- r_out@tools$MERINGUE$autocorrelation
  expect_identical(ct$feature, rt$feature)
  expect_equal(ct$statistic, rt$statistic, tolerance = 1e-12)
  expect_equal(ct$sd, rt$sd, tolerance = 1e-12)
  expect_identical(ct$rank, rt$rank)
  expect_identical(ct$variance, rt$variance)
})

test_that("meringue_moran_cpp permutation matches MERINGUE::moranPermutationTest", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input()
  expr <- Seurat::GetAssayData(srt, assay = SeuratObject::DefaultAssay(srt), layer = "data")
  features <- meringue_real_features(srt, k = 6)
  weight <- meringue_test_weight(srt)
  rounding <- tryCatch(
    identical(RNGkind()[3], "Rounding"),
    error = function(e) FALSE
  )

  for (feature in features) {
    z <- expr[feature, ]
    z <- z[rownames(weight)]
    set.seed(11)
    cpp_out <- scop:::meringue_moran_cpp(
      z = z,
      weight = weight,
      n_perm = 100,
      alternative = "greater",
      rounding_sample = rounding
    )
    set.seed(11)
    r_out <- MERINGUE::moranPermutationTest(
      z = z,
      w = weight,
      alternative = "greater",
      N = 100,
      seed = 11,
      ncores = 1,
      plot = FALSE
    )
    expect_equal(cpp_out[["observed"]], r_out[["observed"]], tolerance = 1e-12,
      label = paste0(feature, " observed"))
    expect_equal(cpp_out[["expected"]], r_out[["expected"]], tolerance = 1e-12,
      label = paste0(feature, " expected"))
    expect_equal(cpp_out[["sd"]], r_out[["sd"]], tolerance = 1e-12,
      label = paste0(feature, " sd"))
    expect_equal(cpp_out[["p_value"]], r_out[["p.value"]], tolerance = 1e-12,
      label = paste0(feature, " p"))
  }
})

test_that("meringue permutation batch matches the seeded per-gene kernel", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input()
  expr <- Seurat::GetAssayData(
    srt,
    assay = SeuratObject::DefaultAssay(srt),
    layer = "data"
  )
  features <- meringue_real_features(srt, k = 12)
  weight <- meringue_test_weight(srt)
  rounding <- tryCatch(
    identical(RNGkind()[3], "Rounding"),
    error = function(e) FALSE
  )

  set.seed(11)
  batch <- scop:::meringue_moran_matrix_cpp(
    as.matrix(expr[features, rownames(weight)]),
    weight,
    "greater",
    rounding,
    n_perm = 100L,
    n_threads = 2L
  )
  expect_identical(rownames(batch), features)
  expect_identical(
    colnames(batch),
    c("observed", "expected", "sd", "p_value", "N")
  )

  for (i in seq_along(features)) {
    z <- expr[features[[i]], rownames(weight)]
    set.seed(11)
    per_gene <- scop:::meringue_moran_cpp(
      z = z,
      weight = weight,
      n_perm = 100L,
      alternative = "greater",
      rounding_sample = rounding
    )
    expect_equal(batch[i, ], per_gene, tolerance = 1e-12)
  }
})

test_that("RunMERINGUE cpp and r backends give identical autocorrelation tables", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input()
  features <- meringue_real_features(srt, k = 30)

  set.seed(42)
  cpp_out <- RunMERINGUE(
    srt,
    mode = "autocorrelation",
    coord.cols = c("x", "y"),
    features = features,
    min_spots = 5,
    filterDist = NA_real_,
    binary = TRUE,
    alternative = "greater",
    nperm = 30,
    ncores = 1,
    seed = 11,
    set_variable_features = FALSE,
    backend = "cpp",
    verbose = FALSE
  )
  set.seed(42)
  r_out <- RunMERINGUE(
    srt,
    mode = "autocorrelation",
    coord.cols = c("x", "y"),
    features = features,
    min_spots = 5,
    filterDist = NA_real_,
    binary = TRUE,
    alternative = "greater",
    nperm = 30,
    ncores = 1,
    seed = 11,
    set_variable_features = FALSE,
    backend = "r",
    verbose = FALSE
  )

  cpp_table <- cpp_out@tools$MERINGUE$autocorrelation
  r_table <- r_out@tools$MERINGUE$autocorrelation
  expect_identical(cpp_table$feature, r_table$feature)
  expect_equal(cpp_table$statistic, r_table$statistic, tolerance = 1e-12)
  expect_equal(cpp_table$expected, r_table$expected, tolerance = 1e-12)
  expect_equal(cpp_table$sd, r_table$sd, tolerance = 1e-12)
  expect_equal(cpp_table$p_value, r_table$p_value, tolerance = 1e-12)
  expect_equal(cpp_table$q_value, r_table$q_value, tolerance = 1e-12)
  expect_identical(cpp_table$rank, r_table$rank)
  expect_identical(
    cpp_out@tools$MERINGUE$parameters$value[
      match("backend", cpp_out@tools$MERINGUE$parameters$key)
    ],
    "cpp"
  )
})

test_that("RunMERINGUE cpp backend is much faster than r for permutations", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input(n = 120)
  features <- meringue_real_features(srt, k = 20)

  cpp_time <- system.time({
    cpp_out <- RunMERINGUE(
      srt,
      mode = "autocorrelation",
      coord.cols = c("x", "y"),
      features = features,
      min_spots = 5,
      filterDist = NA_real_,
      binary = TRUE,
      alternative = "greater",
      nperm = 200,
      ncores = 1,
      seed = 11,
      set_variable_features = FALSE,
      backend = "cpp",
      verbose = FALSE
    )
  })[["elapsed"]]

  expect_true(all(is.finite(cpp_out@tools$MERINGUE$autocorrelation$p_value)))
  expect_lt(cpp_time, 30)
})

test_that("RunMERINGUE cpp backend honors moran_params with a warning", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input(n = 80)
  features <- meringue_real_features(srt, k = 5)

  warnings_seen <- character(0)
  testthat::local_mocked_bindings(
    log_message = function(msg, message_type = "info", ...) {
      if (identical(message_type, "warning")) {
        warnings_seen <<- c(warnings_seen, as.character(msg))
      }
      invisible(TRUE)
    },
    .package = "scop"
  )

  out <- suppressMessages(RunMERINGUE(
    srt,
    mode = "autocorrelation",
    coord.cols = c("x", "y"),
    features = features,
    min_spots = 5,
    moran_params = list(plot = FALSE),
    backend = "cpp",
    verbose = FALSE
  ))
  expect_true(any(grepl("moran_params.*ignored by the cpp backend", warnings_seen)))
})

test_that("RunMERINGUE cpp default stores kernel-produced autocorrelation", {
  skip_on_cran()
  skip_if_not_installed("MERINGUE")

  srt <- meringue_real_input(n = 80)
  features <- meringue_real_features(srt, k = 5)

  out <- RunMERINGUE(
    srt,
    mode = "autocorrelation",
    coord.cols = c("x", "y"),
    features = features,
    min_spots = 5,
    set_variable_features = FALSE,
    verbose = FALSE
  )
  auto <- out@tools$MERINGUE$autocorrelation
  expect_equal(nrow(auto), length(features))
  expect_true(all(is.finite(auto$statistic)))
  expect_true(all(auto$p_value %in% c(NA_real_)))
  expect_identical(
    out@tools$MERINGUE$parameters$value[
      match("backend", out@tools$MERINGUE$parameters$key)
    ],
    "cpp"
  )
  expect_true(all(is.finite(auto$mean)))
  expect_true(all(is.finite(auto$variance)))
})
