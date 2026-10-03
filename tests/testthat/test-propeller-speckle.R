make_propeller_reference_srt <- function(paired = FALSE) {
  tab <- rbind(
    c(100, 0, 0, 0), c(45, 30, 20, 5), c(95, 15, 5, 5), c(40, 20, 15, 5),
    c(15, 55, 20, 10), c(10, 60, 20, 0), c(0, 60, 25, 15), c(20, 40, 35, 5)
  )
  dat <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
    data.frame(CellType = rep(c("T", "B", "Mono", "Rare"), tab[i, ]),
      Condition = if (i <= 4L) "A" else "B",
      Sample = if (paired) paste0("donor", (i - 1L) %% 4L + 1L) else paste0("sample", i))
  }))
  counts <- Matrix::sparseMatrix(i = rep(1L, nrow(dat)), j = seq_len(nrow(dat)), x = 1,
    dims = c(2L, nrow(dat)),
    dimnames = list(c("gene1", "gene2"), paste0("cell", seq_len(nrow(dat)))))
  srt <- SeuratObject::CreateSeuratObject(counts)
  rownames(dat) <- colnames(srt)
  SeuratObject::AddMetaData(srt, dat)
}

test_that("SCOP propeller matches upstream with zeros, ones, and unequal totals", {
  skip_if_not_installed("speckle")
  check <- getFromNamespace("propeller_check_r", "scop")
  check()
  upstream <- thisutils::get_namespace_fun("speckle", "propeller")
  srt <- make_propeller_reference_srt()
  for (transform in c("logit", "asin")) {
    for (robust in c(FALSE, TRUE)) {
      for (trend in c(FALSE, TRUE)) {
        bundle <- RunPropeller(srt, "CellType", "Condition", "Sample",
          comparison = list(c("A", "B")), n_bootstrap = 0L,
          transform = transform, robust = robust, trend = trend, verbose = FALSE)
        ref <- upstream(clusters = factor(srt$CellType), sample = factor(srt$Sample),
          group = factor(srt$Condition, levels = c("A", "B")),
          transform = transform, robust = robust, trend = trend)
        result <- bundle$results$A_vs_B
        ref <- ref[match(result$clusters, rownames(ref)), , drop = FALSE]
        expect_equal(result$pval, ref$P.Value, tolerance = 1e-12)
        expect_equal(result$FDR, ref$FDR, tolerance = 1e-12)
        expect_equal(result$Tstatistic, ref$Tstatistic, tolerance = 1e-12)
        expect_equal(result$PropRatio, ref$PropRatio, tolerance = 1e-12)
        expect_identical(bundle$details$engine, "speckle")
        expect_identical(bundle$parameters$transform, transform)
        expect_true(all(is.na(result$boot_CI_2.5)))
      }
    }
  }
})

test_that("speckle comparison reversal preserves significance and flips direction", {
  skip_if_not_installed("speckle")
  srt <- make_propeller_reference_srt()
  out <- RunProportionTest(srt, group.by = "CellType", split.by = "Condition",
    sample.by = "Sample", proportion_method = "propeller", n_bootstrap = 0L,
    verbose = FALSE)
  result <- out@tools$ProportionTest$results
  a <- result$A_vs_B
  b <- result$B_vs_A[match(a$clusters, result$B_vs_A$clusters), , drop = FALSE]
  expect_equal(a$pval, b$pval, tolerance = 1e-12)
  expect_equal(a$FDR, b$FDR, tolerance = 1e-12)
  expect_equal(a$Tstatistic, -b$Tstatistic, tolerance = 1e-12)
  expect_equal(a$obs_log2FD, -b$obs_log2FD, tolerance = 1e-12)
  expect_gt(a$obs_log2FD[a$clusters == "T"], 0)
  tab <- table(srt$CellType, srt$Sample)
  props <- sweep(tab, 2, colSums(tab), "/")
  expected <- log2((rowMeans(props[, paste0("sample", 1:4)]) + 1e-5) /
    (rowMeans(props[, paste0("sample", 5:8)]) + 1e-5))
  expect_equal(a$obs_log2FD, unname(expected[a$clusters]))
  plots <- ProportionTestPlot(out, proportion_method = "propeller", combine = FALSE)
  expect_match(plots[[1]]$labels$subtitle, "Speckle propeller")
  # Previously saved results still identify the original internal test.
  out@tools$ProportionTest$methods$propeller$results <- lapply(result, function(x) {
    x$engine <- NULL
    x
  })
  plots <- ProportionTestPlot(out, proportion_method = "propeller", combine = FALSE)
  expect_match(plots[[1]]$labels$subtitle, "Internal logit")
})

test_that("independent sample bootstrap is reproducible and leaves speckle tests unchanged", {
  skip_if_not_installed("speckle")
  srt <- make_propeller_reference_srt()
  run <- function(n) RunPropeller(srt, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = n,
    seed = 19L, verbose = FALSE)$results$A_vs_B
  result <- run(30L)
  expect_true(all(is.finite(result$boot_CI_2.5)))
  expect_true(all(is.finite(result$boot_CI_97.5)))
  expect_identical(result$boot_CI_2.5, run(30L)$boot_CI_2.5)
  expect_identical(result$pval, run(0L)$pval)
  expect_identical(result$FDR, run(0L)$FDR)
})

test_that("paired propeller matches an upstream donor-blocked design", {
  skip_if_not_installed("speckle")
  srt <- make_propeller_reference_srt(paired = TRUE)
  bundle <- RunPropeller(srt, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = 30L, verbose = FALSE)
  backend <- bundle$details$backend_results$A_vs_B
  expect_true(backend$paired)
  expect_equal(nrow(backend$design), 8L)
  expect_equal(ncol(backend$design), 5L)
  check <- getFromNamespace("propeller_check_r", "scop")
  check()
  props <- thisutils::get_namespace_fun("speckle", "getTransformedProps")(
    clusters = factor(srt$CellType),
    sample = factor(paste(srt$Sample, srt$Condition, sep = ":")), transform = "logit")
  keys <- colnames(props$Proportions)
  condition <- factor(sub(".*:", "", keys), levels = c("A", "B"))
  donor <- factor(sub(":.*", "", keys))
  design <- stats::model.matrix(~ 0 + condition + donor)
  ref <- thisutils::get_namespace_fun("speckle", "propeller.ttest")(
    props, design, c(1, -1, 0, 0, 0), robust = TRUE, trend = FALSE, sort = FALSE)
  result <- bundle$results$A_vs_B
  ref <- ref[match(result$clusters, rownames(ref)), , drop = FALSE]
  expect_equal(result$pval, ref$P.Value, tolerance = 1e-12)
  expect_equal(result$FDR, ref$FDR, tolerance = 1e-12)
  expect_equal(result$Tstatistic, ref$Tstatistic, tolerance = 1e-12)
  expect_true(all(is.finite(result$boot_CI_2.5)))
  again <- RunPropeller(srt, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = 30L, verbose = FALSE)
  expect_identical(result$boot_CI_2.5, again$results$A_vs_B$boot_CI_2.5)
  shuffled <- srt[, rev(seq_len(ncol(srt))) ]
  again <- RunPropeller(shuffled, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = 30L, verbose = FALSE)
  expect_equal(result, again$results$A_vs_B)
})

test_that("propeller rejects insufficient replication and partial pairing", {
  skip_if_not_installed("speckle")
  srt <- make_propeller_reference_srt(paired = TRUE)
  partial <- srt[, !(srt$Sample == "donor4" & srt$Condition == "B")]
  expect_error(RunPropeller(partial, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = 0L), "fully paired")
  srt <- make_propeller_reference_srt()
  single <- srt[, srt$Sample %in% c("sample1", "sample5", "sample6")]
  expect_error(RunPropeller(single, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = 0L), "at least two biological")
  srt$CellType[1] <- NA_character_
  expect_error(RunPropeller(srt, "CellType", "Condition", "Sample",
    comparison = list(c("A", "B")), n_bootstrap = 0L), "cell type labels")
})

test_that("propeller reports unavailable speckle without substituting another test", {
  testthat::local_mocked_bindings(
    check_r = function(...) list(speckle = FALSE), .package = "scop")
  srt <- make_propeller_reference_srt()
  expect_error(RunPropeller(srt, "CellType", "Condition", "Sample",
    n_bootstrap = 0L), "requires the optional Bioconductor package 'speckle'")
})

test_that("virtual propeller samples mask upstream statistics too", {
  skip_if_not_installed("speckle")
  srt <- make_propeller_reference_srt()
  out <- suppressWarnings(RunProportionTest(srt, group.by = "CellType",
    split.by = "Condition", proportion_method = "propeller",
    allow_pseudo_samples = TRUE, n_bootstrap = 0L, verbose = FALSE))
  bundle <- out@tools$ProportionTest$methods$propeller
  expect_false(bundle$inference_valid)
  expect_true(all(is.na(bundle$results$A_vs_B$pval)))
  expect_true(all(is.na(bundle$results$A_vs_B$Tstatistic)))
  expect_true(all(is.na(bundle$details$backend_results$A_vs_B$results$P.Value)))
})

make_propeller_robust_srt <- function() {
  # Fixed counts avoid RNG/version dependence. Most types vary mildly, while
  # the last two alternate between rare and abundant within both conditions.
  tab <- outer(seq_len(12L), seq_len(30L), function(i, j) {
    40L + 2L * j + ((i * (j + 3L) + i^2) %% 13L)
  })
  tab[7:12, 1:6] <- tab[7:12, 1:6] + 25L
  tab[, 29L] <- rep(c(1L, 400L, 5L, 300L, 2L, 500L), 2L)
  tab[, 30L] <- rep(c(350L, 2L, 450L, 1L, 250L, 3L), 2L)
  dat <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
    data.frame(CellType = rep(sprintf("type%02d", seq_len(ncol(tab))), tab[i, ]),
      Condition = if (i <= 6L) "A" else "B", Sample = sprintf("sample%02d", i))
  }))
  counts <- Matrix::sparseMatrix(i = rep(1L, nrow(dat)), j = seq_len(nrow(dat)), x = 1,
    dims = c(2L, nrow(dat)),
    dimnames = list(c("gene1", "gene2"), paste0("cell", seq_len(nrow(dat)))))
  srt <- SeuratObject::CreateSeuratObject(counts)
  rownames(dat) <- colnames(srt)
  SeuratObject::AddMetaData(srt, dat)
}

test_that("propeller exercises limma robust estimation with variance outliers", {
  skip_if_not_installed("speckle")
  getFromNamespace("propeller_check_r", "scop")()
  srt <- make_propeller_robust_srt()
  expect_equal(length(unique(srt$CellType)), 30L)
  real_fit <- getFromNamespace("fitFDistRobustly", "limma")
  fits <- list()
  testthat::local_mocked_bindings(
    fitFDistRobustly = function(...) {
      # Delegate unchanged to limma and retain its diagnostics. Small-input
      # non-robust fallback has no tail.p.value, even when robust=TRUE.
      fit <- real_fit(...)
      fits[[length(fits) + 1L]] <<- fit
      fit
    }, .package = "limma")
  upstream <- thisutils::get_namespace_fun("speckle", "propeller")
  for (transform in c("logit", "asin")) {
    for (trend in c(FALSE, TRUE)) {
      fits <- list()
      result <- RunPropeller(srt, "CellType", "Condition", "Sample",
        comparison = list(c("A", "B")), n_bootstrap = 0L,
        transform = transform, robust = TRUE, trend = trend,
        verbose = FALSE)$results$A_vs_B
      # RunPropeller also evaluates the reverse comparison.
      expect_length(fits, 2L)
      fit <- fits[[1L]]
      expect_length(fit$tail.p.value, 30L)
      expect_true(all(is.finite(fit$tail.p.value)))
      expect_true(any(fit$tail.p.value < 0.01))
      expect_gt(length(unique(fit$df2.shrunk)), 1L)
      ref <- upstream(clusters = factor(srt$CellType), sample = factor(srt$Sample),
        group = factor(srt$Condition, levels = c("A", "B")),
        transform = transform, robust = TRUE, trend = trend)
      ref <- ref[match(result$clusters, rownames(ref)), , drop = FALSE]
      expect_equal(result$pval, ref$P.Value, tolerance = 1e-12)
      expect_equal(result$FDR, ref$FDR, tolerance = 1e-12)
      expect_equal(result$Tstatistic, ref$Tstatistic, tolerance = 1e-12)
      nonrobust <- RunPropeller(srt, "CellType", "Condition", "Sample",
        comparison = list(c("A", "B")), n_bootstrap = 0L,
        transform = transform, robust = FALSE, trend = trend,
        verbose = FALSE)$results$A_vs_B
      expect_gt(max(abs(result$pval - nonrobust$pval)), 1e-6)
    }
  }
})
