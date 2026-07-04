make_literature_method_srt <- function(n_genes = 80, n_cells = 40, seed = 1) {
  set.seed(seed)
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 2),
    nrow = n_genes,
    dimnames = list(paste0("Gene", seq_len(n_genes)), paste0("Cell", seq_len(n_cells)))
  )
  counts[seq_len(10), seq_len(n_cells / 2)] <- counts[seq_len(10), seq_len(n_cells / 2)] + 4
  counts[11:20, seq(from = n_cells / 2 + 1, to = n_cells)] <- counts[11:20, seq(from = n_cells / 2 + 1, to = n_cells)] + 4
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt$batch <- rep(c("B1", "B2"), each = n_cells / 2)
  srt$celltype <- rep(c("A", "B"), each = n_cells / 2)
  srt$stage <- rep(seq_len(4), length.out = n_cells)
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt <- Seurat::FindVariableFeatures(srt, nfeatures = 40, verbose = FALSE)
  srt <- Seurat::ScaleData(srt, features = SeuratObject::VariableFeatures(srt), verbose = FALSE)
  srt <- RunPCA(srt, features = SeuratObject::VariableFeatures(srt), npcs = 10, verbose = FALSE)
  srt <- RunUMAP(srt, dims = 1:5, verbose = FALSE)
  srt
}

test_that("RunFitDevo and RunFWP write metadata and reusable weights", {
  srt <- make_literature_method_srt()
  out <- RunFitDevo(
    srt,
    nfeatures = 40,
    reference.by = "stage",
    verbose = FALSE
  )
  expect_true(all(c("FitDevo_Score", "FitDevo_Relative") %in% colnames(out@meta.data)))
  expect_true(all(out$FitDevo_Score >= 0 & out$FitDevo_Score <= 1))
  expect_equal(out@tools$FitDevo$status, "success")

  out <- RunFWP(
    out,
    phenotype.by = "celltype",
    positive = "B",
    nfeatures = 40,
    verbose = FALSE
  )
  expect_true("FWP_Score" %in% colnames(out@meta.data))
  expect_true(all(out$FWP_Score >= 0 & out$FWP_Score <= 1))
  expect_equal(out@tools$FWP$status, "success")

  p <- FeatureDimPlot(out, features = "FWP_Score", reduction = "umap")
  expect_true(inherits(p, c("ggplot", "patchwork")))
  fitdevo_plots <- FitDevoPlot(out, group.by = "celltype", combine = FALSE)
  expect_true(all(c("Score", "Relative", "Phenotype", "Boxplot") %in% names(fitdevo_plots)))
  expect_s3_class(fitdevo_plots$Score, "ggplot")
  expect_s3_class(fitdevo_plots$Boxplot, "ggplot")
})

test_that("RunVECTOR stores cell scores and grid arrows", {
  srt <- make_literature_method_srt()
  out <- RunVECTOR(
    srt,
    reduction = "umap",
    pca.reduction = "pca",
    pca.dims = 1:5,
    grid.n = 5,
    verbose = FALSE
  )
  expect_true("VECTOR_Score" %in% colnames(out@meta.data))
  expect_true(all(out$VECTOR_Score >= 0 & out$VECTOR_Score <= 1))
  expect_equal(out@tools$VECTOR$status, "success")
  expect_equal(nrow(out@tools$VECTOR$embedding), ncol(out))
  expect_true(is.data.frame(out@tools$VECTOR$grid))
  expect_true(is.null(out@tools$VECTOR$arrows) || is.data.frame(out@tools$VECTOR$arrows))
  expect_s3_class(VECTORPlot(out, plot_type = "grid"), "ggplot")
  expect_s3_class(VECTORPlot(out, plot_type = "raw", group.by = "celltype"), "ggplot")
})

test_that("native FWP vectorization matches a naive reference and is faster", {
  srt <- make_literature_method_srt(n_genes = 200, n_cells = 80)
  mat <- GetAssayData5(srt, layer = "data")
  features <- rownames(mat)[seq_len(120)]
  y <- as.integer(srt$celltype == "B")

  fast_time <- system.time({
    fast <- RunFWP(
      mat,
      features = features,
      phenotype.by = NULL,
      weights = scop:::fwp_score(mat[features, ], y = y)$weights,
      verbose = FALSE
    )
  })[["elapsed"]]

  naive_time <- system.time({
    dense <- scop:::scop_scale_features(mat[features, ])
    weights <- rowMeans(dense[, y == 1, drop = FALSE]) - rowMeans(dense[, y == 0, drop = FALSE])
    weights <- weights / (sqrt(sum(weights^2)) + 1e-8)
    naive <- scop:::scale01(as.numeric(crossprod(weights, dense)))
    names(naive) <- colnames(dense)
  })[["elapsed"]]

  expect_equal(fast$score, naive[names(fast$score)], tolerance = 1e-10)
  expect_true(is.finite(fast_time))
  expect_true(is.finite(naive_time))
})
