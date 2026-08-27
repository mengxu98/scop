genes <- c(paste0("MT-G", seq_len(5)), paste0("G", seq_len(195)))
cells <- paste0("C", seq_len(120))
counts <- local({
  set.seed(11)
  m <- matrix(
    rpois(length(genes) * length(cells), lambda = 2),
    nrow = length(genes),
    dimnames = list(genes, cells)
  )
  Matrix::Matrix(m, sparse = TRUE)
})
srt <- Seurat::CreateSeuratObject(counts = counts)
srt$percent.mt <- Matrix::colMeans(counts[seq_len(5), ]) /
  Matrix::colSums(counts)
srt$group <- rep(c("A", "B"), length.out = length(cells))

test_that("SCTransform fast path handles default arguments", {
  skip_if_not_installed("glmGamPoi")
  out <- suppressMessages(SCTransform(srt, seed.use = 42))
  expect_s4_class(out[["SCT"]], "SCTAssay")
  expect_true(length(SeuratObject::VariableFeatures(out)) > 0L)
  expect_identical(ncol(out[["SCT"]]), length(cells))
})

test_that("SCTransform regresses variables on the fast path", {
  skip_if_not_installed("glmGamPoi")
  out <- SCTransform(srt, vars.to.regress = "percent.mt", seed.use = 42)
  scale.data <- SeuratObject::GetAssayData(
    out[["SCT"]],
    layer = "scale.data"
  )
  expect_s4_class(out[["SCT"]], "SCTAssay")
  expect_true(max(abs(rowMeans(scale.data))) < 1e-10)
  model <- methods::slot(out[["SCT"]], "SCTModel.list")[[1]]
  expect_identical(model@arguments$sct.latent.vars, "percent.mt")
})

test_that("SCTransform regresses factor covariates on the fast path", {
  skip_if_not_installed("glmGamPoi")
  out <- suppressMessages(SCTransform(
    srt,
    vars.to.regress = c("percent.mt", "group"),
    seed.use = 42
  ))
  expect_s4_class(out[["SCT"]], "SCTAssay")
})

test_that("SCTransform delegates truly unsupported arguments to Seurat", {
  expect_message(
    out <- suppressWarnings(SCTransform(
      srt,
      conserve.memory = TRUE,
      seed.use = 42
    )),
    regexp = "delegating"
  )
  expect_s4_class(out[["SCT"]], "SCTAssay")

  expect_message(
    out2 <- SCTransform(
      srt,
      vst.flavor = "v1",
      seed.use = 42
    ),
    regexp = "delegating"
  )
  expect_s4_class(out2[["SCT"]], "SCTAssay")
})

test_that("SCTransform.default accepts latent.data regression", {
  skip_if_not_installed("glmGamPoi")
  cell.attr <- data.frame(row.names = cells)
  latent <- data.frame(score = rnorm(length(cells)), row.names = cells)
  out <- suppressMessages(SCTransform.default(
    object = counts,
    cell.attr = cell.attr,
    latent.data = latent,
    seed.use = 42
  ))
  expect_identical(out$arguments$sct.latent.vars, "score")
  expect_true(max(abs(rowMeans(out$y))) < 1e-10)
})
