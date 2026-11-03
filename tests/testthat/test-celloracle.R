celloracle_test_object <- function() {
  counts <- matrix(rep(c(1, 3, 2, 5, 4, 6), 20), nrow = 3,
    dimnames = list(c("Tf1", "Gene2", "Gene3"), paste0("c", seq_len(40))))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$cluster <- rep(c("a", "b"), each = 20)
  embedding <- cbind(UMAP_1 = seq_len(40), UMAP_2 = sin(seq_len(40)))
  rownames(embedding) <- colnames(object)
  object[["umap"]] <- SeuratObject::CreateDimReducObject(embedding, key = "UMAP_", assay = "RNA")
  object
}

test_that("CellOracle rejects invalid inputs before Python preparation", {
  local_mocked_bindings(PrepareEnv = function(...) stop("Unexpected preparation"))
  object <- celloracle_test_object()
  run <- function(...) RunCellOracle(object, c(Tf1 = 0), "cluster", ...)
  expect_error(run(), "official TF-info")
  expect_error(run(base_grn = data.frame(TF = "Tf1", target = "Gene2", importance = 1)), "official TF-info")
  expect_error(RunCellOracle(object, c(Tf1 = -1), "cluster"), "non-negative")
  expect_error(RunCellOracle(object, c(Tf1 = 0), "missing"), "metadata column")
  expect_error(run(features = c("Gene2", "Gene3")), "perturbed genes")
  expect_error(run(stage_args = list(simulation = list(perturb_condition = list(Tf1 = 1)))), "SCOP-controlled")
  expect_error(run(stage_args = list(transition = list(calculate_randomized = FALSE))), "SCOP-controlled")
  expect_error(run(stage_args = list(grid = list(smooth = 1))), "SCOP-controlled")
  expect_error(run(stage_args = list(imputation = list(k = 2.5))), "must be an integer")
  object@tools$CellOracle <- list(existing = TRUE)
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster"), "new, non-empty")
  expect_true(object@tools$CellOracle$existing)
})

test_that("CellOracle does not silently normalize expression or overwrite models", {
  object <- celloracle_test_object()
  negative <- -GetAssayData5(object, assay = "RNA", layer = "counts")
  object <- SeuratObject::SetAssayData(object, assay = "RNA", layer = "data", new.data = negative)
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster", layer = "data"), "non-negative")
  prior <- list(Gene2 = "Tf1")
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster", base_grn = prior,
    output_dir = tempdir()), "new directory")
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster", model_dir = tempdir()), "saved Oracle/Links")
  local_mocked_bindings(is_windows = function() TRUE,
    PrepareEnv = function(...) stop("Unexpected preparation"))
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster", base_grn = prior), "native Windows")
  object <- SeuratObject::SetAssayData(object, assay = "RNA", layer = "data", new.data = -negative / 2)
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster", layer = "data", base_grn = prior), "native Windows")
})

test_that("CellOracle analysis requires explicit environment preparation", {
  local_mocked_bindings(is_windows = function() FALSE,
    resolve_conda = function(...) "conda",
    env_exist = function(...) FALSE,
    PrepareEnv = function(...) stop("Unexpected installation"))
  object <- celloracle_test_object()
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster",
    base_grn = list(Gene2 = "Tf1")), "explicitly")
  expect_error(RunCellOracle(object, c(Tf1 = 0), "cluster",
    base_grn = list(Gene2 = "Tf1"), docker_container = "bad container name"), "container name or ID")
})

test_that("CellOracle environment is explicit and isolated", {
  expect_false("celloracle" %in% default_env_modules())
  expect_error(normalize_env_modules(c("celloracle", "scanpy")), "standalone")
  requirements <- env_requirements(modules = "celloracle")
  expect_identical(requirements$python, "3.10-1")
  expect_identical(unname(requirements$packages["celloracle"]), "celloracle==0.20.0")
  expect_identical(unname(requirements$packages["numpy"]), "numpy==1.26.4")
  expect_identical(env_module_environment("celloracle"), "celloracle_env")
})

celloracle_test_result <- function(object) {
  cells <- colnames(object)
  genes <- rownames(object)
  list(cells = cells, genes = genes,
    cell_shift = data.frame(cell = cells, group = object$cluster,
      x = cos(seq_along(cells)), y = sin(seq_along(cells)), dx = 0.1, dy = 0.2,
      dx_random = -0.1, dy_random = -0.2),
    grid_flow = data.frame(x = c(1, 2), y = 0, dx = c(3, 9), dy = 1,
      dx_random = c(-3, -9), dy_random = -1, mass = c(1, 0), filtered = c(FALSE, TRUE)),
    gene_effect = data.frame(gene = genes, mean_delta = c(-1, 2, 0), mean_abs_delta = c(1, 2, 0)),
    filtered_links = list(a = data.frame(source = "Tf1", target = "Gene2", coef_abs = 0.2)),
    expression = list(delta_X = matrix(0, nrow = length(cells), ncol = length(genes))))
}

test_that("CellOracle read-back preserves matrix orientation and rejects reordered cells", {
  object <- celloracle_test_object()
  result <- celloracle_test_result(object)
  validated <- celloracle_validate_result(result, colnames(object), rownames(object))
  expect_identical(dimnames(validated$expression$delta_X), list(colnames(object), rownames(object)))
  invalid <- result
  invalid$cell_shift$dx[1] <- NaN
  expect_error(celloracle_validate_result(invalid, colnames(object), rownames(object)), "invalid numerical")
  result$cell_shift <- result$cell_shift[rev(seq_len(ncol(object))), ]
  expect_error(celloracle_validate_result(result, colnames(object), rownames(object)), "different cell or gene order")
})

test_that("CellOracle plots use upstream vectors and exclusion mask without recomputation", {
  object <- celloracle_test_object()
  object@tools$CellOracle <- celloracle_test_result(object)
  flow <- ggplot2::ggplot_build(CellOraclePlot(object, plot_type = "flow", arrow_scale = 2))
  expect_equal(nrow(flow$data[[2]]), 1L)
  expect_equal(flow$data[[2]]$xend - flow$data[[2]]$x, 6)
  control <- ggplot2::ggplot_build(CellOraclePlot(object, plot_type = "flow", randomized = TRUE))
  expect_equal(control$data[[2]]$xend - control$data[[2]]$x, -3)
  quiver <- ggplot2::ggplot_build(CellOraclePlot(object))
  expect_equal(nrow(quiver$data[[2]]), ncol(object))
  effect <- CellOraclePlot(object, plot_type = "effect", top_n = 2)
  expect_equal(as.character(effect$data$gene), c("Gene2", "Tf1"))
  expect_equal(effect$data$mean_abs_delta, c(2, 1))
  expect_s3_class(CellOraclePlot(object, plot_type = "network"), "ggplot")
  object@tools$CellOracle$filtered_links$a <- data.frame(
    source = c("Tf1", "Gene2"), target = c("Gene2", "Tf1"), coef_abs = c(0.2, 0.2))
  network <- ggplot2::ggplot_build(CellOraclePlot(object, plot_type = "network"))
  expect_equal(nrow(network$data[[1]]), 2L)
  expect_equal(network$data[[1]]$x, rev(network$data[[1]]$xend))
  expect_lt(network$data[[1]]$x[1], network$data[[1]]$xend[1])
  object@tools$CellOracle$filtered_links$b <- object@tools$CellOracle$filtered_links$a
  expect_error(CellOraclePlot(object, plot_type = "network"), "Select one")
})
