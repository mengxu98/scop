make_commot_test_object <- function() {
  counts <- Matrix::sparseMatrix(
    i = rep(seq_len(4), 4), j = rep(seq_len(4), each = 4),
    x = rep(c(5, 3, 2, 1), 4), dims = c(4, 4),
    dimnames = list(c("L1", "R1", "G1", "G2"), paste0("cell", 1:4))
  )
  srt <- SeuratObject::CreateSeuratObject(counts = counts)
  srt$celltype <- c("A", "A", "B", "B")
  srt$col <- c(0, 1, 2, 3)
  srt$row <- c(0, 1, 1, 2)
  srt
}

mock_commot_result <- function(input, cluster = FALSE, direction = FALSE) {
  communication <- data.frame(
    sender = "A", receiver = "B", ligand = "L1", receptor = "R1",
    interaction_name = "L1_R1", pathway_name = "P1", score = 2.5,
    pvalue = NA_real_, method = "COMMOT", score_type = "transport_mass",
    stringsAsFactors = FALSE
  )
  cluster_table <- if (isTRUE(cluster)) data.frame(
    sender = c("A", "B"), receiver = c("B", "A"),
    ligand = NA_character_, receptor = NA_character_,
    interaction_name = NA_character_, pathway_name = "P1",
    score = c(3, 1), pvalue = c(0.01, 0.2), method = "COMMOT",
    score_type = "cluster_communication", key = "CellChat-P1",
    stringsAsFactors = FALSE
  ) else data.frame()
  direction_table <- if (isTRUE(direction)) do.call(rbind, lapply(
    c("sender", "receiver"),
    function(perspective) data.frame(
      cell_id = input$cells, group = input$groups, key = "CellChat-P1",
      perspective = perspective,
      x = input$coordinates$x, y = input$coordinates$y,
      dx = c(0.1, 0.2, 0.1, 0), dy = c(0, 0.1, 0.2, 0.1),
      magnitude = c(0.1, sqrt(0.05), sqrt(0.05), 0.1),
      stringsAsFactors = FALSE
    )
  )) else data.frame()
  list(
    communication = communication,
    cluster_table = cluster_table,
    direction_table = direction_table,
    h5ad = NULL,
    manifest = list(
      n_obs = ncol(input$expression), n_vars = nrow(input$expression),
      cell_ids = as.list(input$cells), feature_ids = as.list(rownames(input$expression)),
      n_communication_rows = nrow(communication),
      n_cluster_rows = nrow(cluster_table), n_direction_rows = nrow(direction_table),
      coordinate_space = "raw", versions = list(commot = "0.0.3")
    )
  )
}

mock_commot_execute <- function(fail = FALSE, empty = FALSE) {
  testthat::local_mocked_bindings(
    commot_execute = function(input, cluster, direction, ...) {
      if (isTRUE(fail)) stop("COMMOT subprocess failed")
      value <- mock_commot_result(input, cluster = cluster, direction = direction)
      if (isTRUE(empty)) {
        value$communication <- value$communication[0, , drop = FALSE]
        value$manifest$n_communication_rows <- 0L
      }
      value
    },
    .package = "scop",
    .env = parent.frame()
  )
}

test_that("COMMOT environment is isolated and pinned to official Python 3.9", {
  modules <- getFromNamespace("normalize_env_modules", "scop")("commot")
  expect_identical(modules, "commot")
  expect_false("commot" %in% getFromNamespace("default_env_modules", "scop")())
  req <- env_requirements(modules = "commot")
  expect_identical(req$python, "3.9-1")
  expect_match(unname(req$packages[["commot"]]), "zcang/COMMOT.git@d117445bc07e")
  expect_identical(unname(req$packages[["numpy"]]), "numpy==1.22.4")
  expect_error(
    getFromNamespace("normalize_env_modules", "scop")(c("commot", "scanpy")),
    "standalone"
  )
})

test_that("COMMOT runtime preflight checks its pinned standalone stack", {
  python <- tempfile(fileext = ".exe")
  file.create(python)
  checked <- NULL
  testthat::with_mocked_bindings(
    commot_runtime("commot_env", FALSE),
    PrepareEnv = function(...) invisible(NULL),
    check_python = function(packages, ...) {
      checked <<- packages
      TRUE
    },
    resolve_conda = function(...) "conda",
    resolve_python_executable = function(...) python,
    .package = "scop"
  )
  expect_identical(unname(checked[["numpy"]]), "numpy==1.22.4")
  expect_identical(unname(checked[["scanpy"]]), "scanpy==1.8.2")
  expect_match(unname(checked[["commot"]]), "zcang/COMMOT.git@d117445bc07e")
})

test_that("RunCOMMOT stores transport mass without inventing p-values", {
  srt <- make_commot_test_object()
  mock_commot_execute()
  out <- RunCOMMOT(
    srt, group.by = "celltype", coord.cols = c("col", "row"),
    distance.threshold = 2, backend = "r", verbose = FALSE
  )
  expect_identical(out@tools$COMMOT$method, "COMMOT")
  expect_identical(out@tools$COMMOT$parameters$coordinate_space, "raw")
  expect_identical(out@tools$COMMOT$parameters$distance_units, "raw coordinate units")
  expect_identical(out@tools$COMMOT$long_table$score_type, "transport_mass")
  expect_true(is.na(out@tools$COMMOT$long_table$pvalue))
  expect_identical(out@tools$COMMOT$provenance$backend_version, "0.0.3")
  expect_true("COMMOT" %in% out@tools$CCC$methods)
  expect_true(inherits(COMMOTPlot(out, plot_type = "network"), c("ggplot", "recordedplot")))
})

test_that("RunCOMMOT retains official cluster p-values and stored directions", {
  srt <- make_commot_test_object()
  mock_commot_execute()
  out <- RunCOMMOT(
    srt, group.by = "celltype", coord.cols = c("col", "row"),
    cluster = TRUE, direction = TRUE,
    cluster.args = list(pathway_name = "P1", n_permutations = 10),
    direction.args = list(pathway_name = "P1"),
    backend = "r", verbose = FALSE
  )
  expect_equal(out@tools$COMMOT$results$default$cluster_table$pvalue, c(0.01, 0.2))
  expect_equal(nrow(out@tools$COMMOT$results$default$direction_table), 8)
  expect_s3_class(COMMOTPlot(out, plot_type = "matrix"), "ggplot")
  expect_s3_class(COMMOTPlot(out, plot_type = "direction"), "ggplot")
})

test_that("RunCOMMOT never silently selects one image", {
  srt <- make_commot_test_object()
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1), y = c(0, 1), row.names = c("cell1", "cell2")),
    type = "centroids", assay = assay, key = "c1_"
  )
  srt[["slice2"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(2, 3), y = c(1, 2), row.names = c("cell3", "cell4")),
    type = "centroids", assay = assay, key = "c2_"
  )
  expect_error(RunCOMMOT(srt, group.by = "celltype", verbose = FALSE), "Multiple spatial images")
})

test_that("RunCOMMOT accepts one explicitly selected spatial image", {
  srt <- make_commot_test_object()
  assay <- SeuratObject::DefaultAssay(srt)
  srt[["slice1"]] <- SeuratObject::CreateFOV(
    data.frame(x = c(0, 1, 2, 3), y = c(0, 1, 1, 2), row.names = colnames(srt)),
    type = "centroids", assay = assay, key = "c1_"
  )
  mock_commot_execute()
  out <- RunCOMMOT(
    srt, group.by = "celltype", image = "slice1",
    backend = "r", verbose = FALSE
  )
  expect_identical(out@tools$COMMOT$parameters$image, "slice1")
  expect_identical(out@tools$COMMOT$cells, colnames(srt))
})

test_that("COMMOT spatial input rejects invalid IDs and coordinates", {
  srt <- make_commot_test_object()
  duplicate <- data.frame(cell_id = c("cell1", "cell1"), x = c(0, 1), y = c(0, 1))
  expect_error(
    testthat::with_mocked_bindings(
      commot_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
      SpatialCoordinates = function(...) list(data = duplicate, source = list()),
      .package = "scop"
    ),
    "valid unique"
  )
  nonfinite <- data.frame(cell_id = "cell1", x = Inf, y = 0)
  expect_error(
    testthat::with_mocked_bindings(
      commot_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
      SpatialCoordinates = function(...) list(data = nonfinite, source = list()),
      .package = "scop"
    ),
    "finite raw"
  )
  unknown <- data.frame(cell_id = "absent", x = 0, y = 0)
  expect_error(
    testthat::with_mocked_bindings(
      commot_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
      SpatialCoordinates = function(...) list(data = unknown, source = list()),
      .package = "scop"
    ),
    "valid unique"
  )
})

test_that("COMMOT expression, labels, and coordinates use the same ID order", {
  srt <- make_commot_test_object()
  reordered <- data.frame(
    cell_id = rev(colnames(srt)), x = 4:1, y = 1:4,
    stringsAsFactors = FALSE
  )
  input <- testthat::with_mocked_bindings(
    commot_input(srt, "celltype", NULL, "counts", NULL, c("col", "row")),
    SpatialCoordinates = function(...) list(data = reordered, source = list()),
    .package = "scop"
  )
  expect_identical(input$cells, reordered$cell_id)
  expect_identical(colnames(input$expression), reordered$cell_id)
  expect_identical(input$groups, as.character(srt[[]][reordered$cell_id, "celltype"]))
})

test_that("COMMOT failures and empty output do not mutate the input object", {
  srt <- make_commot_test_object()
  before <- list(
    metadata = srt[[]], assays = srt@assays, reductions = srt@reductions,
    graphs = srt@graphs, tools = srt@tools
  )
  mock_commot_execute(fail = TRUE)
  expect_error(
    RunCOMMOT(
      srt, group.by = "celltype", coord.cols = c("col", "row"),
      backend = "r", verbose = FALSE
    ),
    "subprocess failed"
  )
  expect_identical(srt[[]], before$metadata)
  expect_identical(srt@assays, before$assays)
  expect_identical(srt@reductions, before$reductions)
  expect_identical(srt@graphs, before$graphs)
  expect_identical(srt@tools, before$tools)
  mock_commot_execute(empty = TRUE)
  expect_error(
    RunCOMMOT(
      srt, group.by = "celltype", coord.cols = c("col", "row"),
      backend = "r", verbose = FALSE
    ),
    "incomplete"
  )
  expect_identical(srt[[]], before$metadata)
  expect_identical(srt@assays, before$assays)
  expect_identical(srt@reductions, before$reductions)
  expect_identical(srt@graphs, before$graphs)
  expect_identical(srt@tools, before$tools)
})

test_that("RunCCC dispatches COMMOT only when explicitly selected", {
  srt <- make_commot_test_object()
  mock_commot_execute()
  out <- RunCCC(
    srt, group.by = "celltype", methods = "commot",
    method_params = list(COMMOT = list(coord.cols = c("col", "row"))),
    backend = "r", verbose = FALSE
  )
  expect_identical(out@tools$RunCCC$completed_methods, "COMMOT")
  expect_true("COMMOT" %in% out@tools$CCC$methods)
  expect_false("COMMOT" %in% eval(formals(RunCCC)$methods))
})

test_that("COMMOTPlot requires an explicit result when several are stored", {
  srt <- make_commot_test_object()
  mock_commot_execute()
  out <- RunCOMMOT(
    srt, group.by = "celltype", coord.cols = c("col", "row"),
    result.name = "first", backend = "r", verbose = FALSE
  )
  out <- RunCOMMOT(
    out, group.by = "celltype", coord.cols = c("col", "row"),
    result.name = "second", backend = "r", verbose = FALSE
  )
  expect_error(COMMOTPlot(out), "select.*result.name")
  expect_true(inherits(
    COMMOTPlot(out, result.name = "first"),
    c("ggplot", "recordedplot")
  ))
})

test_that("COMMOT H5AD publication is external and overwrite-aware", {
  source <- tempfile(fileext = ".h5ad")
  writeBin(charToRaw("commot"), source)
  destination <- file.path(tempfile("commot-artifact-"), "result.h5ad")
  path <- getFromNamespace("commot_copy_h5ad", "scop")(source, destination, FALSE)
  expect_true(file.exists(path))
  expect_error(
    getFromNamespace("commot_copy_h5ad", "scop")(source, destination, FALSE),
    "already exists"
  )
})

test_that("bundled COMMOT runner uses only official API entry points", {
  path <- getFromNamespace("runner_script_path", "scop")("commot_runner.py", "COMMOT")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(identical(basename(path), "commot.py"))
  expect_match(source, "ct.tl.spatial_communication", fixed = TRUE)
  expect_match(source, "ct.tl.cluster_communication", fixed = TRUE)
  expect_match(source, "ct.tl.communication_direction", fixed = TRUE)
  expect_false(grepl("reticulate", source, fixed = TRUE))
})
