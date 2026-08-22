# End-to-end consistency between the scop SpatialCellChat wrapper and the
# original SpatialCellChat pipeline on real Visium data.

spatialcellchat_real_input <- function(n = 80, seed = 42) {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  srt$sample <- ifelse(srt$y > stats::median(srt$y), "slice_a", "slice_b")
  set.seed(seed)
  srt <- srt[, sample(colnames(srt), n)]
  srt
}

run_original_spatialcellchat <- function(
  srt,
  group.by,
  species,
  database,
  analysis.level = "spot",
  composition = NULL,
  ratio,
  tol,
  interaction.range = 250,
  scale.distance = 0.2,
  min.cells = 5,
  min.links = 1,
  avg.type = "avg",
  do.permutation = TRUE,
  nboot = 2,
  seed.use = 1
) {
  cells <- colnames(srt)
  metric <- getFromNamespace("spatialcellchat_metric_coordinates", "scop")(
    srt = srt,
    cells = cells,
    image = NULL,
    coord.cols = c("x", "y"),
    technology = "visium",
    coordinate.unit = "pixel",
    ratio = ratio,
    tol = tol
  )
  expression <- Seurat::GetAssayData(srt, assay = SeuratObject::DefaultAssay(srt), layer = "data")
  metadata <- data.frame(
    labels = as.character(srt@meta.data[cells, group.by]),
    samples = "ALL",
    row.names = cells,
    stringsAsFactors = FALSE
  )
  database.use <- getFromNamespace("spatialcellchat_database", "scop")(
    species = species,
    database = database
  )
  scc_call <- getFromNamespace("spatialcellchat_call", "scop")
  set.seed(as.integer(seed.use))
  chat <- scc_call("createSpatialCellChat", list(
    object = expression,
    meta = metadata,
    group.by = "labels",
    datatype = "spatial",
    coordinates = as.matrix(metric$data[, c("x", "y"), drop = FALSE]),
    spatial.factors = metric$spatial.factors
  ), analysis.level)
  chat@DB <- database.use
  chat <- scc_call("subsetData", list(object = chat), analysis.level)
  chat <- scc_call("preProcessing", list(object = chat), analysis.level)
  chat <- scc_call("identifyOverExpressedGenes", list(
    object = chat,
    selection.method = "meringue",
    do.grid = FALSE
  ), analysis.level)
  chat <- scc_call("identifyOverExpressedInteractions", list(
    object = chat,
    variable.both = FALSE
  ), analysis.level)
  probability_args <- list(
    object = chat,
    raw.use = TRUE,
    distance.use = TRUE,
    scale.distance = scale.distance,
    interaction.range = interaction.range,
    contact.dependent = FALSE
  )
  chat <- scc_call("computeCommunProb", probability_args, analysis.level)
  chat <- scc_call("filterProbability", list(
    object = chat,
    nboot = nboot,
    seed.use = seed.use,
    thresh = 0.05
  ), analysis.level)
  chat <- scc_call("filterCommunication", list(
    object = chat,
    min.cells = NULL,
    min.links = min.links,
    min.cells.sr = min.cells
  ), analysis.level)
  if (identical(analysis.level, "composition")) {
    chat <- scc_call("computeAvgCommunProb_Visium", list(
      object = chat,
      cell.type.decomposition = composition,
      avg.type = avg.type,
      do.permutation = do.permutation,
      nboot = nboot,
      seed.use = seed.use
    ), analysis.level)
  } else {
    chat <- scc_call("computeAvgCommunProb", list(
      object = chat,
      group.by = "labels",
      avg.type = avg.type,
      min.cells.sr = min.cells,
      do.permutation = do.permutation,
      nboot = nboot,
      seed.use = seed.use
    ), analysis.level)
  }
  chat <- scc_call("filterCommunication", list(
    object = chat,
    min.cells = min.cells,
    min.links = NULL,
    min.cells.sr = NULL
  ), analysis.level)
  chat <- scc_call("computeCommunProbPathway", list(object = chat), analysis.level)
  scc_call("aggregateNet", list(object = chat), analysis.level)
}

test_that("RunSpatialCellChat matches the original SpatialCellChat pipeline", {
  skip_on_cran()
  skip_if_not_installed("SpatialCellChat")

  srt <- spatialcellchat_real_input()

  # Group-level permutation tests are disabled to keep the suite fast; the
  # wrapped and original runs still share identical inputs and settings, so
  # net@prob/@weight/@count and netP$prob remain strong consistency checks.
  set.seed(42)
  wrapped <- RunSpatialCellChat(
    srt,
    group.by = "sample",
    technology = "visium",
    species = "Homo_sapiens",
    database = "protein",
    analysis.level = "spot",
    coordinate.unit = "pixel",
    ratio = 0.5,
    tol = 5,
    min.cells = 5,
    min.links = 1,
    nboot = 2,
    seed.use = 1,
    do.permutation = FALSE,
    store.object = "full",
    verbose = FALSE
  )
  chat_w <- wrapped@tools$SpatialCellChat$results$default$ALL$native_object
  expect_s4_class(chat_w, "SpatialCellChat")

  set.seed(42)
  chat_m <- run_original_spatialcellchat(
    srt = srt,
    group.by = "sample",
    species = "Homo_sapiens",
    database = "protein",
    ratio = 0.5,
    tol = 5,
    do.permutation = FALSE
  )

  for (slot in c("prob", "pval", "weight", "count")) {
    expect_equal(chat_w@net[[slot]], chat_m@net[[slot]], label = paste0("net@", slot))
  }
  expect_equal(chat_w@netP$prob, chat_m@netP$prob)
  expect_equal(chat_w@netP$pval, chat_m@netP$pval)
  expect_identical(as.character(chat_w@idents), as.character(chat_m@idents))

  # With permutation disabled, upstream subsetCommunication is unavailable for
  # ligand-receptor records, so the bundled long table falls back to
  # pathway-level records; assert it is built and structurally consistent.
  long_table <- wrapped@tools$SpatialCellChat$long_table
  expect_true(nrow(long_table) > 0L)
  expect_identical(unique(long_table$result_level), "pathway")
  expect_true(all(c("sender", "receiver", "pathway_name", "score") %in% colnames(long_table)))
})
