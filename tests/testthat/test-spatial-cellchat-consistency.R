# End-to-end consistency between the scop SpatialCellChat wrapper and the
# original SpatialCellChat pipeline on real Visium data.

spatialcellchat_real_input <- function(n = 120, seed = 42) {
  data(visium_mouse_brain_slices_sub)
  srt <- visium_mouse_brain_slices_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
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

  set.seed(42)
  wrapped <- RunSpatialCellChat(
    srt,
    group.by = "sample",
    technology = "visium",
    species = "Mus_musculus",
    database = "all",
    analysis.level = "spot",
    coordinate.unit = "pixel",
    ratio = 0.5,
    tol = 5,
    min.cells = 5,
    min.links = 1,
    nboot = 2,
    seed.use = 1,
    store.object = "full",
    verbose = FALSE
  )
  chat_w <- wrapped@tools$SpatialCellChat$results$default$ALL$native_object
  expect_s4_class(chat_w, "SpatialCellChat")

  set.seed(42)
  chat_m <- run_original_spatialcellchat(
    srt = srt,
    group.by = "sample",
    species = "Mus_musculus",
    database = "all",
    ratio = 0.5,
    tol = 5
  )

  for (slot in c("prob", "pval", "weight", "count")) {
    expect_equal(chat_w@net[[slot]], chat_m@net[[slot]], label = paste0("net@", slot))
  }
  expect_equal(chat_w@netP$prob, chat_m@netP$prob)
  expect_equal(chat_w@netP$pval, chat_m@netP$pval)
  expect_identical(as.character(chat_w@idents), as.character(chat_m@idents))
})

test_that("RunSpatialCellChat long table matches subsetCommunication extraction", {
  skip_on_cran()
  skip_if_not_installed("SpatialCellChat")

  srt <- spatialcellchat_real_input()

  set.seed(42)
  wrapped <- RunSpatialCellChat(
    srt,
    group.by = "sample",
    technology = "visium",
    species = "Mus_musculus",
    database = "all",
    analysis.level = "spot",
    coordinate.unit = "pixel",
    ratio = 0.5,
    tol = 5,
    min.cells = 5,
    min.links = 1,
    nboot = 2,
    seed.use = 1,
    store.object = "full",
    verbose = FALSE
  )
  chat_w <- wrapped@tools$SpatialCellChat$results$default$ALL$native_object
  long_table <- wrapped@tools$SpatialCellChat$long_table
  expect_true(nrow(long_table) > 0L)

  native <- SpatialCellChat::subsetCommunication(
    chat_w,
    slot.name = "net",
    thresh = 0.05
  )
  native <- as.data.frame(native)
  expect_true(nrow(native) > 0L)

  key <- paste(native$source, native$target, native$interaction_name)
  long_key <- paste(long_table$sender, long_table$receiver, long_table$interaction_name)
  expect_gt(length(intersect(long_key, key)) / length(unique(key)), 0.9)

  merged <- merge(
    data.frame(key = key, prob = native$prob),
    data.frame(key = long_key, score = long_table$score),
    by = "key"
  )
  expect_gt(nrow(merged), 0L)
  expect_equal(merged$score, merged$prob, tolerance = 1e-12)
})
