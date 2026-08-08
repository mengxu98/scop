# End-to-end consistency between the scop wrapper and the original backend
# method, run with real data. These tests do not mock the backend: they run
# `RunCellChat` and the CellChat package's own vignette pipeline on the same
# input and require the stored results to be identical.

cellchat_consistency_input <- function(n_cells = 120) {
  data(pancreas_sub)
  set.seed(42)
  cells <- sample(colnames(pancreas_sub), n_cells)
  srt <- pancreas_sub[, cells]
  srt$CellType <- factor(srt$CellType)
  Seurat::NormalizeData(srt, verbose = FALSE)
}

run_original_cellchat <- function(
  srt,
  species = "Mus_musculus",
  thresh = 0.05,
  min.cells = 5,
  do.fast = TRUE
) {
  data.input <- Seurat::GetAssayData(srt, assay = "RNA", layer = "data")
  meta <- data.frame(label = srt$CellType, row.names = colnames(srt))
  object <- CellChat::createCellChat(
    object = data.input,
    meta = meta,
    group.by = "label"
  )
  object@DB <- CellChat::CellChatDB.mouse
  object <- CellChat::subsetData(object)
  object <- CellChat::identifyOverExpressedGenes(
    object,
    thresh.p = thresh,
    min.cells = min.cells,
    do.fast = do.fast
  )
  object <- CellChat::identifyOverExpressedInteractions(object)
  object <- CellChat::computeCommunProb(object)
  object <- CellChat::filterCommunication(object, min.cells = min.cells)
  object <- CellChat::computeCommunProbPathway(object, thresh = thresh)
  object <- CellChat::aggregateNet(object, thresh = thresh)
  object <- CellChat::netAnalysis_computeCentrality(object, thresh = thresh)
  object
}

test_that("RunCellChat stores a CellChat object identical to the original pipeline", {
  skip_if_not_installed("CellChat")
  skip_if_not_installed("presto")

  srt <- cellchat_consistency_input()

  set.seed(42)
  wrapped <- RunCellChat(
    srt,
    group.by = "CellType",
    species = "Mus_musculus",
    thresh = 0.05,
    min.cells = 5,
    do.fast = TRUE,
    backend = "cpp",
    verbose = FALSE
  )
  cc_wrap <- wrapped@tools$CellChat$results$ALL$cellchat_object

  set.seed(42)
  cc_orig <- run_original_cellchat(srt)

  for (slot in c("prob", "pval", "weight", "count")) {
    expect_equal(cc_wrap@net[[slot]], cc_orig@net[[slot]], label = paste0("net@", slot))
  }
  expect_equal(cc_wrap@netP$prob, cc_orig@netP$prob)
  expect_equal(cc_wrap@netP$pval, cc_orig@netP$pval)
  expect_equal(cc_wrap@data.signaling, cc_orig@data.signaling)
  expect_identical(as.character(cc_wrap@idents), as.character(cc_orig@idents))
})

test_that("RunCellChat long table matches CellChat::subsetCommunication", {
  skip_if_not_installed("CellChat")
  skip_if_not_installed("presto")

  srt <- cellchat_consistency_input()

  set.seed(42)
  wrapped <- RunCellChat(
    srt,
    group.by = "CellType",
    species = "Mus_musculus",
    thresh = 0.05,
    min.cells = 5,
    do.fast = TRUE,
    backend = "cpp",
    verbose = FALSE
  )
  cc_wrap <- wrapped@tools$CellChat$results$ALL$cellchat_object

  native <- CellChat::subsetCommunication(
    cc_wrap,
    slot.name = "net",
    thresh = 0.05
  )
  native <- as.data.frame(native)
  long <- wrapped@tools$CellChat$long_table
  long <- long[!is.na(long$sender) & !is.na(long$receiver), , drop = FALSE]

  expect_equal(nrow(long), nrow(native))
  key <- paste(native$source, native$target, native$interaction_name)
  long_key <- paste(long$sender, long$receiver, long$interaction_name)
  expect_setequal(long_key, key)

  merged <- merge(
    data.frame(key = key, prob = native$prob, pval = native$pval),
    data.frame(key = long_key, score = long$score, pvalue = long$pvalue),
    by = "key"
  )
  expect_equal(merged$score, merged$prob, tolerance = 1e-12)
  expect_equal(merged$pvalue, merged$pval, tolerance = 1e-12)
})

test_that("RunCellChat cpp and r backends produce identical pair tables", {
  skip_if_not_installed("CellChat")
  skip_if_not_installed("presto")

  srt <- cellchat_consistency_input()

  set.seed(42)
  cpp_out <- RunCellChat(
    srt,
    group.by = "CellType",
    species = "Mus_musculus",
    thresh = 0.05,
    min.cells = 5,
    do.fast = TRUE,
    backend = "cpp",
    verbose = FALSE
  )
  set.seed(42)
  r_out <- RunCellChat(
    srt,
    group.by = "CellType",
    species = "Mus_musculus",
    thresh = 0.05,
    min.cells = 5,
    do.fast = TRUE,
    backend = "r",
    verbose = FALSE
  )

  pt_cpp <- cpp_out@tools$CellChat$pair_table
  pt_r <- r_out@tools$CellChat$pair_table
  pt_cpp <- pt_cpp[order(pt_cpp$sender, pt_cpp$receiver), , drop = FALSE]
  pt_r <- pt_r[order(pt_r$sender, pt_r$receiver), , drop = FALSE]
  rownames(pt_cpp) <- NULL
  rownames(pt_r) <- NULL
  expect_equal(pt_cpp, pt_r, tolerance = 1e-12)
})

test_that("RunCellChat group comparison mode matches manual mergeCellChat", {
  skip_if_not_installed("CellChat")
  skip_if_not_installed("presto")

  srt <- cellchat_consistency_input(n_cells = 200)
  srt$batch <- rep(c("g1", "g2"), length.out = ncol(srt))

  set.seed(42)
  wrapped <- RunCellChat(
    srt,
    group.by = "CellType",
    species = "Mus_musculus",
    group_column = "batch",
    thresh = 0.05,
    min.cells = 5,
    do.fast = TRUE,
    backend = "cpp",
    verbose = FALSE
  )

  comparisons <- wrapped@tools$CellChat$comparisons
  expect_true("g1_vs_g2" %in% names(comparisons))
  merged <- comparisons[["g1_vs_g2"]]$comparison_object
  expect_s4_class(merged, "CellChat")

  objects <- comparisons[["g1_vs_g2"]]$object.list
  expect_identical(names(objects), c("g1", "g2"))

  manual <- CellChat::mergeCellChat(
    object.list = objects,
    add.names = names(objects)
  )
  expect_equal(merged@meta, manual@meta)
  expect_equal(merged@idents, manual@idents)
})
