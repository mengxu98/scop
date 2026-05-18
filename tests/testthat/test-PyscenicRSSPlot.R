make_pyscenic_rss_test_object <- function() {
  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 2, 0, 3, 1,
        0, 1, 0, 2, 0, 3,
        2, 2, 1, 1, 0, 0,
        1, 1, 1, 1, 1, 1
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt@meta.data[["cluster"]] <- factor(
    c("A", "A", "B", "B", "C", "C"),
    levels = c("A", "B", "C")
  )

  scores_cells_by_regulon <- matrix(
    c(
      8, 7, 1, 1, 1, 1,
      1, 1, 8, 7, 1, 1,
      1, 1, 1, 1, 8, 7,
      3, 2, 3, 2, 3, 2
    ),
    nrow = 6,
    byrow = FALSE
  )
  rownames(scores_cells_by_regulon) <- colnames(srt)
  colnames(scores_cells_by_regulon) <- c("Areg(+)", "Breg(+)", "Creg(+)", "Shared(+)")
  srt@tools[["Pyscenic"]] <- list(scores_cells_by_regulon = scores_cells_by_regulon)

  srt
}

test_that("PyscenicRSSPlot computes RSS from tools slot", {
  srt <- make_pyscenic_rss_test_object()

  out <- PyscenicRSSPlot(
    srt = srt,
    group.by = "cluster",
    top_n = 2,
    combine = FALSE,
    verbose = FALSE
  )

  expect_equal(dim(out$rss_matrix), c(4, 3))
  expect_equal(colnames(out$rss_matrix), c("A", "B", "C"))
  expect_equal(length(out$plots), 3)
  expect_true(all(c("group", "regulon", "TF", "specificity_score", "rank", "is_top", "is_highlight") %in% colnames(out$rank_table)))
  expect_equal(as.integer(table(out$top_table$group)), c(2L, 2L, 2L))
  expect_true(inherits(out$plots[[1]], "ggplot"))
})

test_that("PyscenicRSSPlot rank table is sorted within each group", {
  srt <- make_pyscenic_rss_test_object()

  out <- PyscenicRSSPlot(
    srt = srt,
    group.by = "cluster",
    top_n = 2,
    combine = FALSE,
    verbose = FALSE
  )

  rank_by_group <- split(out$rank_table, out$rank_table$group)
  expect_true(all(vapply(rank_by_group, function(one_group) {
    all(diff(one_group$specificity_score) <= 0)
  }, logical(1))))
  expect_equal(out$top_table$rank[out$top_table$group == "A"], 1:2)
})

test_that("PyscenicRSSPlot falls back to pyscenic assay", {
  srt <- make_pyscenic_rss_test_object()
  auc_mat <- t(srt@tools[["Pyscenic"]][["scores_cells_by_regulon"]])
  srt@tools[["Pyscenic"]] <- NULL
  suppressWarnings({
    srt[["pyscenic"]] <- Seurat::CreateAssayObject(data = auc_mat)
  })

  out <- PyscenicRSSPlot(
    srt = srt,
    group.by = "cluster",
    top_n = 1,
    combine = FALSE,
    verbose = FALSE
  )

  expect_equal(dim(out$rss_matrix), c(4, 3))
  expect_equal(as.integer(table(out$top_table$group)), c(1L, 1L, 1L))
})

test_that("PyscenicRSSPlot highlights TF and regulon names", {
  srt <- make_pyscenic_rss_test_object()

  out_tf <- PyscenicRSSPlot(
    srt = srt,
    group.by = "cluster",
    top_n = 1,
    highlight_tf = "Breg",
    combine = FALSE,
    verbose = FALSE
  )

  expect_true(any(out_tf$rank_table$is_highlight))
  expect_true(all(out_tf$rank_table$regulon[out_tf$rank_table$is_highlight] == "Breg(+)"))
  expect_equal(length(out_tf$plots), 3)

  out_regulon <- PyscenicRSSPlot(
    srt = srt,
    group.by = "cluster",
    top_n = 1,
    highlight_tf = "Breg(+)",
    combine = FALSE,
    verbose = FALSE
  )

  expect_true(any(out_regulon$rank_table$is_highlight))
  expect_true(all(out_regulon$rank_table$TF[out_regulon$rank_table$is_highlight] == "Breg"))
})

test_that("PyscenicRSSPlot keeps plots when highlight TF is missing", {
  srt <- make_pyscenic_rss_test_object()

  out <- PyscenicRSSPlot(
    srt = srt,
    group.by = "cluster",
    top_n = 1,
    highlight_tf = "MissingTF",
    combine = FALSE,
    verbose = FALSE
  )

  expect_false(any(out$rank_table$is_highlight))
  expect_equal(length(out$plots), 3)
  expect_true(inherits(out$plots[[1]], "ggplot"))
})

test_that("PyscenicRSSPlot reports missing metadata and cell mismatches", {
  srt <- make_pyscenic_rss_test_object()

  expect_error(
    PyscenicRSSPlot(srt, group.by = "missing_cluster", verbose = FALSE),
    "group.by"
  )

  rownames(srt@tools[["Pyscenic"]][["scores_cells_by_regulon"]]) <- paste0(
    "other_cell",
    seq_len(ncol(srt))
  )
  expect_error(
    PyscenicRSSPlot(srt, group.by = "cluster", verbose = FALSE),
    "No shared cells"
  )
})

test_that("PyscenicRSSPlot reports all missing group annotations", {
  srt <- make_pyscenic_rss_test_object()
  srt@meta.data[["all_na"]] <- NA_character_

  expect_error(
    PyscenicRSSPlot(srt, group.by = "all_na", verbose = FALSE),
    "All cells have missing values"
  )
})
