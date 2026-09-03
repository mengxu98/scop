make_unbalanced_pseudobulk_srt <- function(seed = 20260903) {
  set.seed(seed)
  n_bg <- 80L
  genes <- c("marker", paste0("bg", seq_len(n_bg)))

  spec <- list(
    list(sample = "CTRL_1", group = "CTRL", type = "TypeA", n = 40L, marker = 20),
    list(sample = "CTRL_2", group = "CTRL", type = "TypeA", n = 40L, marker = 20),
    list(sample = "CTRL_3", group = "CTRL", type = "TypeA", n = 40L, marker = 20),
    list(sample = "Rp_1", group = "Rp", type = "TypeA", n = 40L, marker = 1),
    list(sample = "Rp_2", group = "Rp", type = "TypeA", n = 40L, marker = 1),
    list(sample = "Rp_tiny", group = "Rp", type = "TypeA", n = 2L, marker = 400),
    list(sample = "CTRL_1", group = "CTRL", type = "TypeB", n = 20L, marker = 20),
    list(sample = "CTRL_2", group = "CTRL", type = "TypeB", n = 20L, marker = 20),
    list(sample = "CTRL_3", group = "CTRL", type = "TypeB", n = 20L, marker = 20),
    list(sample = "Rp_1", group = "Rp", type = "TypeB", n = 20L, marker = 1),
    list(sample = "Rp_2", group = "Rp", type = "TypeB", n = 20L, marker = 1),
    list(sample = "Rp_tiny", group = "Rp", type = "TypeB", n = 20L, marker = 1)
  )

  blocks <- lapply(spec, function(s) {
    n_cells <- s$n
    mat <- matrix(
      stats::rpois(length(genes) * n_cells, lambda = 8),
      nrow = length(genes),
      ncol = n_cells
    )
    rownames(mat) <- genes
    mat["marker", ] <- stats::rpois(n_cells, lambda = s$marker)
    colnames(mat) <- paste(s$sample, s$type, seq_len(n_cells), sep = "_")
    list(
      counts = mat,
      meta = data.frame(
        orig.ident = rep(s$sample, n_cells),
        group = rep(s$group, n_cells),
        subcelltype = rep(s$type, n_cells),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
      )
    )
  })

  counts <- do.call(cbind, lapply(blocks, `[[`, "counts"))
  meta <- do.call(rbind, lapply(blocks, `[[`, "meta"))
  srt <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  srt$orig.ident <- meta[colnames(srt), "orig.ident"]
  srt$group <- meta[colnames(srt), "group"]
  srt$subcelltype <- meta[colnames(srt), "subcelltype"]
  Seurat::NormalizeData(srt, verbose = FALSE)
}

marker_logfc <- function(markers, group = NULL) {
  rows <- markers[markers$gene == "marker", , drop = FALSE]
  if (!is.null(group)) {
    rows <- rows[as.character(rows$group1) == group, , drop = FALSE]
  }
  as.numeric(rows$avg_log2FC[[1L]])
}

run_sample_de <- function(srt, min.cells.sample = NULL) {
  args <- list(
    object = srt,
    assay = "RNA",
    layer = "counts",
    group.by = "subcelltype",
    group1 = "CTRL",
    group2 = "Rp",
    sample_col = "orig.ident",
    condition_col = "group",
    test.use = "edgeR",
    only.pos = FALSE,
    fc.threshold = 1,
    p.adjust.method = "BH",
    cores = 1,
    verbose = FALSE
  )
  if (!is.null(min.cells.sample)) {
    args$min.cells.sample <- min.cells.sample
  }
  suppressWarnings(do.call(RunDEtest, args))
}

test_that("a low-cell sample can reverse sample-level DE direction versus cell-level tests", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")

  srt <- make_unbalanced_pseudobulk_srt()
  srt_a <- srt[, srt$subcelltype == "TypeA"]
  cell_out <- suppressWarnings(RunDEtest(
    srt_a,
    group.by = "group",
    group1 = "CTRL",
    group2 = "Rp",
    test.use = "wilcox",
    only.pos = FALSE,
    fc.threshold = 1,
    min.pct = 0,
    verbose = FALSE
  ))
  cell_fc <- -marker_logfc(cell_out@tools$DEtest_custom$AllMarkers_wilcox)
  sample_out <- run_sample_de(srt, min.cells.sample = 1)
  sample_markers <- sample_out@tools$DEtest_subcelltype$AllMarkers_edgeR
  typea_fc <- marker_logfc(sample_markers, "TypeA")
  typeb_fc <- marker_logfc(sample_markers, "TypeB")

  expect_lt(cell_fc, 0)
  expect_gt(typea_fc, 0)
  expect_lt(typeb_fc, 0)
  expect_true(cell_fc * typea_fc < 0)
})

test_that("min.cells.sample drops low-cell samples and restores DE direction", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")

  srt <- make_unbalanced_pseudobulk_srt()
  srt_a <- srt[, srt$subcelltype == "TypeA"]
  cell_out <- suppressWarnings(RunDEtest(
    srt_a,
    group.by = "group",
    group1 = "CTRL",
    group2 = "Rp",
    test.use = "wilcox",
    only.pos = FALSE,
    fc.threshold = 1,
    min.pct = 0,
    verbose = FALSE
  ))
  cell_fc <- -marker_logfc(cell_out@tools$DEtest_custom$AllMarkers_wilcox)
  sample_out <- run_sample_de(srt)
  tools_out <- sample_out@tools$DEtest_subcelltype
  sample_markers <- tools_out$AllMarkers_edgeR
  typea <- sample_markers[
    sample_markers$gene == "marker" & as.character(sample_markers$group1) == "TypeA",
    ,
    drop = FALSE
  ]
  typeb <- sample_markers[
    sample_markers$gene == "marker" & as.character(sample_markers$group1) == "TypeB",
    ,
    drop = FALSE
  ]
  inventory <- tools_out$sample_inventory
  typea_tiny <- inventory[
    inventory$group == "TypeA" & inventory$sample == "Rp_tiny",
    ,
    drop = FALSE
  ]
  typeb_tiny <- inventory[
    inventory$group == "TypeB" & inventory$sample == "Rp_tiny",
    ,
    drop = FALSE
  ]

  expect_lt(cell_fc, 0)
  expect_lt(marker_logfc(sample_markers, "TypeA"), 0)
  expect_true(cell_fc * marker_logfc(sample_markers, "TypeA") > 0)
  expect_lt(marker_logfc(sample_markers, "TypeB"), 0)
  expect_equal(typea$sample_number1[[1L]], 3L)
  expect_equal(typea$sample_number2[[1L]], 2L)
  expect_equal(typeb$sample_number1[[1L]], 3L)
  expect_equal(typeb$sample_number2[[1L]], 3L)
  expect_equal(tools_out$min.cells.sample, 10)
  expect_equal(typea_tiny$n_cells[[1L]], 2L)
  expect_false(typea_tiny$kept[[1L]])
  expect_equal(typeb_tiny$n_cells[[1L]], 20L)
  expect_true(typeb_tiny$kept[[1L]])
})

test_that("min.cells.sample skips a group that loses a condition", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")

  srt <- make_unbalanced_pseudobulk_srt()
  sample_out <- run_sample_de(srt, min.cells.sample = 30)
  markers <- sample_out@tools$DEtest_subcelltype$AllMarkers_edgeR
  inventory <- sample_out@tools$DEtest_subcelltype$sample_inventory

  expect_gt(nrow(markers), 0)
  expect_true(all(as.character(markers$group1) == "TypeA"))
  expect_true(any(inventory$group == "TypeB"))
  expect_false(any(inventory$kept[inventory$group == "TypeB"]))
  expect_true(any(inventory$kept[inventory$group == "TypeA"]))
  expect_lt(marker_logfc(markers, "TypeA"), 0)
})

test_that("min.cells.sample must be a single number >= 1", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")

  srt <- make_unbalanced_pseudobulk_srt()
  expect_error(
    run_sample_de(srt, min.cells.sample = 0),
    "min.cells.sample"
  )
  expect_error(
    suppressWarnings(RunDEtest(
      srt,
      group.by = "subcelltype",
      test.use = "wilcox",
      min.cells.sample = -1,
      verbose = FALSE
    )),
    "min.cells.sample"
  )
})
