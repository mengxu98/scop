test_that("CytoSPACE singleton sampling preserves reference cell identities", {
  ref_expr <- matrix(
    seq_len(12),
    nrow = 4,
    dimnames = list(paste0("gene", 1:4), c("other1", "other2", "singleton"))
  )
  labels <- factor(c("other", "other", "singleton"))

  for (desired in c(1L, 4L)) {
    for (seed in 1:5) {
      sampled <- cytospace_sample_reference_cells(
        ref_expr = ref_expr,
        labels = labels,
        cell_types = c("other", "singleton"),
        target_cell_type_counts = c(other = 0L, singleton = desired),
        scRNA_max_transcripts_per_cell = NULL,
        seed = seed
      )

      expect_identical(sampled$source_index, rep(3L, desired))
      expect_identical(sampled$cell_ids, rep("singleton", desired))
      expect_identical(sampled$labels, rep("singleton", desired))
      expect_equal(
        unname(sampled$expr),
        unname(ref_expr[, rep(3L, desired), drop = FALSE])
      )
      expect_identical(anyDuplicated(colnames(sampled$expr)), 0L)
    }
  }
})

test_that("CytoSPACE multi-cell sampling retains type counts and seed replay", {
  ref_expr <- matrix(
    seq_len(24),
    nrow = 4,
    dimnames = list(paste0("gene", 1:4), paste0("reference", 1:6))
  )
  labels <- c("B", "A", "B", "A", "C", "A")
  sample_cells <- function() {
    cytospace_sample_reference_cells(
      ref_expr = ref_expr,
      labels = labels,
      cell_types = c("A", "B", "C"),
      target_cell_type_counts = c(A = 2L, B = 5L, C = 0L),
      scRNA_max_transcripts_per_cell = NULL,
      seed = 23L
    )
  }
  sampled <- sample_cells()

  expect_identical(sampled, sample_cells())
  expect_identical(sampled$labels, c(rep("A", 2L), rep("B", 5L)))
  expect_identical(sampled$cell_ids, colnames(ref_expr)[sampled$source_index])
  expect_identical(sampled$labels, labels[sampled$source_index])
  expect_length(unique(sampled$source_index[1:2]), 2L)
  expect_true(all(c(1L, 3L) %in% sampled$source_index[3:7]))
})

test_that("CytoSPACE summaries preserve reserved and sanitized type labels", {
  cell_types <- c("A-B", "A.B", "A.B.1", "Total cells", "Total cells.1")
  assignments <- data.frame(
    SpotID = c("spot1", "spot1", "spot1", "spot2", "spot2", "spot2"),
    CellType = c("A-B", "A-B", "A.B", "A.B.1", "Total cells", "Total cells.1")
  )
  summary <- cytospace_build_spot_summary(
    assignments,
    c("spot1", "spot2", "empty"),
    cell_types
  )

  expect_identical(colnames(summary$fractions), cell_types)
  expect_identical(summary$total_cells_column, "Total cells.2")
  expect_identical(colnames(summary$counts), c(cell_types, "Total cells.2"))
  expect_equal(summary$counts[["Total cells"]], c(0L, 1L, 0L))
  expect_equal(summary$counts[["Total cells.1"]], c(0L, 1L, 0L))
  expect_equal(summary$counts[[summary$total_cells_column]], c(3L, 3L, 0L))
  expect_equal(summary$fractions[["Total cells"]], c(0, 1 / 3, 0))
  expect_equal(unname(rowSums(summary$fractions)), c(1, 1, 0))
  expect_identical(names(summary$metadata_cell_type_suffixes), cell_types)
  expect_identical(
    unname(summary$metadata_cell_type_suffixes),
    c("A.B.2", "A.B", "A.B.1", "Total.cells", "Total.cells.1")
  )
})

test_that("CytoSPACE empty assignments preserve spots and expression dimensions", {
  spot_ids <- c("spot1", "spot2")
  coords <- data.frame(x = c(0, 1), y = c(1, 0), row.names = spot_ids)
  assignments <- cytospace_build_assignment_table(
    result = list(cell_index = integer(), spot_index = integer(), score = numeric()),
    sampled_cells = character(),
    sampled_labels = character(),
    spot_ids = spot_ids,
    coords = coords
  )
  summary <- cytospace_build_spot_summary(assignments, spot_ids, c("A", "B"))
  expr <- cytospace_build_assigned_expression(
    matrix(numeric(), nrow = 3, ncol = 0, dimnames = list(paste0("gene", 1:3), NULL)),
    assignments
  )

  expect_identical(dim(assignments), c(0L, 7L))
  expect_identical(assignments$UniqueCID, character())
  expect_identical(rownames(summary$counts), spot_ids)
  expect_equal(unname(as.matrix(summary$counts)), matrix(0, 2, 3))
  expect_equal(unname(as.matrix(summary$fractions)), matrix(0, 2, 2))
  expect_identical(summary$total_cells_column, "Total cells")
  expect_identical(dim(expr), c(3L, 0L))
})

test_that("CytoSPACE metadata preserves colliding types and leaves empty spots unlabeled", {
  counts <- matrix(
    seq_len(9),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), c("empty", "spot2", "spot1"))
  )
  spatial <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  assignments <- data.frame(
    SpotID = c("spot1", "spot1", "spot1", "spot2"),
    CellType = c("A-B", "A-B", "A.B", "Total cells")
  )
  summary <- cytospace_build_spot_summary(
    assignments,
    c("spot1", "spot2", "empty"),
    c("A-B", "A.B", "Total cells")
  )
  out <- cytospace_add_metadata(spatial, summary, prefix = "test")
  metadata <- out[[]]

  expect_identical(rownames(metadata), colnames(spatial))
  expect_equal(metadata$test_total_cells, c(0L, 1L, 3L))
  expect_equal(metadata$test_count_A.B, c(0L, 0L, 1L))
  expect_equal(metadata$test_count_A.B.1, c(0L, 0L, 2L))
  expect_equal(metadata$test_count_Total.cells, c(0L, 1L, 0L))
  expect_equal(metadata$test_frac_A.B, c(0, 0, 1 / 3))
  expect_equal(metadata$test_frac_A.B.1, c(0, 0, 2 / 3))
  expect_equal(metadata$test_frac_Total.cells, c(0, 1, 0))
  expect_identical(metadata$test_dominant_type, c(NA_character_, "Total cells", "A-B"))
})

make_cytospace_alignment_objects <- function() {
  genes <- paste0("gene", 1:5)
  counts <- matrix(
    c(12, 1, 3, 2, 4, 1, 13, 2, 3, 5, 3, 2, 14, 1, 6, 2, 4, 1, 15, 3),
    nrow = 5,
    dimnames = list(genes, paste0("ref", 1:4))
  )
  reference <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  reference$cell_type <- c("A-B", "A.B", "Total cells", "Total cells.1")
  spatial_counts <- counts[, 1:3]
  colnames(spatial_counts) <- c("empty", "spot2", "spot1")
  spatial <- suppressWarnings(SeuratObject::CreateSeuratObject(spatial_counts))
  spatial$col <- c(0, 1, 2)
  spatial$row <- c(1, 0, 1)
  list(spatial = spatial, reference = reference)
}

test_that("CytoSPACE public workflow preserves cell types and zero-capacity spots", {
  objects <- make_cytospace_alignment_objects()
  cell_types <- objects$reference$cell_type
  out <- RunCytoSPACE(
    objects$spatial,
    reference = objects$reference,
    reference_label = "cell_type",
    features = rownames(objects$spatial),
    cell_fractions = stats::setNames(rep(0.25, 4), cell_types),
    n_cells_per_spot = c(spot1 = 2L, empty = 0L, spot2 = 2L),
    scRNA_max_transcripts_per_cell = NULL,
    seed = 1L,
    prefix = "test",
    verbose = FALSE
  )
  results <- out@tools$CytoSPACE
  assignments <- results$assigned_locations
  metadata <- out[[]]
  original_labels <- stats::setNames(cell_types, colnames(objects$reference))

  expect_setequal(assignments$OriginalCID, colnames(objects$reference))
  expect_identical(assignments$CellType, unname(original_labels[assignments$OriginalCID]))
  expect_identical(as.integer(table(factor(assignments$CellType, levels = cell_types))), rep(1L, 4))
  expect_identical(rownames(results$cell_type_assignments_by_spot), colnames(out))
  expect_identical(results$total_cells_column, "Total cells.2")
  expect_identical(results$metadata_cell_type_suffixes[["A-B"]], "A.B.1")
  expect_identical(results$metadata_cell_type_suffixes[["A.B"]], "A.B")
  expect_equal(metadata$test_total_cells, c(0L, 2L, 2L))
  expect_true(is.na(metadata["empty", "test_dominant_type"]))
  expect_false(anyNA(metadata[c("spot2", "spot1"), "test_dominant_type"]))
  for (ct in cell_types) {
    suffix <- results$metadata_cell_type_suffixes[[ct]]
    expect_equal(
      metadata[[paste0("test_count_", suffix)]],
      results$cell_type_assignments_by_spot[[ct]]
    )
    expect_equal(
      metadata[[paste0("test_frac_", suffix)]],
      results$fractional_abundances_by_spot[[ct]]
    )
    expect_equal(sum(metadata[[paste0("test_count_", suffix)]]), 1)
  }
  expect_identical(colnames(results$assigned_expression), assignments$UniqueCID)
  expect_equal(
    unname(results$assigned_expression),
    unname(as.matrix(GetAssayData5(objects$reference, layer = "counts"))[
      , assignments$OriginalCID, drop = FALSE
    ])
  )
})

test_that("CytoSPACE public workflow supports an entirely empty allocation", {
  objects <- make_cytospace_alignment_objects()
  out <- RunCytoSPACE(
    objects$spatial,
    reference = objects$reference,
    reference_label = "cell_type",
    features = rownames(objects$spatial),
    cell_fractions = stats::setNames(rep(0.25, 4), objects$reference$cell_type),
    n_cells_per_spot = rep(0L, ncol(objects$spatial)),
    scRNA_max_transcripts_per_cell = NULL,
    verbose = FALSE
  )

  expect_identical(nrow(out@tools$CytoSPACE$assigned_locations), 0L)
  expect_identical(dim(out@tools$CytoSPACE$assigned_expression), c(5L, 0L))
  expect_equal(
    out$CytoSPACE_total_cells,
    stats::setNames(rep(0L, ncol(out)), colnames(out))
  )
  expect_true(all(is.na(out$CytoSPACE_dominant_type)))
  expect_true(all(out@tools$CytoSPACE$fractional_abundances_by_spot == 0))
})
