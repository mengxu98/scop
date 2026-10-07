# These tests capture wrapper inputs; they do not replace real-backend validation.
make_deconv_alignment_pair <- function(reverse_order = FALSE) {
  st <- matrix(c(10, 20, 970, 20, 10, 1970, 15, 15, 2970), 3,
    dimnames = list(c("g1", "g2", "other"), c("s1", "s2", "s3")))
  ref <- matrix(c(9, 1, 990, 8, 2, 990, 1, 9, 1990, 2, 8, 1990, 5, 5, 90), 3,
    dimnames = list(rownames(st), c("A1", "A2", "B1", "B2", "C1")))
  st <- Matrix::Matrix(st, sparse = TRUE)
  ref <- Matrix::Matrix(ref, sparse = TRUE)
  if (isTRUE(reverse_order)) {
    st <- st[, rev(colnames(st))]
    ref <- ref[, rev(colnames(ref))]
  }
  spatial <- SeuratObject::CreateSeuratObject(st)
  spatial$x <- as.integer(sub("s", "", colnames(spatial)))
  spatial$y <- c(s1 = 1, s2 = 2, s3 = 1)[colnames(spatial)]
  spatial[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = st[, c("s1", "s3")])
  reference <- SeuratObject::CreateSeuratObject(ref)
  reference$celltype <- c(A1 = "A", A2 = "A", B1 = "B", B2 = "B", C1 = NA)[colnames(reference)]
  reference$sample <- c(A1 = "sampleA", A2 = "sampleA", B1 = "sampleB", B2 = "sampleB", C1 = NA)[colnames(reference)]
  reference[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = ref[, c("A1", "B1")])
  list(spatial = spatial, reference = reference)
}

capture_deconv_inputs <- function(method, pair, ...) {
  captured <- NULL
  capture <- function(...) {
    captured <<- list(...)
    stop("deconvolution input capture complete")
  }
  testthat::local_mocked_bindings(
    rctd_run_spacexr = capture,
    card_run_backend = capture,
    spotlight_run_backend = capture,
    check_r = function(...) TRUE,
    .package = "scop"
  )
  args <- c(list(object = pair$spatial, reference = pair$reference,
    reference_label = "celltype", verbose = FALSE), list(...))
  if (identical(method, "RunRCTD")) args$min_cells <- args$min_cells %||% 1
  if (identical(method, "RunCARD")) args$sample_varname <- "sample"
  expect_error(do.call(get(method), args), "deconvolution input capture complete")
  captured
}

test_that("deconvolution wrappers align subset assay cells by ID", {
  pair <- make_deconv_alignment_pair()
  original <- pair
  for (method in c("RunRCTD", "RunCARD", "RunSPOTlight")) {
    captured <- capture_deconv_inputs(method, pair,
      assay = "ALT", reference_assay = "ALT", features = c("g2", "g1"))
    expect_identical(colnames(captured$ref_counts), c("A1", "B1"))
    expect_identical(colnames(captured$st_counts), c("s1", "s3"))
    expect_identical(rownames(captured$ref_counts), c("g2", "g1"))
    expect_identical(rownames(captured$st_counts), c("g2", "g1"))
    if (identical(method, "RunCARD")) {
      expect_identical(rownames(captured$ref_meta), c("A1", "B1"))
      expect_identical(captured$ref_meta$.scop_cell_type, c("A", "B"))
      expect_identical(captured$ref_meta$.scop_sample, c("sampleA", "sampleB"))
      expect_identical(captured$ct_select, c("A", "B"))
    } else {
      labels <- if (identical(method, "RunRCTD")) captured$ref_labels else captured$labels
      expect_identical(names(labels), c("A1", "B1"))
      expect_identical(as.character(labels), c("A", "B"))
    }
    if (!identical(method, "RunSPOTlight")) {
      expect_identical(rownames(captured$coords), c("s1", "s3"))
      expect_equal(captured$coords$x, c(1, 3))
    }
  }
  expect_identical(pair, original)
})

test_that("RCTD preserves full selected-assay library depth before feature filtering", {
  pair <- make_deconv_alignment_pair()
  captured <- capture_deconv_inputs("RunRCTD", pair,
    assay = "ALT", reference_assay = "ALT", features = c("g2", "g1"))
  expect_equal(captured$st_numi, c(s1 = 1000, s3 = 3000))
  expect_equal(captured$ref_numi, c(A1 = 1000, B1 = 2000))
  expect_equal(Matrix::colSums(captured$st_counts), c(s1 = 30, s3 = 30))
  expect_equal(Matrix::colSums(captured$ref_counts), c(A1 = 10, B1 = 10))

  # Non-shared genes also contribute to full library depth when features = NULL.
  ref <- GetAssayData5(pair$reference, assay = "ALT")
  rownames(ref)[3] <- "reference-only"
  pair$reference[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = ref)
  captured <- capture_deconv_inputs("RunRCTD", pair, assay = "ALT", reference_assay = "ALT")
  expect_identical(rownames(captured$ref_counts), c("g1", "g2"))
  expect_equal(captured$st_numi, c(s1 = 1000, s3 = 3000))
  expect_equal(captured$ref_numi, c(A1 = 1000, B1 = 2000))
})

test_that("missing labels and selected-feature zero counts are filtered after alignment", {
  pair <- make_deconv_alignment_pair()
  # Full library depth is positive, but B1 and s3 have no selected-gene counts.
  for (which in c("reference", "spatial")) {
    counts <- GetAssayData5(pair[[which]], assay = "ALT")
    counts[c("g1", "g2"), if (which == "reference") "B1" else "s3"] <- 0
    pair[[which]][["ALT"]] <- SeuratObject::CreateAssay5Object(counts = counts)
  }
  for (method in c("RunRCTD", "RunCARD", "RunSPOTlight")) {
    captured <- capture_deconv_inputs(method, pair,
      assay = "ALT", reference_assay = "ALT", features = c("g1", "g2"))
    expect_identical(colnames(captured$ref_counts), "A1")
    expect_identical(colnames(captured$st_counts), "s1")
    if (identical(method, "RunRCTD")) {
      expect_identical(names(captured$ref_labels), "A1")
      expect_equal(captured$ref_numi, c(A1 = 1000))
      expect_equal(captured$st_numi, c(s1 = 1000))
    } else if (identical(method, "RunCARD")) {
      expect_identical(captured$ct_select, "A")
      expect_identical(captured$ref_meta$.scop_cell_type, "A")
    } else {
      expect_identical(names(captured$labels), "A1")
      expect_identical(levels(captured$labels), "A")
    }
  }

  pair <- make_deconv_alignment_pair()
  pair$reference$celltype[c(1, 5)] <- NA
  for (method in c("RunRCTD", "RunCARD", "RunSPOTlight")) {
    captured <- capture_deconv_inputs(method, pair, reference_assay = "ALT")
    expect_identical(colnames(captured$ref_counts), "B1")
  }
  pair$reference$celltype[3] <- NA
  for (method in c("RunRCTD", "RunCARD", "RunSPOTlight")) {
    expect_error(do.call(get(method), c(list(pair$spatial, pair$reference,
      reference_label = "celltype", reference_assay = "ALT", verbose = FALSE),
      if (identical(method, "RunRCTD")) list(min_cells = 1))),
      "non-missing class")
  }
})

test_that("RCTD applies min_cells to selected-assay and zero-count survivors", {
  pair <- make_deconv_alignment_pair()
  expect_error(RunRCTD(pair$spatial, pair$reference, reference_label = "celltype",
    reference_assay = "ALT", min_cells = 2, verbose = FALSE), "min_cells")
  ref <- GetAssayData5(pair$reference, assay = "RNA")
  ref[c("g1", "g2"), "B2"] <- 0
  pair$reference[["ALT"]] <- NULL
  pair$reference[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = ref)
  captured <- capture_deconv_inputs("RunRCTD", pair,
    reference_assay = "ALT", min_cells = 2, features = c("g1", "g2"))
  expect_identical(colnames(captured$ref_counts), c("A1", "A2"))
  expect_identical(names(captured$ref_labels), c("A1", "A2"))
  expect_equal(captured$ref_numi, c(A1 = 1000, A2 = 1000))
})

test_that("ordinary full-assay inputs retain their order and library depths", {
  pair <- make_deconv_alignment_pair()
  for (method in c("RunRCTD", "RunCARD", "RunSPOTlight")) {
    captured <- capture_deconv_inputs(method, pair)
    expect_identical(colnames(captured$ref_counts), c("A1", "A2", "B1", "B2"))
    expect_identical(colnames(captured$st_counts), c("s1", "s2", "s3"))
    if (identical(method, "RunRCTD")) {
      expect_equal(captured$ref_numi, Matrix::colSums(captured$ref_counts))
      expect_equal(captured$st_numi, Matrix::colSums(captured$st_counts))
    }
  }
})

test_that("deconvolution rejects ambiguous or unmatched identities", {
  pair <- make_deconv_alignment_pair()
  counts <- GetAssayData5(pair$reference, assay = "ALT")
  for (ids in list(NULL, c("A1", "A1"), c("A1", NA), c("A1", ""))) {
    bad <- counts
    colnames(bad) <- ids
    expect_error(deconv_align_reference(bad, pair$reference, "celltype"),
      "unique, non-missing cell IDs")
  }
  colnames(counts) <- c("B1", "unknown")
  expect_error(deconv_align_reference(counts, pair$reference, "celltype"),
    "missing cell IDs")
  expect_error(deconv_check_cell_order(counts, rev(colnames(counts)), "Reference labels"),
    "cell IDs and order must exactly match")
})


test_that("reference alignment follows matrix order independently of metadata order", {
  pair <- make_deconv_alignment_pair()
  counts <- GetAssayData5(pair$reference, assay = "ALT")[, c("B1", "A1")]
  aligned <- deconv_align_reference(counts, pair$reference, "celltype", verbose = FALSE)
  expect_identical(colnames(aligned$counts), c("B1", "A1"))
  expect_identical(names(aligned$labels), c("B1", "A1"))
  expect_identical(as.character(aligned$labels), c("B", "A"))

  # Seurat's assay setter canonicalizes assay order to the object. Reorder the
  # initial count matrices, rather than assuming a supplied assay permutation survives.
  pair <- make_deconv_alignment_pair(reverse_order = TRUE)
  for (method in c("RunRCTD", "RunCARD", "RunSPOTlight")) {
    captured <- capture_deconv_inputs(method, pair, assay = "ALT", reference_assay = "ALT")
    expect_identical(colnames(captured$ref_counts), c("B1", "A1"))
    expect_identical(colnames(captured$st_counts), c("s3", "s1"))
    if (identical(method, "RunCARD")) {
      expect_identical(captured$ref_meta$.scop_cell_type, c("B", "A"))
      expect_identical(captured$ref_meta$.scop_sample, c("sampleB", "sampleA"))
    } else {
      labels <- if (identical(method, "RunRCTD")) captured$ref_labels else captured$labels
      expect_identical(names(labels), c("B1", "A1"))
      expect_identical(as.character(labels), c("B", "A"))
    }
  }
})


test_that("RCTD depth uses the selected assay and its documented rounding policy", {
  pair <- make_deconv_alignment_pair()
  for (which in c("reference", "spatial")) {
    counts <- GetAssayData5(pair[[which]], assay = "ALT")
    counts["other", ] <- counts["other", ] + 500.2
    pair[[which]][["ALT"]] <- SeuratObject::CreateAssay5Object(counts = counts)
  }
  captured <- capture_deconv_inputs("RunRCTD", pair,
    assay = "ALT", reference_assay = "ALT", features = c("g1", "g2"))
  expect_equal(captured$st_numi, c(s1 = 1500, s3 = 3500))
  expect_equal(captured$ref_numi, c(A1 = 1500, B1 = 2500))
  expect_error(RunRCTD(pair$spatial, pair$reference, reference_label = "celltype",
    assay = "ALT", reference_assay = "ALT", features = c("g1", "g2"),
    min_cells = 1, round_counts = FALSE, verbose = FALSE), "non-integer")
})

test_that("RCTD records each min_cells removal once across filtering stages", {
  pair <- make_deconv_alignment_pair()
  pair$reference$celltype[5] <- "C"
  counts <- GetAssayData5(pair$reference, assay = "RNA")
  counts[c("g1", "g2"), "B2"] <- 0
  pair$reference[["ALT"]] <- NULL
  pair$reference[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = counts)
  testthat::local_mocked_bindings(rctd_run_spacexr = function(st_counts, ref_labels, ...) {
    list(weights = matrix(1, ncol(st_counts), 1,
      dimnames = list(colnames(st_counts), unique(as.character(ref_labels)))),
      metadata = NULL, api = "input contract adapter", object = NULL)
  }, .package = "scop")
  out <- RunRCTD(pair$spatial, pair$reference, reference_label = "celltype",
    reference_assay = "ALT", features = c("g1", "g2"), min_cells = 2, verbose = FALSE)
  expect_equal(out@tools$RCTD$dropped_cell_types,
    data.frame(cell_type = c("C", "B"), n_cells = c(1L, 1L)))
  expect_identical(rownames(out@tools$RCTD$proportions), colnames(pair$spatial))
})
