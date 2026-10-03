make_irea_reference <- function() {
  mat <- matrix(0,
    nrow = 4, ncol = 12,
    dimnames = list(c("G1", "G2", "G3", "G4"), paste0("cell", seq_len(12)))
  )
  mat["G1", ] <- c(rep(1, 4), rep(30, 4), rep(2, 4))
  mat["G2", ] <- c(rep(10, 4), rep(1, 4), rep(10, 4))
  mat["G3", ] <- 2
  mat["G4", ] <- c(rep(2, 8), rep(30, 4))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(mat, sparse = TRUE))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object$sample <- rep(c("PBS", "IL12", "IL15"), each = 4)
  object$polarization <- c(rep("None", 8), rep("NK-e", 4))
  structure(list(
    object = object,
    cytokine = data.frame(Cytokine_Str = c("IL12", "IL15"), Gene = c("G1", "G4"), FDR = c(0.001, 0.001)),
    polarization = data.frame(Polarization = "NK-e", Gene = "G4", P_adj = 0.001),
    cell_type = "NK_cell", species = "Mouse", paths = character(), checksum_md5 = character(),
    provenance = "synthetic test"
  ), class = "irea_reference")
}

test_that("gene scores and projection keep their response direction", {
  r <- make_irea_reference()
  a <- RunIREA("G1", reference = r)$table
  expect_gt(a$effect[a$term == "IL12"], 0)
  expect_lt(a$effect[a$term == "IL15"], a$effect[a$term == "IL12"])
  p <- RunIREA(c(G1 = 1, G2 = -1), reference = r, gene_diff_cutoff = 0)$table
  expect_gt(p$effect[p$term == "IL12"], 0)
  n <- RunIREA(c(G1 = -1, G2 = 1), reference = r, gene_diff_cutoff = 0)$table
  expect_lt(n$effect[n$term == "IL12"], 0)
  expect_equal(p$fdr, stats::p.adjust(p$p_value, "BH"))
})

test_that("hypergeometric and polarization retain separate meanings", {
  r <- make_irea_reference()
  h <- RunIREA("G1", reference = r, method = "hypergeometric")$table
  expect_equal(h$effect[h$term == "IL12"], 1)
  expect_equal(h$effect[h$term == "IL15"], 0)
  s <- RunIREA("G4", reference = r, analysis = "cell_polarization")$table
  expect_gt(s$effect, 0)
  expect_s3_class(IREAPlot(s <- RunIREA("G4", reference = r, analysis = "cell_polarization"), "radar"), "ggplot")
})

test_that("invalid input fails instead of returning a biological zero", {
  r <- make_irea_reference()
  expect_error(RunIREA("UNKNOWN", reference = r), "No input genes match")
  expect_error(RunIREA(c(G1 = 0, G2 = 0), reference = r), "nonzero")
  expect_error(RunIREA("G1", reference = r, matrix = c(G1 = 1)), "Unused")
  expect_error(RunIREA(c(G1 = 1, G1 = 2), reference = r), "duplicate")
})

test_that("Seurat contrasts use the same projection core", {
  r <- make_irea_reference()
  o <- r$object
  out <- RunIREA(o,
    reference = r, group.by = "sample", case = "IL12", control = "PBS",
    gene_diff_cutoff = 0
  )
  expect_s4_class(out, "Seurat")
  expect_s3_class(out@tools$IREA, "irea_result")
  expect_true(all(is.finite(out@tools$IREA$table$effect)))
})

test_that("human orthologue mapping is explicit and conflicts are omitted", {
  r <- make_irea_reference()
  r$species <- "Human"
  r$cytokine$Gene_Human <- c("H1", "H4")
  r$polarization$Gene_Human <- "H4"
  out <- RunIREA("H1", reference = r)
  expect_equal(out$matched_genes, "G1")
  expect_equal(out$parameters$species, "Human")
  expect_error(RunIREA("BAD", reference = r), "No unambiguous")
})

test_that("portal-style matrix file selects a named contrast", {
  r <- make_irea_reference()
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data.frame(
    gene = c("G1", "G2"), case1 = c(1, -1),
    case2 = c(-1, 1)
  ), path, row.names = FALSE)
  expect_error(RunIREA(path, reference = r, gene_diff_cutoff = 0), "Select a contrast")
  a <- RunIREA(path, reference = r, contrast = "case1", gene_diff_cutoff = 0)
  b <- RunIREA(path, reference = r, contrast = "case2", gene_diff_cutoff = 0)
  expect_equal(a$table$effect, -b$table$effect)
})

test_that("no positively enriched polarization has a zero radar score", {
  r <- make_irea_reference()
  out <- RunIREA("G3", reference = r, analysis = "cell_polarization")
  expect_true(all(out$table$radar_score == 0))
  expect_true(all(out$table$effect <= 0))
  expect_error(PrepareIREA(db = "unknown", data_dir = tempdir()), "Unsupported")
})

test_that("multiple results can be compared without dropping terms", {
  r <- make_irea_reference()
  a <- RunIREA("G1", reference = r)
  b <- RunIREA("G4", reference = r)
  plot <- IREAPlot(list(first = a, second = b), "heatmap")
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), nrow(a$table) + nrow(b$table))
})

test_that("sparse and dense reference layers give the same scores", {
  sparse <- make_irea_reference()
  dense <- sparse
  data <- as.matrix(SeuratObject::LayerData(dense$object, assay = "RNA", layer = "data"))
  SeuratObject::LayerData(dense$object, assay = "RNA", layer = "data") <- data
  a <- RunIREA(c(G1 = 1, G2 = -1), reference = sparse, gene_diff_cutoff = 0)
  b <- RunIREA(c(G1 = 1, G2 = -1), reference = dense, gene_diff_cutoff = 0)
  expect_equal(a$table$effect, b$table$effect)
  expect_equal(a$table$p_value, b$table$p_value)
})

test_that("polarization compares against PBS rather than unpolarized treated cells", {
  r <- make_irea_reference()
  d <- SeuratObject::LayerData(r$object, assay = "RNA", layer = "data")
  score <- as.numeric(d["G1", ])
  actual <- RunIREA("G1", reference = r, analysis = "cell_polarization")$table
  target <- score[r$object$polarization == "NK-e"]
  control <- score[r$object$sample == "PBS"]
  expect_equal(actual$effect, mean(target) - mean(control))
  expect_equal(actual$p_value, stats::wilcox.test(target, control, exact = FALSE)$p.value)
  expect_equal(actual$n_control, 4L)
  expect_false(isTRUE(all.equal(
    actual$effect,
    mean(target) - mean(score[r$object$polarization == "None"])
  )))
  r$object$sample[r$object$sample == "PBS"] <- "Other"
  expect_error(RunIREA("G1", reference = r, analysis = "cell_polarization"), "PBS baseline")
})

test_that("Seurat group labels are aligned to cells in the requested layer", {
  r <- make_irea_reference()
  o <- r$object
  counts <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts")
  selected <- c(1, 2, 5, 6)
  o[["partial"]] <- SeuratObject::CreateAssay5Object(counts = counts[, selected])
  actual <- RunIREA(o, reference = r, group.by = "sample", case = "IL12", control = "PBS", assay = "partial", layer = "counts", gene_diff_cutoff = 0)@tools$IREA
  expected <- Matrix::rowMeans(counts[, 5:6]) - Matrix::rowMeans(counts[, 1:2])
  expect_equal(actual$table, RunIREA(expected, reference = r, gene_diff_cutoff = 0)$table)
  expect_error(
    RunIREA(o, reference = r, group.by = "sample", case = "IL15", control = "PBS", assay = "partial", layer = "counts"),
    "selected assay/layer"
  )
  o$sample[3] <- NA_character_
  actual <- RunIREA(o, reference = r, group.by = "sample", case = "IL12", control = "PBS", layer = "counts", gene_diff_cutoff = 0)@tools$IREA
  expected <- Matrix::rowMeans(counts[, 5:8]) - Matrix::rowMeans(counts[, c(1, 2, 4)])
  expect_equal(actual$table, RunIREA(expected, reference = r, gene_diff_cutoff = 0)$table)
  expect_error(RunIREA(o, reference = r, group.by = "sample", case = "PBS", control = "PBS"), "distinct")
})

test_that("comparison plots reject incompatible effect scales", {
  r <- make_irea_reference()
  a <- RunIREA("G1", reference = r)
  b <- RunIREA("G1", reference = r, method = "hypergeometric")
  expect_error(IREAPlot(list(a = a, b = b), "heatmap"), "same method")
  expect_error(IREAPlot(list(a = a, b = b), "dotplot"), "same method")
})

test_that("radar thresholds are recomputed without changing result statistics", {
  r <- RunIREA("G4", reference = make_irea_reference(), analysis = "cell_polarization")
  r$table <- data.frame(
    term = c("A", "B", "C"), effect = c(2, 1, -1),
    fdr = c(0.03, 0.2, 0.8), radar_score = c(1, 0.5, 0)
  )
  before <- r$table
  loose <- IREAPlot(r, "radar", padjustCutoff = 0.05)$layers[[3]]$data
  strict <- IREAPlot(r, "radar", padjustCutoff = 0.01)$layers[[3]]$data
  expect_equal(IREAPlot(r, plot_type = "radar", padjustCutoff = 0.05)$layers[[3]]$data, loose)
  expect_error(IREAPlot(r, plot_type = "radar", padjustCutoff = NA_real_), "padjustCutoff")
  expect_equal(loose$radar_score, c(1, 0.5, 0))
  expect_equal(strict$radar_score, c(0, 0, 0))
  expect_equal(r$table, before)
  expect_equal(irea_radar_score(c(NA, 1), c(NA, 0.03)), c(0, 1))
})

test_that("missing list genes and nonfinite cutoffs are handled explicitly", {
  r <- make_irea_reference()
  expect_equal(RunIREA(c(NA, "G1", ""), reference = r)$matched_genes, "G1")
  expect_error(RunIREA(NA_character_, reference = r), "No usable")
  expect_error(RunIREA("G1", reference = r, gene_diff_cutoff = Inf), "nonnegative")
})

test_that("split layer names cannot silently select the first batch", {
  r <- make_irea_reference()
  counts <- SeuratObject::LayerData(r$object, assay = "RNA", layer = "counts")
  o <- SeuratObject::CreateSeuratObject(counts = list(first = counts[, 1:6], second = counts[, 7:12]))
  o$sample <- r$object$sample
  expect_error(
    RunIREA(o, reference = r, group.by = "sample", case = "IL12", control = "PBS", layer = "counts"),
    "exact name"
  )
})

test_that("multiple matrix columns share the specified BH family", {
  r <- make_irea_reference()
  input <- data.frame(gene = c("G1", "G2", "G3"), a = c(1, -1, 0), b = c(0, 0, 1))
  a <- RunIREA(input, reference = r, contrast = "a", gene_diff_cutoff = 0)
  b <- RunIREA(input, reference = r, contrast = "b", gene_diff_cutoff = 0)
  expected <- stats::p.adjust(c(a$table$p_value, b$table$p_value), "BH")
  expect_equal(c(a$table$fdr, b$table$fdr), expected)
  expect_equal(a$parameters$contrasts_adjusted, c("a", "b"))
  single <- RunIREA(input,
    reference = r, contrast = "a", gene_diff_cutoff = 0,
    fdr_scope = "selected_contrast"
  )
  expect_equal(single$table$fdr, stats::p.adjust(single$table$p_value, "BH"))
  expect_equal(single$table$p_value, a$table$p_value)
  expect_equal(single$table$effect, a$table$effect)
  expect_equal(a$parameters$fdr_family, "reference_terms_across_all_supplied_contrasts")
})

test_that("Wilcoxon extreme tails remain positive and symmetric", {
  a <- 101:200
  b <- 1:100
  expected <- 2 * stats::pnorm(-(5000 - 0.5) / sqrt(100 * 100 * 201 / 12))
  actual <- irea_wilcox(a, b)
  expect_gt(actual, 0)
  expect_equal(actual / expected, 1, tolerance = 1e-12)
  expect_equal(irea_wilcox(b, a) / expected, 1, tolerance = 1e-12)
  expect_equal(
    irea_wilcox(c(1, 1, 2, 3), c(1, 2, 2, 4)),
    2 * stats::wilcox.test(c(1, 1, 2, 3), c(1, 2, 2, 4),
      alternative = "less", exact = FALSE, correct = TRUE, digits.rank = Inf
    )$p.value
  )
  expect_equal(irea_wilcox(rep(1, 5), rep(1, 8)), 1)
  expect_true(is.na(irea_wilcox(numeric(), 1:3)))
})


test_that("factor contrasts cannot silently change numerical values", {
  r <- make_irea_reference()
  numeric <- data.frame(gene = c("G1", "G2"), comparison = c(1, -1))
  factored <- numeric
  factored$comparison <- factor(c("1", "-1"))
  expect_error(RunIREA(factored, reference = r), "must be numeric")
  text <- numeric
  text$comparison <- c("1", "-1")
  expect_error(RunIREA(text, reference = r), "must be numeric")
  expected <- RunIREA(c(G1 = 1, G2 = -1), reference = r, gene_diff_cutoff = 0)
  expect_equal(RunIREA(numeric, reference = r, gene_diff_cutoff = 0)$table, expected$table)
})

test_that("selected contrasts and Seurat design remain traceable", {
  r <- make_irea_reference()
  input <- data.frame(gene = c("G1", "G2"), a = c(1, -1), b = c(-1, 1))
  a <- RunIREA(input, reference = r, contrast = "a", gene_diff_cutoff = 0)
  b <- RunIREA(input, reference = r, contrast = "b", gene_diff_cutoff = 0)
  expect_identical(a$parameters[["contrast", exact = TRUE]], "a")
  expect_identical(b$parameters[["contrast", exact = TRUE]], "b")
  expect_false(identical(a$parameters, b$parameters))
  o <- r$object
  o@tools$Other <- list(marker = TRUE)
  result <- RunIREA(o,
    reference = r, group.by = "sample", case = "IL12", control = "PBS",
    assay = "RNA", layer = "data", gene_diff_cutoff = 0
  )
  expect_identical(result@tools$Other, o@tools$Other)
  expect_identical(
    result@tools$IREA$parameters[c("group.by", "case", "control", "assay", "layer")],
    list(group.by = "sample", case = "IL12", control = "PBS", assay = "RNA", layer = "data")
  )
})

test_that("plots use custom colours and themes without changing numerical results", {
  r <- RunIREA("G1", reference = make_irea_reference())
  plot <- IREAPlot(r,
    palcolor = c(Positive = "red", Negative = "blue", `Not significant` = "grey"),
    theme_use = ggplot2::theme_classic()
  )
  built <- ggplot2::ggplot_build(plot)
  expect_true(all(built$data[[1]]$fill %in% c("red", "blue", "grey")))
  expect_equal(plot$data$effect, r$table$effect)
})

test_that("IREA preparation reads cached sources and detects changed contents", {
  skip_if_not_installed("R.cache")
  skip_if_not_installed("readxl")
  skip_if_not_installed("openxlsx")
  r <- make_irea_reference()
  directory <- tempfile("irea-cache-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  saveRDS(r$object, file.path(directory, "downloadableData_ligands-seurat-NK_cell.RDS"))
  openxlsx::write.xlsx(
    list(NK_cell = r$cytokine),
    file.path(directory, "dataFiles_SuppTable3_Cytokine_Signatures.xlsx")
  )
  openxlsx::write.xlsx(
    list(`NK cell` = r$polarization),
    file.path(directory, "dataFiles_SuppTable7_Polarization_Signatures.xlsx")
  )
  before <- getOption("timeout")
  first <- PrepareDB(
    db = "IREA_NK_cell", species = "Mus_musculus",
    data_dir = directory, verbose = FALSE
  )[["Mus_musculus"]][["IREA_NK_cell"]]
  expect_s3_class(first, "irea_reference")
  expect_equal(getOption("timeout"), before)
  expect_equal(first$cytokine, r$cytokine)
  expect_equal(first$polarization, r$polarization)
  manual <- PrepareIREA(species = "Mus_musculus", data_dir = directory, db = "IREA_NK_cell")[["Mus_musculus"]][["IREA_NK_cell"]]
  expect_equal(RunIREA("G1", reference = first)$table, RunIREA("G1", reference = manual)$table)
  selected <- PrepareDB(
    db = "IREA_NK_cell", species = "Mus_musculus",
    data_dir = list(IREA_NK_cell = directory, IREA = file.path(directory, "unused")),
    verbose = FALSE
  )[["Mus_musculus"]][["IREA_NK_cell"]]
  expect_equal(selected$paths, first$paths)
  manifest <- utils::read.csv(file.path(directory, "reference_manifest.csv"))
  expect_equal(nrow(manifest), 3L)
  expect_equal(manifest$checksum_md5, unname(tools::md5sum(first$paths)))
  second <- PrepareDB(
    db = "IREA_NK_cell", species = "Mus_musculus",
    data_dir = list(IREA = directory), verbose = FALSE
  )[["Mus_musculus"]][["IREA_NK_cell"]]
  expect_equal(RunIREA("G1", reference = first)$table, RunIREA("G1", reference = second)$table)
  annotation <- list(
    TERM2GENE = data.frame(Term = "T", symbol = "G1"),
    TERM2NAME = data.frame(Term = "T", Name = "Example"), version = "test"
  )
  local_mocked_bindings(
    list_db_cache_entries = function(species, db, exact_db = FALSE) {
      expect_true(all(db %in% c("CachedAnnotation", "IREA_NK_cell")))
      cached <- data.frame(DB = "CachedAnnotation", Species = species, timestamp = as.POSIXct("2026-01-01", tz = "UTC"), file = "cached-annotation")
      cached[cached$DB %in% db, , drop = FALSE]
    },
    .package = "scop"
  )
  local_mocked_bindings(
    readCacheHeader = function(pathname, ...) {
      list(comment = "test|Mus_musculus-CachedAnnotation", timestamp = as.POSIXct("2026-01-01", tz = "UTC"))
    },
    loadCache = function(pathname, ...) annotation,
    .package = "R.cache"
  )
  mixed <- PrepareDB(
    db = c("CachedAnnotation", "IREA_NK_cell"), species = "Mus_musculus",
    db_IDtypes = "symbol", data_dir = list(IREA = directory), verbose = FALSE
  )[["Mus_musculus"]]
  expect_setequal(names(mixed), c("CachedAnnotation", "IREA_NK_cell"))
  expect_equal(mixed[["CachedAnnotation"]], annotation)
  expect_equal(mixed[["IREA_NK_cell"]]$paths, first$paths)
  expect_equal(RunIREA("G1", reference = mixed[["IREA_NK_cell"]])$table, RunIREA("G1", reference = first)$table)
  model <- list(data_dir = tempdir(), files = "model_parameters.rds", version = "1.1.0")
  local_mocked_bindings(PrepareCytoTRACE2 = function(db_update = FALSE, verbose = TRUE) list(CytoTRACE2 = model), .package = "scop")
  combined <- PrepareDB(
    db = c("CytoTRACE2", "CachedAnnotation", "IREA_NK_cell"), species = "Mus_musculus",
    db_IDtypes = "symbol", data_dir = list(IREA = directory), verbose = FALSE
  )
  expect_identical(names(combined), c("CytoTRACE2", "Mus_musculus"))
  expect_identical(combined$CytoTRACE2, model)
  expect_equal(combined$Mus_musculus, mixed)
  nested <- tempfile("irea-directory-")
  dir.create(file.path(nested, "IREA_NK_cell"), recursive = TRUE)
  nested <- normalizePath(nested, winslash = "/", mustWork = TRUE)
  on.exit(unlink(nested, recursive = TRUE), add = TRUE)
  file.copy(first$paths, nested)
  file.copy(first$paths[["object"]], file.path(nested, "IREA_NK_cell"))
  nested_reference <- PrepareDB(
    db = "IREA_NK_cell", species = "Mus_musculus", data_dir = nested, verbose = FALSE
  )[["Mus_musculus"]][["IREA_NK_cell"]]
  expect_equal(dirname(nested_reference$paths[["object"]]), file.path(nested, "IREA_NK_cell"))
  expect_equal(dirname(nested_reference$paths[["cytokine"]]), nested)
  expect_equal(nested_reference$cytokine, first$cytokine)
  expect_equal(RunIREA("G1", reference = nested_reference)$table, RunIREA("G1", reference = first)$table)
  manifest <- utils::read.csv(file.path(directory, "reference_manifest.csv"))
  r$object$sample[1] <- "Changed"
  saveRDS(r$object, first$paths[["object"]])
  expect_error(PrepareDB(
    db = "IREA_NK_cell", species = "Mus_musculus",
    data_dir = directory, verbose = FALSE
  ), "checksum mismatch")
  expect_equal(getOption("timeout"), before)
  expect_equal(utils::read.csv(file.path(directory, "reference_manifest.csv")), manifest)
})

test_that("IREA database selectors reject ambiguous or unused options", {
  skip_if_not_installed("R.cache")
  expect_error(PrepareDB(db = "IREA", verbose = FALSE), "IREA_<cell type>")
  expect_error(PrepareDB(db = "IREA_unknown", verbose = FALSE), "Unsupported IREA database selector")
  expect_error(PrepareDB(db = "IREA_NK_cell", cell_type = "NK_cell", verbose = FALSE), "Unused")
  expect_error(PrepareDB(db = c("IREA_NK_cell", NA_character_), verbose = FALSE), "missing")
  expect_error(PrepareDB(db = "IREA_NK_cell", species = "Danio_rerio", verbose = FALSE), "supports")
  expect_error(PrepareDB(
    db = "IREA_NK_cell", data_dir = list(Other = tempdir()),
    verbose = FALSE
  ), "path for IREA_NK_cell")
})

test_that("IREA selectors preserve complete identifiers and named reference results", {
  skip_if_not_installed("R.cache")
  skip_if_not_installed("readxl")
  skip_if_not_installed("openxlsx")
  r <- make_irea_reference()
  directory <- tempfile("irea-selection-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  for (type in c("NK_cell", "T_cell_CD4")) {
    saveRDS(r$object, file.path(directory, paste0("downloadableData_ligands-seurat-", type, ".RDS")))
  }
  other <- r$cytokine
  other$Gene <- c("G2", "G3")
  openxlsx::write.xlsx(
    list(NK_cell = r$cytokine, T_cell_CD4 = other),
    file.path(directory, "dataFiles_SuppTable3_Cytokine_Signatures.xlsx")
  )
  openxlsx::write.xlsx(
    list(`NK cell` = r$polarization, `CD4+ T cell` = r$polarization),
    file.path(directory, "dataFiles_SuppTable7_Polarization_Signatures.xlsx")
  )
  databases <- c("IREA_T_cell_CD4", "IREA_NK_cell")
  prepared <- PrepareDB(
    db = databases, species = "Mus_musculus",
    data_dir = list(IREA = directory), verbose = FALSE
  )[["Mus_musculus"]]
  expect_identical(names(prepared), databases)
  expect_identical(prepared[["IREA_T_cell_CD4"]]$cell_type, "T_cell_CD4")
  expect_identical(prepared[["IREA_NK_cell"]]$cell_type, "NK_cell")
  expect_equal(prepared[["IREA_T_cell_CD4"]]$cytokine, other)
  expect_equal(prepared[["IREA_NK_cell"]]$cytokine, r$cytokine)
  expect_equal(
    RunIREA("G1", reference = prepared[["IREA_NK_cell"]])$table,
    RunIREA("G1", reference = r)$table
  )
})
