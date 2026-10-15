make_atac_motif_fixture <- function(n = 40L) {
  skip_if_not_installed("Signac")
  skip_if_not_installed("GenomicRanges")
  counts <- Matrix::Matrix(matrix(1L, n, 3L, dimnames = list(
    paste0("input", seq_len(n)), c("A-1", "B-1", "C-1")
  )), sparse = TRUE)
  ranges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = "chr1",
    start = seq(1L, by = 200L, length.out = n), end = seq(100L, by = 200L, length.out = n)
  ))
  names(ranges) <- rownames(counts)
  object <- ReadATACData(counts, ranges = ranges, genome = "hg38", verbose = FALSE)
  incidence <- Matrix::Matrix(matrix(0, n, 3L,
    dimnames = list(rownames(object), c("M1", "M2", "M3"))
  ), sparse = TRUE)
  incidence[c(1, 2, 3, 6), "M1"] <- 1
  incidence[, "M2"] <- 1
  pwm <- stats::setNames(rep(list(matrix(1, 4, 3,
    dimnames = list(c("A", "C", "G", "T"), NULL)
  )), 3), colnames(incidence))
  motif <- Signac::CreateMotifObject(
    data = incidence, pwm = pwm,
    motif.names = list(M1 = "toy_enriched", M2 = "toy_all", M3 = "toy_none")
  )
  Signac::Motifs(object[["peaks"]]) <- motif
  metadata <- data.frame(
    GC.percent = seq(25, 70, length.out = n),
    sequence.length = rep(100, n), row.names = rownames(object)
  )
  object[["peaks"]] <- SeuratObject::AddMetaData(object[["peaks"]], metadata)
  object
}

test_that("motif population and native hypergeometric denominators are explicit", {
  object <- make_atac_motif_fixture(20L)
  foreground <- rownames(object)[1:5]
  population <- rownames(object)
  expect_warning(actual <- RunMotifEnrichment(object, foreground,
    assay = "peaks",
    universe = population, background = population, verbose = FALSE
  ), "small number")
  expect_warning(native <- Signac::FindMotifs(object, foreground,
    assay = "peaks",
    background = population, p.adjust.method = "BH", verbose = FALSE
  ), "small number")
  expect_equal(actual$results, native)
  # M=20 population peaks; m=4 motif-positive; n=16 negative; k=5 queried;
  # observed=3. This is an exact finite-population enumeration, not another fit.
  exact <- (choose(4, 3) * choose(16, 2) + choose(4, 4) * choose(16, 1)) / choose(20, 5)
  table <- actual$results
  expect_equal(table$pvalue[table$motif == "M1"], exact, tolerance = 1e-14)
  expect_equal(table$observed[table$motif == "M1"], 3)
  expect_equal(table$background[table$motif == "M1"], 4)
  expect_equal(table$pvalue[table$motif %in% c("M2", "M3")], rep(1, 2))
  expect_equal(table$p.adjust, stats::p.adjust(table$pvalue, "BH"))
  expect_identical(actual$receipt$foreground, foreground)
  expect_identical(actual$receipt$background, population)
  expect_identical(actual$receipt$reference, population[-(1:5)])
})

test_that("native matched references are recorded separately from full population", {
  skip_if_not_installed("withr")
  object <- make_atac_motif_fixture()
  foreground <- rownames(object)[11:20]
  universe <- rownames(object)[1:35]
  metadata <- object[["peaks"]][[]]
  seed <- 17L
  expected <- withr::with_seed(seed, Signac::MatchRegionStats(
    meta.feature = metadata[setdiff(universe, foreground), , drop = FALSE],
    query.feature = metadata[foreground, , drop = FALSE],
    features.match = c("GC.percent", "sequence.length"), n = 15, verbose = FALSE
  ))
  set.seed(29)
  state <- .Random.seed
  before <- serialize(object, NULL)
  actual <- RunMotifEnrichment(object, foreground,
    assay = "peaks", universe = universe,
    n.background = 15, seed = seed, verbose = FALSE
  )
  expect_identical(.Random.seed, state)
  expect_identical(serialize(object, NULL), before)
  expect_identical(actual$receipt$reference, expected)
  expect_identical(actual$receipt$background, c(foreground, expected))
  expect_true(all(actual$receipt$background %in% universe))
  expect_equal(actual$results, Signac::FindMotifs(object, foreground,
    assay = "peaks",
    background = c(foreground, expected), p.adjust.method = "BH", verbose = FALSE
  ))
})

test_that("motif enrichment rejects ambiguous universes and stale incidence order", {
  object <- make_atac_motif_fixture()
  ids <- rownames(object)
  expect_error(RunMotifEnrichment(object, rep(ids[1], 2), background = ids), "unique foreground")
  expect_error(RunMotifEnrichment(object, ids[1:5], universe = ids[-1]), "within universe")
  expect_error(RunMotifEnrichment(object, ids[1:5], background = ids[-1]), "include every foreground")
  expect_error(RunMotifEnrichment(object, ids[1:5], background = c(ids, ids[1])), "unique population")
  object[["peaks"]]@motifs@data <- object[["peaks"]]@motifs@data[rev(ids), , drop = FALSE]
  expect_error(RunMotifEnrichment(object, ids[1:5], background = ids), "aligned exactly")
})

test_that("motif matching does not silently drop invalid metadata", {
  object <- make_atac_motif_fixture()
  ids <- rownames(object)
  metadata <- object[["peaks"]][[]]
  metadata$GC.percent[30] <- NA_real_
  object[["peaks"]] <- SeuratObject::AddMetaData(object[["peaks"]], metadata)
  expect_error(RunMotifEnrichment(object, ids[1:10], n.background = 10), "finite numeric")
  expect_error(RunMotifEnrichment(object, ids[1:10], match.features = "missing"), "metadata are missing")
  expect_error(RunMotifEnrichment(object, ids[1:10], seed = NA_real_), "[Ss]eed")
  expect_error(RunMotifEnrichment(object, ids[1:10], n.background = 0), "positive integer")
})

test_that("native peak-name caches do not invalidate aligned motif incidence", {
  skip_if_not_installed("fastmatch")
  object <- make_atac_motif_fixture()
  ids <- as.vector(rownames(object))
  cached <- paste0(ids, "")
  fastmatch::fmatch(cached[1], cached)
  expect_false(is.null(attr(cached, ".match.hash")))
  object[["peaks"]]@motifs@data@Dimnames[[1L]] <- cached
  expect_false(identical(rownames(object[["peaks"]]@motifs@data), rownames(object)))
  actual <- RunMotifEnrichment(object, ids[1:10], background = ids, verbose = FALSE)
  expect_equal(actual$results, Signac::FindMotifs(object, ids[1:10],
    background = ids, p.adjust.method = "BH", verbose = FALSE
  ))
  object[["peaks"]]@motifs@data@Dimnames[[1L]] <- rev(ids)
  expect_error(RunMotifEnrichment(object, ids[1:10], background = ids), "aligned exactly")
})

test_that("pattern sparse motif incidence retains native enrichment semantics", {
  object <- make_atac_motif_fixture()
  object[["peaks"]]@motifs@data <- as(object[["peaks"]]@motifs@data, "nMatrix")
  ids <- rownames(object)
  expected <- Signac::FindMotifs(object, ids[1:10], background = ids, verbose = FALSE)
  actual <- RunMotifEnrichment(object, ids[1:10], background = ids, verbose = FALSE)
  expect_equal(actual$results, expected)
})

test_that("matching restores an absent RNG seed and preserves zero-hit native output", {
  withr::local_preserve_seed()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  object <- make_atac_motif_fixture()
  ids <- rownames(object)
  actual <- RunMotifEnrichment(object, ids[1:10], n.background = 20L, verbose = FALSE)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  zero <- actual$results[actual$results$motif == "M3", , drop = FALSE]
  expect_equal(zero$pvalue, 1)
  expect_true(is.nan(zero$fold.enrichment))
  expect_equal(actual$results$p.adjust, stats::p.adjust(actual$results$pvalue, "BH"))
})
