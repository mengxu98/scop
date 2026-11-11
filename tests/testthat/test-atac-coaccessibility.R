test_that("coaccessibility validates reference, embedding and result ownership before fitting", {
  x <- make_atac_test_object()
  expect_error(RunATACCoaccessibility(x, c(chr2 = 1000000), "toy", dims = 1:2, k = 2), "cover source peaks")
  expect_error(RunATACCoaccessibility(x, c(chr1 = 1000000), "wrong", dims = 1:2, k = 2), "genome|reference")
  expect_error(RunATACCoaccessibility(x, c(chr1 = 1000000), "toy", dims = 3, k = 2), "dims")
  expect_error(RunATACCoaccessibility(x, c(chr1 = 1000000), "toy", dims = 1:2, k = 4), "k < cell")
  expect_error(RunATACCoaccessibility(x, c(chr1 = 1000000), "toy", dims = 1:2, k = 2, window = 100000), "window >")
  expect_error(RunATACCoaccessibility(x, c(chr1 = 1000000, `1` = 1000000), "toy", dims = 1:2, k = 2), "collide")
  SeuratObject::Misc(x[["peaks"]], "scop_coaccessibility") <- list(existing = TRUE)
  expect_error(RunATACCoaccessibility(x, c(chr1 = 1000000), "toy", dims = 1:2, k = 2), "unused")
})

test_that("coaccessibility dispatch preserves sparse binary data, order and native NA scores", {
  x <- make_atac_test_object()
  original <- serialize(x, NULL)
  calls <- new.env()
  fake_new <- function(expression_data, cell_metadata, gene_metadata) {
    calls$counts <- expression_data
    calls$genes <- gene_metadata
    expression_data
  }
  fake_aggregate <- function(cds, reduced_coordinates, k, size_factor_normalize, silent, return_agg_info) {
    calls$coords <- reduced_coordinates
    expect_true(return_agg_info)
    list(cds[, 1:2, drop = FALSE], data.frame(cell = colnames(cds), bin = c(1, 1, 2, 2)))
  }
  fake_run <- function(cds, genomic_coords, window, sample_num, silent) {
    data.frame(Peak1 = rownames(cds)[1], Peak2 = rownames(cds)[2], coaccess = NA_real_)
  }
  old <- get_namespace_fun
  local_mocked_bindings(check_r = function(...) invisible(TRUE), get_namespace_fun = function(package, fun) {
    if (package == "monocle3") {
      return(fake_new)
    }
    if (package == "cicero" && fun == "make_cicero_cds") {
      return(fake_aggregate)
    }
    if (package == "cicero" && fun == "run_cicero") {
      return(fake_run)
    }
    old(package, fun)
  })
  original_description <- utils::packageDescription
  original_version <- utils::packageVersion
  local_mocked_bindings(
    packageDescription = function(pkg, ...) if (pkg == "cicero") list(Imports = "monocle3") else original_description(pkg, ...),
    packageVersion = function(pkg, ...) if (pkg %in% c("cicero", "monocle3")) numeric_version("1.0.0") else original_version(pkg, ...),
    .package = "utils"
  )
  set.seed(123)
  rng <- .Random.seed
  out <- RunATACCoaccessibility(x, c(chr1 = 1000000), "toy", dims = 1:2, k = 2)
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(x, NULL), original)
  expect_true(inherits(calls$counts, "sparseMatrix"))
  expect_true(all(calls$counts@x == 1))
  expect_identical(rownames(calls$coords), colnames(x))
  expect_identical(calls$genes$bp1, c(1000L, 2000L, 3000L))
  receipt <- SeuratObject::Misc(out[["peaks"]], "scop_coaccessibility")
  expect_true(is.na(receipt$connections$coaccess))
  expect_identical(SeuratObject::DefaultAssay(out), SeuratObject::DefaultAssay(x))
  expect_identical(Signac::Links(out[["peaks"]]), Signac::Links(x[["peaks"]]))
})
