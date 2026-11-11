test_that("fragment plotting checks files and preserves native length data", {
  skip_if_not_installed("Rsamtools")
  object <- make_atac_plot_test_object()
  expect_error(ATACFragmentPlot(object), "explicit non-empty")
  expect_error(ATACFragmentPlot(object, region = "chr1-1-1000"), "mapped fragment files")
  directory <- tempfile("atac_plot_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  plain <- file.path(directory, "fragments.tsv")
  writeLines(c(
    "chr1\t10\t50\tA\t4", "chr1\t20\t200\tB\t1",
    "chr1\t30\t400\tC\t2", "chr1\t40\t90\tD\t1"
  ), plain)
  path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(path, format = "bed")
  Signac::Fragments(object[["peaks"]]) <- Signac::CreateFragmentObject(path,
    cells = stats::setNames(colnames(object), colnames(object)), verbose = FALSE
  )
  before <- serialize(object, NULL)
  actual <- ATACFragmentPlot(object, region = "chr1-1-1000", group.by = "group")
  expect_identical(serialize(object, NULL), before)
  native <- Signac::FragmentHistogram(object, region = "chr1-1-1000", group.by = "group")
  expect_equal(actual$data, native$data)
  expect_equal(sort(actual$data$length), c(40, 50, 180, 370))
  expect_s3_class(ATACFragmentPlot(object, region = "chr1-1-1000", cells = c("A", "C")), "ggplot")
  expect_error(ATACFragmentPlot(object, region = "chr2-1-1000"), "No mapped fragments")
  expect_error(ATACFragmentPlot(object, region = "chr1-1-1000", log.scale = NA), "TRUE or FALSE")
})


test_that("fragment histograms preserve per-file barcode maps and ignore inactive resources", {
  skip_if_not_installed("Rsamtools")
  object <- make_atac_plot_test_object()
  directory <- tempfile("atac_plot_maps_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  maps <- list(c(A = "raw1", B = "raw2"), c(C = "raw1", D = "raw2"))
  records <- list(
    c("chr1\t10\t50\traw1\t7", "chr1\t20\t200\traw2\t1"),
    c("chr1\t30\t400\traw1\t5", "chr1\t40\t90\traw2\t1")
  )
  fragments <- lapply(seq_along(maps), function(i) {
    plain <- file.path(directory, paste0("sample", i, ".tsv"))
    writeLines(records[[i]], plain)
    path <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"), overwrite = TRUE)
    Rsamtools::indexTabix(path, format = "bed")
    Signac::CreateFragmentObject(path, cells = maps[[i]], verbose = FALSE)
  })
  Signac::Fragments(object[["peaks"]]) <- fragments
  p <- ATACFragmentPlot(object, "chr1-1-1000", group.by = "group")
  expect_equal(sort(p$data$length), c(40, 50, 180, 370))
  expect_equal(sort(unique(as.character(p$data$cell))), c("A", "B", "C", "D"))
  object[["peaks"]]@fragments[[2]]@path <- file.path(directory, "absent.tsv.gz")
  p <- ATACFragmentPlot(object, "chr1-1-1000", cells = c("B", "A"))
  expect_equal(sort(p$data$length), c(40, 180))
  expect_error(ATACFragmentPlot(object, "chr1-1-1000"), "accessible local")
})

test_that("fragment plot input errors are actionable", {
  expect_error(ATACFragmentPlot(list()), "Seurat")
  object <- make_atac_plot_test_object()
  expect_error(ATACFragmentPlot(object, assay = "missing"), "ChromatinAssay")
  expect_error(ATACFragmentPlot(object, group.by = "missing"), "metadata column")
  object$group[1] <- NA
  expect_error(ATACFragmentPlot(object, group.by = "group"), "non-missing group")
})
