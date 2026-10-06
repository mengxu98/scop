make_mistyr_audit_input <- function() {
  x <- matrix(1:15, 3, dimnames = list(c("Gene-1", "Gene.1", "Gene3"), paste0("c", 1:5)))
  s <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(x, sparse = TRUE))
  s$x <- c(0, 1, 2, 0, 1); s$y <- c(0, 0, 0, 1, 1)
  SeuratObject::SetAssayData(s, layer = "data", new.data = GetAssayData5(s, layer = "counts"))
}
test_that("MISTy maps original target names without name collisions", {
  map <- setNames(c("Gene-1", "Gene.1", "Other"), c("Gene.1", "Gene.1.1", "Other"))
  expect_identical(mistyr_resolve_targets(c("Gene-1", "Gene.1"), map), c("Gene.1", "Gene.1.1"))
  expect_null(mistyr_resolve_targets(NULL, map))
  expect_identical(mistyr_resolve_targets(c(2L, 1L), map), c(2L, 1L))
  expect_error(mistyr_resolve_targets("absent", map), "target_subset")
})

test_that("MISTy rejects nonfinite observations instead of imputing zeros", {
  srt <- make_mistyr_audit_input()
  for (bad in c(NA_real_, Inf, -Inf)) {
    expr <- GetAssayData5(srt, assay = "RNA", layer = "data")
    expr[1, 1] <- bad
    srt <- SeuratObject::SetAssayData(srt, assay = "RNA", layer = "data", new.data = expr)
    expect_error(mistyr_prepare_input(srt, "RNA", "data", features = rownames(srt)[1:3],
      coord.cols = c("x", "y"), coordinate_space = "raw"), "finite")
  }
})

test_that("MISTy optional sample boundary rejects combined coordinate fields", {
  srt <- make_mistyr_audit_input()
  srt$sample <- c("A", "A", "A", "B", "B")
  expect_error(RunMistyR(srt, features = rownames(srt)[1:3], coord.cols = c("x", "y"),
    sample.by = "sample", verbose = FALSE), "single sample")
})


test_that("RunMistyR passes mapped targets to the backend and retains original provenance", {
  srt <- make_mistyr_audit_input()
  testthat::local_mocked_bindings(
    check_r = function(...) invisible(TRUE),
    get_namespace_fun = function(package, name) switch(name,
      create_initial_view = function(expression) list(intraview = list(data = expression)),
      add_paraview = function(...) stop("unexpected spatial view"),
      add_juxtaview = function(...) stop("unexpected spatial view"),
      run_misty = function(views, target.subset, results.folder, ...) {
        expect_identical(target.subset, c("Gene.1.1", "Gene.1"))
        expect_true(all(target.subset %in% colnames(views$intraview$data)))
        results.folder
      },
      collect_results = function(...) list()
    )
  )
  srt$sample <- "one"
  out <- RunMistyR(srt, features = rownames(srt), coord.cols = c("x", "y"),
    views = NULL, target_subset = c("Gene-1", "Gene.1"), sample.by = "sample", verbose = FALSE)
  expect_identical(out@tools$MistyR$parameters$target_subset, c("Gene-1", "Gene.1"))
})
