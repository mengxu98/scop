# Proposed hosted acceptance entry. Run only after this candidate has been
# installed by the existing pkgdown local::. dependency step into this R library.
# This script does not install dependencies or substitute mocked acceptance.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[[1]], mustWork = TRUE)
output <- normalizePath(args[[2]], mustWork = TRUE)
setwd(repo)
required <- c(
  "scop", "Seurat", "SeuratObject", "Signac", "Rsamtools", "GenomicRanges",
  "IRanges", "Matrix", "withr", "testthat", "pkgload"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing acceptance dependencies: ", paste(missing, collapse = ", "))
versions <- data.frame(
  package = required,
  version = vapply(required, function(x) as.character(utils::packageVersion(x)), ""),
  library = vapply(required, find.package, "")
)
print(versions)
utils::write.csv(versions, file.path(output, "atac-package-versions.csv"), row.names = FALSE)
# Verify local::. installed the candidate rather than a same-version cached
# package. Compare function syntax without evaluating repository expressions.
checked <- character()
for (file in c("R/StandardWorkflow.R", "R/integration-wnn.R", "R/integration.R", "R/CallATACPeaks.R")) {
  for (expression in parse(file.path(repo, file), keep.source = FALSE)) {
    if (!is.call(expression) || !identical(expression[[1L]], as.name("<-")) ||
      !is.symbol(expression[[2L]]) || !is.call(expression[[3L]]) ||
      !identical(expression[[3L]][[1L]], as.name("function"))) next
    name <- as.character(expression[[2L]])
    installed <- get(name, envir = asNamespace("scop"), inherits = FALSE)
    expected <- expression[[3L]]
    stopifnot(
      identical(deparse(formals(installed)), deparse(expected[[2L]])),
      identical(deparse(body(installed)), deparse(expected[[3L]]))
    )
    checked <- c(checked, name)
  }
}
stopifnot(length(checked) > 0L)
writeLines(checked, file.path(output, "installed-candidate-functions.txt"))
results <- testthat::test_local(
  path = repo,
  filter = "^atac-(workflow-completion|wnn|peak-calling)$",
  reporter = "summary",
  load_package = "installed",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)
summary <- as.data.frame(results)
utils::write.csv(summary[, setdiff(names(summary), "result")],
  file.path(output, "atac-focused-tests.csv"), row.names = FALSE
)
stopifnot(
  nrow(summary) > 0L,
  !any(summary$failed > 0L), !any(summary$error), !any(summary$skipped),
  setequal(c("test-atac-workflow-completion.R", "test-atac-wnn.R", "test-atac-peak-calling.R"), basename(summary$file))
)
# The unit peak-call test intentionally mocks the external process. A separate
# real MACS2 check is required before claiming that boundary has been exercised.
# This extra case tests the cross-PR representation-record/retention contract.
input <- new.env(parent = emptyenv())
utils::data("pbmcmultiome_sub", package = "scop", envir = input)
paired <- subset(input$pbmcmultiome_sub, cells = colnames(input$pbmcmultiome_sub)[seq_len(80)])
result <- scop::WNN_integrate(
  srt_merge = paired, assay = c("RNA", "peaks"),
  linear_reduction_dims = 10, linear_reduction_dims_use = 2:8,
  nHVF = 100, neighbor_k = 10, nonlinear_reduction_dims = 2, verbose = FALSE
)
stopifnot(all(result@misc$WNN_reduction_list %in% SeuratObject::Reductions(result)))
receipt <- SeuratObject::Misc(result[["ATAClsi"]], "scop_atac_representation")
embedding <- SeuratObject::Embeddings(result[["ATAClsi"]])
stopifnot(
  is.list(receipt), length(receipt) > 0L,
  identical(receipt$dimensions, 2:8),
  identical(receipt$cells, rownames(embedding)),
  isTRUE(all.equal(unname(receipt$depth_correlation),
    as.vector(stats::cor(embedding, receipt$depth))))
)
saveRDS(list(receipt = receipt, reductions = SeuratObject::Reductions(result)),
  file.path(output, "atac-wnn-combined-contract.rds")
)
cat("Cross-PR WNN retention and non-null ATAC representation record: PASS\n")
cat("Focused unit tests include mocked MACS2; real-caller acceptance is separate.\n")
