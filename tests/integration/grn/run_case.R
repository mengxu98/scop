arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) %% 2L != 0L || !length(arguments)) stop("Use --case NAME --output DIRECTORY, optionally --package scop --library DIRECTORY --request request.rds.")
options <- stats::setNames(as.list(arguments[seq.int(2L, length(arguments), 2L)]), sub("^--", "", arguments[seq.int(1L, length(arguments), 2L)]))
if (is.null(options$case) || is.null(options$output)) stop("--case and --output are required.")
if (any(!names(options) %in% c("case", "output", "package", "library", "request"))) stop("Unknown integration driver option.")
if (!is.null(options$library)) .libPaths(c(options$library, .libPaths()))
package <- if (is.null(options$package)) "scop" else options$package
if (!identical(package, "scop")) stop("This driver validates the scop GRN interface.")
suppressPackageStartupMessages(library(package, character.only = TRUE))
if (dir.exists(options$output) && length(list.files(options$output, all.files = TRUE, no.. = TRUE))) stop("Choose a new empty evidence directory.")
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
result <- tryCatch({
  if (!is.null(options$request)) {
    request <- readRDS(options$request)
    if (!is.list(request) || is.null(request$object) || !is.list(request$args)) stop("A request must contain object and args.")
    request$args$return_data <- TRUE
    value <- do.call(scop::RunGRN, c(list(object = request$object), request$args))
  } else {
    case_files <- c("legacy-contract" = "legacy-contract.R", "rna-tutorial" = "rna-tutorial.R")
    if (!options$case %in% names(case_files)) stop("This case needs an explicit --request RDS file.")
    script_argument <- grep("^--file=", commandArgs(), value = TRUE)
    script <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
    environment <- new.env(parent = globalenv())
    sys.source(file.path(dirname(script), "cases", case_files[[options$case]]), envir = environment)
    value <- environment$run_case(options$output)
  }
  if (!inherits(value, "GRNResult")) stop("The case did not return GRNResult.")
  saveRDS(value, file.path(options$output, "result.rds"))
  for (group in names(value$groups)) {
    layers <- scop::GetGRN(value, group = group, what = "layers")
    for (layer in names(layers)) if (is.data.frame(layers[[layer]])) {
      utils::write.table(layers[[layer]], file.path(options$output, paste0(sprintf("group-%04d", match(group, names(value$groups))), "-", layer, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  list(status = "passed", method = value$method, implementation = value$implementation,
    groups = names(value$groups), result_sha256 = digest::digest(file = file.path(options$output, "result.rds"), algo = "sha256"))
}, error = function(error) list(status = "failed", error = conditionMessage(error)))
result$case <- options$case
result$started_at <- as.character(started)
result$finished_at <- as.character(Sys.time())
result$package <- list(version = as.character(utils::packageVersion("scop")), path = find.package("scop"))
result$session <- utils::capture.output(utils::sessionInfo())
saveRDS(result, file.path(options$output, "receipt.rds"))
jsonlite::write_json(result, file.path(options$output, "receipt.json"), auto_unbox = TRUE, pretty = TRUE)
if (!identical(result$status, "passed")) stop(result$error, call. = FALSE)
cat("GRN case passed:", options$case, "\n")
