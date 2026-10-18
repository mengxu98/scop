#!/usr/bin/env Rscript
# Proposed real-process acceptance. This script has not been executed here.
# Usage: Rscript --vanilla acceptance-macs2.R /absolute/path/to/macs2 NEW_OUTPUT_DIR
# Uses an already-installed candidate scop package. Installs nothing.
# All input data are synthetic toy fixtures. A pass is software/process evidence,
# not validation of PBMC data, biological reproducibility, or biological claims.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: acceptance-macs2.R MACS2_EXECUTABLE NEW_OUTPUT_DIR", call. = FALSE)
}
out <- arguments[[2L]]
if (!nzchar(out) || file.exists(out) ||
    isTRUE(nzchar(Sys.readlink(out), keepNA = TRUE)) ||
    !dir.exists(dirname(out))) {
  stop("Output must be a new path within an existing directory.", call. = FALSE)
}
if (!dir.create(out)) stop("Cannot create output directory.", call. = FALSE)
out <- normalizePath(out, mustWork = TRUE)

summary <- list(
  status = "RUNNING", stage = "preflight",
  scope = "Synthetic toy fragments; real MACS2 process and software consistency only",
  biological_validation = FALSE, skipped = 0L,
  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  expected_versions = c(Signac = "1.17.1", MACS2 = "2.2.9.1"),
  versions = character(), positive_calls = list(), negative_call = NULL,
  warnings = character(), error = NULL, git = character()
)

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

# Snapshot substantive input values, excluding native fastmatch lookup caches.
snapshot <- function(object) {
  assay <- object[["peaks"]]
  matrix_value <- function(layer) {
    matrix <- SeuratObject::GetAssayData(assay, layer = layer)
    list(class = class(matrix), dim = dim(matrix),
         dimnames = lapply(dimnames(matrix), as.vector),
         i = matrix@i, p = matrix@p, x = matrix@x)
  }
  metadata <- object[[]]
  metadata_values <- lapply(metadata, function(x) {
    if (is.character(x)) as.vector(x) else x
  })
  feature_metadata <- assay[[]]
  rownames(feature_metadata) <- as.vector(rownames(feature_metadata))
  list(
    counts = matrix_value("counts"), data = matrix_value("data"),
    cells = as.vector(colnames(object)), assays = as.vector(names(object)),
    default_assay = SeuratObject::DefaultAssay(object),
    metadata = metadata_values, metadata_rows = as.vector(rownames(metadata)),
    ranges = as.data.frame(Signac::granges(assay)),
    genome = GenomeInfoDb::genome(Signac::granges(assay)),
    seqlengths = GenomeInfoDb::seqlengths(Signac::granges(assay)),
    feature_metadata = feature_metadata, misc = SeuratObject::Misc(assay),
    fragments = lapply(Signac::Fragments(assay), function(fragment) {
      cells <- Signac::GetFragmentData(fragment, slot = "cells")
      list(path = Signac::GetFragmentData(fragment, slot = "path"),
           names = as.vector(names(cells)), barcodes = as.vector(cells),
           hash = as.vector(Signac::GetFragmentData(fragment, slot = "hash")))
    })
  )
}

run_acceptance <- function() {
  revision <- system2("git", c("rev-parse", "HEAD", "HEAD^{tree}"),
                      stdout = TRUE, stderr = TRUE)
  check(is.null(attr(revision, "status")) && length(revision) == 2L &&
          all(grepl("^[0-9a-f]{40}$", revision)),
        "Cannot identify the checked-out source commit and tree.")
  summary$git <<- setNames(revision, c("commit", "tree"))
  required <- c("scop", "Signac", "Seurat", "SeuratObject", "Matrix",
                "Rsamtools", "GenomicRanges", "GenomeInfoDb", "S4Vectors")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  check(all(available), paste("Missing required packages:",
                              paste(required[!available], collapse = ", ")))
  summary$versions <<- vapply(required, function(package) {
    as.character(utils::packageVersion(package))
  }, character(1))
  summary$scop_library <<- find.package("scop")
  check(identical(unname(summary$versions[["Signac"]]), "1.17.1"),
        "This acceptance requires Signac 1.17.1; no version substitution is allowed.")
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

  check(file.exists(arguments[[1L]]) && !dir.exists(arguments[[1L]]),
        "MACS2_EXECUTABLE must identify an existing executable file.")
  macs <- normalizePath(arguments[[1L]], mustWork = TRUE)
  check(file.access(macs, mode = 1L) == 0L, "MACS2 executable is not executable.")
  version <- system2(macs, "--version", stdout = TRUE, stderr = TRUE)
  writeLines(version, file.path(out, "macs2-version.txt"))
  summary$macs2_executable <<- macs
  summary$macs2_executable_md5 <<- unname(tools::md5sum(macs))
  check(is.null(attr(version, "status")) &&
          identical(trimws(paste(version, collapse = " ")), "macs2 2.2.9.1"),
        "This acceptance requires a working MACS2 2.2.9.1 executable.")

  cells <- c("toycap1-A", "toycap1-B", "toycap2-A", "toycap2-B")
  metadata <- data.frame(
    capture = rep(c("toy/captureA", "toy_captureA"), each = 2L),
    group = rep(c("toy/group1", "toy_group1"), 2L),
    row.names = cells, stringsAsFactors = FALSE
  )
  mappings <- list(setNames(c("A", "B"), cells[1:2]),
                   setNames(c("A", "B"), cells[3:4]))
  ranges <- GenomicRanges::GRanges(c("toyChr1:9901-10400", "toyChr1:49901-50400"))
  GenomeInfoDb::seqlengths(ranges) <- c(toyChr1 = 1000000L)
  counts <- Matrix::sparseMatrix(
    i = c(1L, 1L, 2L, 2L), j = 1:4, x = 100,
    dims = c(2L, 4L), dimnames = list(c("toyPeak1", "toyPeak2"), cells)
  )
  write.table(metadata, file.path(out, "toy-cell-metadata.tsv"), sep = "\t",
              quote = FALSE, col.names = NA)
  writeLines(c(
    "Synthetic toy fixture only; no human or PBMC data.",
    "Reference: scop-toy-1M; one contig toyChr1, length 1000000.",
    "Two capture files, each with 200 records and reused raw barcodes A/B.",
    "Four capture-by-group calls, each with 100 distinct start positions.",
    "Fifth-column values alternate 1 and 7; input counts are record counts.",
    "BED; --nomodel; --extsize 200; --shift -100; -g 1000000; -q 0.05.",
    "No --keep-dup override: MACS2 default duplicate handling is retained.",
    "Negative fixture adds INVALID as column six to request a real parser error.",
    "If Signac rejects the negative fixture before MACS2 runs, acceptance fails."
  ), file.path(out, "fixture-manifest.txt"))

  make_fixture <- function(directory, invalid_strand = FALSE) {
    check(dir.create(directory), "Cannot create a fresh fixture directory.")
    paths <- character(2L)
    original <- setNames(vector("list", 2L), unique(metadata$capture))
    for (i in 1:2) {
      start <- c(10000L, 50000L)[i] + 0:199
      table <- data.frame("toyChr1", start, start + 50L,
                          rep(c("A", "B"), 100L), rep(c(1L, 7L), 100L),
                          stringsAsFactors = FALSE)
      if (invalid_strand) table[[6L]] <- "INVALID"
      plain <- file.path(directory, paste0("toy-capture", i, ".tsv"))
      write.table(table, plain, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      original[[i]] <- readLines(plain)
      paths[[i]] <- Rsamtools::bgzip(plain, dest = paste0(plain, ".gz"))
      Rsamtools::indexTabix(paths[[i]], format = "bed")
    }
    source_paths <- c(paths, paste0(paths, ".tbi"))
    check(all(file.exists(source_paths)), "Fixture fragment or index is missing.")
    checksums <- tools::md5sum(source_paths)
    write.table(data.frame(path = names(checksums), md5 = unname(checksums)),
                file.path(directory, "input-md5.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    object <- scop::ReadATACData(
      counts, ranges = ranges, metadata = metadata, genome = "scop-toy-1M",
      fragments = paths, fragment.cells = mappings, verbose = FALSE
    )
    list(object = object, original = original,
         source_paths = source_paths, checksums = checksums)
  }

  call_arguments <- function(bed, id, directory) {
    c("callpeak", "-t", bed, "-f", "BED", "--nomodel",
      "--extsize", 200L, "--shift", -100L, "-g", 1e6, "-q", 0.05,
      "-n", id, "--outdir", directory)
  }

  summary$stage <<- "positive fixture construction"
  fixture <- make_fixture(file.path(out, "positive toy fixture"))
  before <- snapshot(fixture$object)
  saveRDS(before, file.path(out, "positive-input-before.rds"))
  summary$stage <<- "positive real MACS2 wrapper calls"
  called <- scop::CallATACPeaks(
    fixture$object, sample.by = "capture", group.by = "group", cells = cells,
    macs2.path = macs, outdir = file.path(out, "wrapped toy calls"),
    effective.genome.size = 1e6, qvalue = 0.05,
    extsize = 200L, shift = -100L, verbose = FALSE
  )
  after <- snapshot(fixture$object)
  saveRDS(after, file.path(out, "positive-input-after.rds"))
  check(identical(before, after), "Positive call changed substantive input values.")
  check(identical(fixture$checksums, tools::md5sum(fixture$source_paths)),
        "Positive call changed input fragment/index files.")

  expected_groups <- data.frame(call_id = sprintf("unit%04d", 1:4),
                                sample = metadata$capture, group = metadata$group,
                                stringsAsFactors = FALSE)
  check(identical(called$groups, expected_groups), "Capture/group identities differ.")
  check(identical(names(called$peaks), expected_groups$call_id) &&
          identical(names(called$receipt$commands), expected_groups$call_id),
        "Peak or command identifiers differ from expected call identifiers.")
  check(identical(as.vector(called$receipt$cells), cells) &&
          identical(called$receipt$sample.by, "capture") &&
          identical(called$receipt$group.by, "group"),
        "Receipt does not retain selected cells and grouping fields.")
  check(identical(unname(called$receipt$backend[["Signac"]]), "1.17.1") &&
          identical(unname(called$receipt$backend[["MACS2"]]), "macs2 2.2.9.1"),
        "Receipt backend versions differ from preflight.")

  summary$stage <<- "independent split membership and direct caller comparison"
  for (i in 1:4) {
    group <- expected_groups[i, ]
    id <- group$call_id
    selected <- rownames(metadata)[metadata$capture == group$sample &
                                   metadata$group == group$group]
    check(identical(as.vector(names(called$receipt$membership)[
      called$receipt$membership == id]), selected),
      paste("Wrong cell membership for", id))
    bed <- file.path(called$receipt$outdir, "fragments", paste0(id, ".bed"))
    original <- fixture$original[[group$sample]]
    barcode <- if (group$group == "toy/group1") "A" else "B"
    raw_barcodes <- vapply(strsplit(original, "\t", fixed = TRUE),
                           function(row) row[[4L]], character(1))
    expected_lines <- original[raw_barcodes == barcode]
    check(length(expected_lines) == 100L && file.exists(bed) &&
            identical(readLines(bed), expected_lines),
          paste("Split BED differs from original capture/barcode records for", id))
    command <- called$receipt$commands[[id]]
    expected_arguments <- call_arguments(bed, id, called$receipt$outdir)
    check(identical(command$executable, macs) && command$status == 0L &&
            identical(command$arguments, expected_arguments),
          paste("Unexpected executable, arguments or wrapper status for", id))
    check(file.exists(command$log) && file.exists(command$output),
          paste("Missing wrapper log or narrowPeak for", id))
    check(identical(readLines(file.path(called$receipt$outdir,
                                       paste0(id, ".arguments.txt"))),
                    c(macs, expected_arguments)),
          paste("Argument evidence differs for", id))

    direct_directory <- file.path(out, paste0("direct toy ", id))
    check(dir.create(direct_directory), "Cannot create direct-call directory.")
    direct_arguments <- call_arguments(bed, id, direct_directory)
    writeLines(c(macs, direct_arguments), file.path(direct_directory, "arguments.txt"))
    direct_log <- file.path(direct_directory, "caller.log")
    direct_status <- system2(macs, shQuote(direct_arguments),
                             stdout = direct_log, stderr = direct_log)
    writeLines(as.character(direct_status), file.path(direct_directory, "exit-status.txt"))
    direct_output <- file.path(direct_directory, paste0(id, "_peaks.narrowPeak"))
    summary$positive_calls[[id]] <<- list(
      wrapper_status = command$status, direct_status = direct_status,
      split_records = length(expected_lines), wrapper_output = command$output,
      direct_output = direct_output
    )
    check(direct_status == 0L && file.exists(direct_output),
          paste("Independent MACS2 call failed or omitted output for", id))
    check(file.info(command$output)$size > 0 && file.info(direct_output)$size > 0,
          paste("Empty positive fixture output is not acceptance for", id))
    check(identical(readLines(command$output), readLines(direct_output)),
          paste("Wrapper and independent MACS2 narrowPeak differ for", id))

    raw <- utils::read.delim(command$output, header = FALSE, stringsAsFactors = FALSE)
    peaks <- called$peaks[[id]]
    check(ncol(raw) == 10L && !anyNA(raw) && nrow(raw) == length(peaks) &&
            length(peaks) > 0L, paste("Malformed or empty positive peak set for", id))
    check(all(as.character(GenomicRanges::seqnames(peaks)) == raw[[1L]]) &&
            all(GenomicRanges::start(peaks) == raw[[2L]] + 1L) &&
            all(GenomicRanges::end(peaks) == raw[[3L]]) &&
            all(GenomeInfoDb::genome(peaks) == "scop-toy-1M"),
          paste("Imported coordinates or genome differ for", id))
    columns <- c(peak_name = 4L, score = 5L, fold_change = 7L,
                 neg_log10pvalue = 8L, neg_log10qvalue = 9L, summit_offset = 10L)
    for (column in names(columns)) {
      check(isTRUE(all.equal(S4Vectors::mcols(peaks)[[column]], raw[[columns[[column]]]],
                             check.attributes = FALSE)),
            paste("Imported narrowPeak field differs:", id, column))
    }
    summary$positive_calls[[id]]$peaks <<- length(peaks)
    summary$positive_calls[[id]]$narrowPeak_md5 <<- unname(tools::md5sum(command$output))
  }
  check(file.exists(file.path(called$receipt$outdir, "commands.rds")) &&
          file.exists(file.path(called$receipt$outdir, "peak_calls.rds")),
        "Successful call did not retain commands.rds and peak_calls.rds.")

  summary$stage <<- "negative toy fixture construction (must reach real caller)"
  negative <- make_fixture(file.path(out, "negative toy fixture"), invalid_strand = TRUE)
  negative_before <- snapshot(negative$object)
  saveRDS(negative_before, file.path(out, "negative-input-before.rds"))
  failed_directory <- file.path(out, "failed toy calls")
  summary$stage <<- "negative real MACS2 wrapper call"
  failure <- tryCatch({
    scop::CallATACPeaks(
      negative$object, sample.by = "capture", group.by = "group", cells = cells,
      macs2.path = macs, outdir = failed_directory,
      effective.genome.size = 1e6, qvalue = 0.05,
      extsize = 200L, shift = -100L, verbose = FALSE
    )
    NULL
  }, error = identity)
  negative_after <- snapshot(negative$object)
  saveRDS(negative_after, file.path(out, "negative-input-after.rds"))
  error_text <- if (inherits(failure, "error")) conditionMessage(failure) else "NO ERROR"
  writeLines(error_text, file.path(out, "negative-error.txt"))
  check(identical(negative_before, negative_after),
        "Negative call changed substantive input values.")
  check(identical(negative$checksums, tools::md5sum(negative$source_paths)),
        "Negative call changed input fragment/index files.")
  check(inherits(failure, "error") && grepl("MACS2 failed", error_text, fixed = TRUE),
        "Negative fixture did not produce a real MACS2-call failure; preflight rejection is not a pass.")
  commands_path <- file.path(failed_directory, "commands.rds")
  check(file.exists(commands_path), "Negative caller exit evidence is missing.")
  failed_commands <- readRDS(commands_path)
  attempted <- Filter(Negate(is.null), failed_commands)
  check(length(attempted) == 1L && identical(names(attempted), "unit0001"),
        "Expected exactly the first negative caller attempt before stopping.")
  failed <- attempted[[1L]]
  failed_bed <- file.path(failed_directory, "fragments", "unit0001.bed")
  check(identical(failed$executable, macs) &&
          identical(failed$arguments, call_arguments(failed_bed, "unit0001", failed_directory)) &&
          length(failed$status) == 1L && !is.na(failed$status) && failed$status != 0L &&
          file.exists(failed$log) && file.exists(failed_bed),
        "Negative test lacks the expected executable, command, real nonzero exit or files.")
  check(any(grepl("StrandFormatError", readLines(failed$log), fixed = TRUE)),
        "Negative MACS2 log does not contain the expected real BED parser failure.")
  check(!file.exists(file.path(failed_directory, "peak_calls.rds")),
        "Failed caller incorrectly retained a successful result receipt.")
  summary$negative_call <<- list(status = failed$status, log = failed$log,
                                 error = error_text, real_parser_error = TRUE)
  summary$stage <<- "all required toy real-process assertions completed"
  invisible(NULL)
}

failure <- tryCatch(
  withCallingHandlers(run_acceptance(), warning = function(warning) {
    summary$warnings <<- c(summary$warnings, conditionMessage(warning))
  }),
  error = identity
)
summary$status <- if (inherits(failure, "error")) "FAIL" else "PASS"
summary$error <- if (inherits(failure, "error")) conditionMessage(failure) else NULL
summary$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
diagnostic_failure <- tryCatch({
  saveRDS(summary, file.path(out, "summary.rds"))
  capture.output(dput(summary), file = file.path(out, "summary.txt"))
  NULL
}, error = identity)
cat(summary$status, ": ", summary$scope, "\n", sep = "")
cat("Stage: ", summary$stage, "\nEvidence: ", out, "\n", sep = "")
if (!is.null(summary$error)) cat("Error: ", summary$error, "\n", sep = "")
if (inherits(failure, "error")) {
  if (inherits(diagnostic_failure, "error")) {
    message("Could not save acceptance diagnostics: ", conditionMessage(diagnostic_failure))
  }
  stop(failure)
}
if (inherits(diagnostic_failure, "error")) stop(diagnostic_failure)
# A successful run reaches ordinary Rscript shutdown without an explicit exit.
