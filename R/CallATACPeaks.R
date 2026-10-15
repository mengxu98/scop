#' Call MACS2 peaks for explicitly selected capture groups
#'
#' @md
#' @description
#' Split selected cells with native Signac, run the installed MACS2 executable,
#' and retain separate peak sets, command arguments and logs. Peak inference is
#' performed by MACS2; no caller is installed or silently substituted.
#' @param object A validated Seurat object containing a ChromatinAssay and local
#'   indexed fragments with explicit cell mappings.
#' @param sample.by Metadata column identifying capture/library IDs.
#' @param group.by Optional metadata column subdividing each capture, for example
#'   a cell type. Calls are capture-by-group, never pooled across captures.
#' @param cells Explicit selected assay cell IDs, or NULL for all assay cells.
#' @param assay Source ChromatinAssay; NULL selects the active assay.
#' @param macs2.path Path to an installed MACS2 2.x executable in the current
#'   runtime. On Windows use a prepared Linux container with explicit mounted
#'   paths. This function does not manage containers or translate host paths.
#' @param outdir New output directory. Existing paths are rejected. Split BEDs,
#'   native outputs and logs remain here, including when a call fails.
#' @param effective.genome.size Explicit effective genome size passed to MACS2.
#' @param qvalue MACS2 q-value cutoff, strictly between zero and one.
#' @param extsize,shift MACS2 fragment extension size and shift. The supported
#'   path uses BED input and `--nomodel`; MACS3 FRAG arguments are not used.
#' @param verbose Whether Signac reports fragment splitting progress.
#' @details
#' Missing input, nonzero caller exit, missing output and malformed output are
#' errors. A successful zero-byte narrowPeak file is a valid empty peak set;
#' this is distinct from a failed call. No matrix or existing result is changed.
#' Combine the returned sets with [CombineATACPeaks()] and recount the resulting
#' ranges with [RecountATAC()] before comparing previously different universes.
#' Signac's splitting API is reused. The process is invoked directly so its
#' exit status and both output streams can be retained, rather than relying on
#' a caller wrapper that ignores the process exit status.
#' @return A list with named `peaks` (`GRanges` per call), `groups` mapping safe
#'   output IDs to original capture/group labels, and `receipt` containing actual
#'   cell membership, reference, fragment resources, backend versions, command
#'   arguments, logs, exit status and output files. Counts represent the supplied
#'   fragment records, not fifth-column read-support weights. Peak significance
#'   is a caller result, not a biological-replicate condition comparison.
#' @seealso [RunATACQC()], [CombineATACPeaks()], [RecountATAC()]
#' @export
#' @examples
#' \dontrun{
#' called <- CallATACPeaks(atac,
#'   sample.by = "capture", group.by = "celltype",
#'   macs2.path = "/usr/local/bin/macs2", outdir = "atac_peak_calls",
#'   effective.genome.size = 2.7e9
#' )
#' common <- CombineATACPeaks(called$peaks, genome = "hg38")
#' }
CallATACPeaks <- function(object, sample.by, group.by = NULL, cells = NULL,
                          assay = NULL, macs2.path, outdir,
                          effective.genome.size, qvalue = 0.05,
                          extsize = 200L, shift = -extsize / 2, verbose = TRUE) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) log_message("{paste(report$errors, collapse = ' ')}", message_type = "error")
  assay <- report$assay
  if (length(report$genome) != 1L || !nzchar(report$genome)) {
    log_message("Peak calling requires an explicit source genome label.", message_type = "error")
  }
  if (is.null(cells)) cells <- colnames(object[[assay]])
  if (!atac_unique_ids(cells) || !length(cells) || !all(cells %in% colnames(object[[assay]]))) {
    log_message("cells must contain unique selected-assay cell IDs.", message_type = "error")
  }
  metadata <- object[[]][cells, , drop = FALSE]
  fields <- c(sample.by, group.by)
  if (!is.character(sample.by) || length(sample.by) != 1L ||
    (!is.null(group.by) && (!is.character(group.by) || length(group.by) != 1L)) ||
    anyNA(fields) || any(!nzchar(fields)) || anyDuplicated(fields) ||
    !all(fields %in% colnames(metadata))) {
    log_message("Supply distinct sample.by and optional group.by metadata columns.", message_type = "error")
  }
  grouping <- as.data.frame(lapply(metadata[fields], as.character), stringsAsFactors = FALSE)
  if (anyNA(grouping) || any(vapply(grouping, function(x) any(!nzchar(x)), logical(1)))) {
    log_message("Selected cells need complete capture/group labels.", message_type = "error")
  }
  groups <- unique(grouping)
  rownames(groups) <- NULL
  ids <- sprintf("unit%04d", seq_len(nrow(groups)))
  membership <- rep(NA_character_, length(cells))
  for (i in seq_len(nrow(groups))) {
    keep <- Reduce(`&`, lapply(fields, function(field) grouping[[field]] == groups[[field]][i]))
    membership[keep] <- ids[i]
  }
  groups <- data.frame(call_id = ids, groups, check.names = FALSE)
  if (!is.character(macs2.path) || length(macs2.path) != 1L ||
    is.na(macs2.path) || !file.exists(macs2.path)) {
    log_message("macs2.path must identify an installed MACS2 executable.", message_type = "error")
  }
  macs2.path <- normalizePath(macs2.path, mustWork = TRUE)
  version <- system2(macs2.path, "--version", stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(version, "status")) || !any(grepl("macs2[[:space:]]+2[.]", version, ignore.case = TRUE))) {
    log_message("The supported peak caller is MACS2 2.x; inspect macs2.path --version.", message_type = "error")
  }
  if (!is.numeric(effective.genome.size) || length(effective.genome.size) != 1L ||
    !is.finite(effective.genome.size) || effective.genome.size <= 0 ||
    !is.numeric(qvalue) || length(qvalue) != 1L || !is.finite(qvalue) || qvalue <= 0 || qvalue >= 1 ||
    !is.numeric(extsize) || length(extsize) != 1L || !is.finite(extsize) || extsize < 1 || extsize %% 1 != 0 ||
    !is.numeric(shift) || length(shift) != 1L || !is.finite(shift) || shift %% 1 != 0) {
    log_message("Supply a positive genome size/extension, an integer shift and a qvalue between zero and one.", message_type = "error")
  }
  resources <- atac_validate_fragment_inputs(object, assay, cells, verbose)
  if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir) ||
    file.exists(outdir) || nzchar(Sys.readlink(outdir)) || !dir.exists(dirname(outdir))) {
    log_message("outdir must be a new path within an existing directory.", message_type = "error")
  }
  if (!dir.create(outdir)) log_message("Cannot create outdir.", message_type = "error")
  outdir <- normalizePath(outdir, mustWork = TRUE)
  splitdir <- file.path(outdir, "fragments")
  dir.create(splitdir)
  temporary <- subset(object, cells = cells)
  temporary[[assay]]@fragments <- resources$fragments
  temporary[["scop_peak_call_unit"]] <- stats::setNames(membership, cells)[colnames(temporary)]
  Signac::SplitFragments(temporary,
    assay = assay, group.by = "scop_peak_call_unit",
    outdir = splitdir, append = TRUE, verbose = verbose
  )
  peaks <- commands <- stats::setNames(vector("list", length(ids)), ids)
  for (id in ids) {
    bed <- file.path(splitdir, paste0(id, ".bed"))
    if (!file.exists(bed) || file.info(bed)$size == 0) {
      log_message("{paste0('No split fragments for ', id, '; retained output: ', outdir)}", message_type = "error")
    }
    log <- file.path(outdir, paste0(id, ".log"))
    args <- c(
      "callpeak", "-t", bed, "-f", "BED", "--nomodel", "--extsize", extsize,
      "--shift", shift, "-g", effective.genome.size, "-q", qvalue, "-n", id, "--outdir", outdir
    )
    writeLines(c(macs2.path, args), file.path(outdir, paste0(id, ".arguments.txt")))
    status <- system2(macs2.path, shQuote(as.character(args)), stdout = log, stderr = log)
    output <- file.path(outdir, paste0(id, "_peaks.narrowPeak"))
    commands[[id]] <- list(
      executable = macs2.path, arguments = args,
      status = status, log = log, output = output
    )
    saveRDS(commands, file.path(outdir, "commands.rds"))
    if (status != 0L) log_message("{paste0('MACS2 failed for ', id, ' (exit ', status, '). See ', log)}", message_type = "error")
    peaks[[id]] <- atac_read_narrowpeaks(output, report$genome)
  }
  receipt <- list(
    cells = cells, membership = stats::setNames(membership, cells),
    assay = assay, genome = report$genome, fragments = resources$resources,
    commands = commands, outdir = outdir,
    backend = c(Signac = as.character(utils::packageVersion("Signac")), MACS2 = paste(version, collapse = " "))
  )
  result <- list(peaks = peaks, groups = groups, receipt = receipt)
  saveRDS(result, file.path(outdir, "peak_calls.rds"))
  result
}

atac_read_narrowpeaks <- function(file, genome) {
  if (!file.exists(file)) log_message("{paste0('MACS2 did not produce ', file)}", message_type = "error")
  if (file.info(file)$size == 0) {
    return(GenomicRanges::GRanges())
  }
  table <- utils::read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(table) != 10L || anyNA(table) ||
    !is.numeric(table[[2]]) || !is.numeric(table[[3]]) ||
    any(!is.finite(table[[2]])) || any(!is.finite(table[[3]])) ||
    any(table[[2]] < 0 | table[[3]] <= table[[2]]) ||
    any(table[[2]] %% 1 != 0 | table[[3]] %% 1 != 0)) {
    log_message("{paste0('Malformed MACS2 narrowPeak output: ', file)}", message_type = "error")
  }
  names(table) <- c(
    "seqnames", "start", "end", "peak_name", "score", "strand",
    "fold_change", "neg_log10pvalue", "neg_log10qvalue", "summit_offset"
  )
  ranges <- GenomicRanges::makeGRangesFromDataFrame(table,
    keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
  )
  atac_backend_fun("GenomeInfoDb", "genome<-")(ranges, value = genome)
}
