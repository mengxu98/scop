#' Count strand-oriented insertion endpoints around ATAC regions
#'
#' @md
#' @description Build group-specific region-by-position matrices from local
#' indexed fragments, without subsampling cells or regions. Repeated endpoints
#' are counted with their multiplicity. This is an explicit SCOP endpoint
#' counter, not a call to Signac's RegionMatrix.
#' @param object A Seurat object containing a validated ChromatinAssay with
#'   local indexed fragments and explicit cell-to-barcode mappings.
#' @param regions Non-empty GRanges in 1-based closed coordinates. Each region
#'   is centered at start + floor((width - 1)/2). Minus-strand regions are
#'   reversed; unstranded regions use genomic orientation. Duplicate genomic
#'   intervals are rejected even when their names differ. Unique names are
#'   used as row IDs when provided, otherwise genomic coordinate IDs are used.
#' @param assay Chromatin assay; defaults to the active assay.
#' @param group.by Metadata column for grouping. NULL uses active identities.
#' @param cells Unique source-assay cells; NULL includes all assay cells.
#' @param upstream,downstream Non-negative integer distances from the center.
#'   Windows must fit entirely within known chromosome boundaries; no trimming
#'   or reference/chromosome conversion is performed.
#' @param bin.size Positive integer width of positional bins. The last bin may
#'   be shorter; actual widths and centers are returned.
#' @param max.entries Maximum number of dense output entries across all groups.
#'   Exceeding this explicit memory guard errors rather than sampling. This
#'   guards output size, not transient memory for a single dense fragment query.
#' @param verbose Whether to print progress messages.
#' @details A BED fragment [start,end) contributes two endpoints, start + 1
#'   and end, in 1-based coordinates. Only endpoints inside each region window
#'   are counted. Read-support column 5 is not a weight. There is no extra Tn5
#'   shift, bias correction, fragment-length weighting or depth normalization.
#'   Overlapping windows can count the same endpoint. Grouping is descriptive,
#'   not a replicate-level statistical test; library-size and cell-composition
#'   confounding remain. Declared reference contigs absent from a file produce
#'   zero counts for that file, not evidence of biological inaccessibility.
#'   The caller's object is not modified.
#' @return A scop_atac_regions list with matrices (one raw count matrix per
#'   group), regions, bin positions/widths, cells.per.group, and a receipt.
#' @seealso [ATACRegionHeatmap()], [ATACRegionPlot()]
#' @export
RunATACRegionMatrix <- function(
  object, regions, assay = NULL, group.by = NULL, cells = NULL,
  upstream = 3000L, downstream = 3000L, bin.size = 50L,
  max.entries = 5e7, verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) stop(paste(report$errors, collapse = " "), call. = FALSE)
  assay <- report$assay
  chromatin <- object[[assay]]
  if (length(report$genome) != 1L) stop("Region counting requires one verified assay genome label.", call. = FALSE)
  problems <- atac_ranges_problems(regions)
  if (length(problems)) stop(paste(problems, collapse = " "), call. = FALSE)
  if (!length(regions)) stop("regions must not be empty.", call. = FALSE)
  atac_qc_ranges(regions, chromatin, "Regions")
  for (x in list(upstream, downstream)) atac_integer(x, "upstream/downstream", 0)
  atac_integer(bin.size, "bin.size", 1)
  if (!is.numeric(max.entries) || length(max.entries) != 1L ||
    !is.finite(max.entries) || max.entries < 1) {
    stop("max.entries must be positive.", call. = FALSE)
  }
  if (as.double(upstream) + downstream + 1 > .Machine$integer.max) stop("Region window is too wide.", call. = FALSE)
  ids <- names(regions)
  if (is.null(ids)) ids <- atac_range_ids(regions)
  if (!atac_unique_ids(ids)) stop("Regions require unique non-empty IDs or coordinates.", call. = FALSE)
  if (is.null(cells)) cells <- colnames(chromatin)
  if (!atac_unique_ids(cells) || !all(cells %in% colnames(chromatin))) {
    stop("cells must be unique non-empty source-assay cell IDs.", call. = FALSE)
  }
  if (is.null(group.by)) {
    groups <- as.character(SeuratObject::Idents(object)[cells])
  } else {
    if (!is.character(group.by) || length(group.by) != 1L || is.na(group.by) ||
      !group.by %in% colnames(object[[]])) {
      stop("group.by must identify one metadata column.", call. = FALSE)
    }
    groups <- as.character(object[[]][cells, group.by])
  }
  if (anyNA(groups) || any(!nzchar(groups))) stop("Selected cells must have non-missing, non-empty groups.", call. = FALSE)
  names(groups) <- cells
  levels <- unique(groups)
  window.width <- as.double(upstream) + downstream + 1
  nbins <- ceiling(window.width / bin.size)
  if (length(regions) * nbins * length(levels) > max.entries) {
    stop("Output exceeds max.entries; use wider bins, explicit fewer regions, or raise the guard knowingly.", call. = FALSE)
  }
  first <- seq.int(0, window.width - 1, by = bin.size)
  widths <- pmin(bin.size, window.width - first)
  positions <- -upstream + first + (widths - 1) / 2
  centers <- GenomicRanges::start(regions) + floor((GenomicRanges::width(regions) - 1) / 2)
  minus <- as.character(GenomicRanges::strand(regions)) == "-"
  window.start <- centers - ifelse(minus, downstream, upstream)
  window.end <- centers + ifelse(minus, upstream, downstream)
  chromosomes <- as.character(GenomicRanges::seqnames(regions))
  source.lengths <- GenomeInfoDb::seqlengths(chromatin)
  supplied.lengths <- GenomeInfoDb::seqlengths(regions)
  shared <- source.lengths[names(supplied.lengths)]
  if (any(shared != supplied.lengths, na.rm = TRUE)) stop("Region and assay chromosome lengths conflict.", call. = FALSE)
  limits <- source.lengths[chromosomes]
  missing <- is.na(limits)
  limits[missing] <- supplied.lengths[chromosomes[missing]]
  if (any(window.start < 1) || any(window.end > limits, na.rm = TRUE)) {
    stop("Extended region windows cross chromosome boundaries; reduce the flanks or explicitly select valid regions.", call. = FALSE)
  }
  windows <- GenomicRanges::GRanges(chromosomes, IRanges::IRanges(window.start, window.end))
  fragments <- atac_indexed_fragments(chromatin, cells, verbose)
  known <- union(names(source.lengths), unlist(lapply(fragments, `[[`, "contigs"), use.names = FALSE))
  if (length(setdiff(chromosomes, known))) stop("Regions contain unknown reference contigs.", call. = FALSE)
  matrices <- lapply(levels, function(x) {
    matrix(0,
      nrow = length(regions), ncol = nbins,
      dimnames = list(ids, as.character(positions))
    )
  })
  names(matrices) <- levels
  for (file in fragments) {
    connection <- Rsamtools::TabixFile(file$path, index = file$index)
    open(connection)
    # Each connection is closed even when parsing or counting fails.
    matrices <- tryCatch(
      {
        for (i in which(chromosomes %in% file$contigs)) {
          records <- Rsamtools::scanTabix(connection, param = windows[i])[[1L]]
          if (!length(records)) next
          parts <- strsplit(records, "\t", fixed = TRUE)
          if (any(lengths(parts) < 4L)) stop("Malformed fragment record: at least four fields required.", call. = FALSE)
          barcodes <- vapply(parts, `[[`, character(1), 4L)
          selected <- match(barcodes, unname(file$map))
          keep <- !is.na(selected)
          if (!any(keep)) next
          parts <- parts[keep]
          mapped.groups <- groups[names(file$map)[selected[keep]]]
          starts <- suppressWarnings(as.numeric(vapply(parts, `[[`, character(1), 2L)))
          ends <- suppressWarnings(as.numeric(vapply(parts, `[[`, character(1), 3L)))
          if (any(!is.finite(starts) | !is.finite(ends) | starts < 0 | ends <= starts |
            starts %% 1 != 0 | ends %% 1 != 0 | ends > .Machine$integer.max)) {
            stop("Invalid BED fragment coordinates.", call. = FALSE)
          }
          for (group in unique(mapped.groups)) {
            index <- mapped.groups == group
            endpoints <- c(starts[index] + 1, ends[index])
            offsets <- if (minus[i]) centers[i] - endpoints else endpoints - centers[i]
            offsets <- offsets[offsets >= -upstream & offsets <= downstream]
            bins <- floor((offsets + upstream) / bin.size) + 1L
            matrices[[group]][i, ] <- matrices[[group]][i, ] + tabulate(bins, nbins = nbins)
          }
        }
        matrices
      },
      finally = close(connection)
    )
  }
  names(regions) <- ids
  result <- list(
    matrices = matrices, regions = regions, positions = positions,
    bin.widths = widths, cells.per.group = table(factor(groups, levels = levels)),
    receipt = list(
      source_assay = assay, cells = cells, groups = groups,
      upstream = upstream, downstream = downstream, bin.size = bin.size,
      genome = report$genome, coordinate_system = "1-based closed endpoints: BED start+1 and end",
      endpoint_shift = 0L, column5_weighted = FALSE, bias_corrected = FALSE,
      fragments = lapply(fragments, function(x) x[c("path", "index", "hash")]),
      backend = "SCOP indexed endpoint tabulation", signac_version = as.character(utils::packageVersion("Signac"))
    )
  )
  class(result) <- c("scop_atac_regions", "list")
  result
}

#' Plot positional ATAC region signal
#'
#' @md
#' @param x Output of [RunATACRegionMatrix()].
#' @param normalize Divide counts by the number of selected cells in each group.
#' @param per.base Divide counts by actual bin widths (including a shorter last
#'   bin). TRUE makes differently sized bins comparable.
#' @param groups Optional unique group names to plot, in display order.
#' @param log1p Apply log1p to displayed heatmap values only.
#' @param order.rows Order heatmap rows by their total displayed signal across
#'   all selected groups; the same order is used in every panel.
#' @param low,high Heatmap color endpoints.
#' @return A ggplot with its summarized plotting data. The line plot shows the
#'   mean across supplied regions, not a confidence interval or enrichment test.
#' @seealso [RunATACRegionMatrix()]
#' @export
ATACRegionHeatmap <- function(
  x, normalize = TRUE, per.base = TRUE,
  groups = NULL, log1p = FALSE, order.rows = FALSE, low = "white", high = "#2166AC"
) {
  data <- atac_region_plot_data(x, normalize, per.base, groups)
  for (value in list(log1p, order.rows)) if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("Plot flags must be TRUE or FALSE.", call. = FALSE)
  ids <- rownames(x$matrices[[1L]])
  if (order.rows) {
    totals <- tapply(data$value, data$region, sum)
    ids <- ids[order(totals[ids])]
  }
  data$region <- factor(data$region, levels = ids)
  if (log1p) data$value <- base::log1p(data$value)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$position, y = .data$region, fill = .data$value)) +
    ggplot2::geom_tile(ggplot2::aes(width = .data$bin.width), height = 1) +
    ggplot2::facet_wrap(~group) +
    ggplot2::scale_fill_gradient(low = low, high = high) +
    ggplot2::labs(
      x = "Distance from region center (bp)", y = "Region",
      fill = paste0(if (log1p) "log1p " else "", atac_region_units(normalize, per.base))
    ) +
    ggplot2::theme_classic()
}

#' @rdname ATACRegionHeatmap
#' @export
ATACRegionPlot <- function(x, normalize = TRUE, per.base = TRUE, groups = NULL) {
  data <- atac_region_plot_data(x, normalize, per.base, groups)
  data <- stats::aggregate(value ~ position + group, data = data, FUN = mean)
  ggplot2::ggplot(data, ggplot2::aes(x = .data$position, y = .data$value, color = .data$group, group = .data$group)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Distance from region center (bp)",
      y = paste("Mean across regions:", atac_region_units(normalize, per.base)), color = "Group"
    ) +
    ggplot2::theme_classic()
}

atac_region_units <- function(normalize, per.base) {
  paste0("Endpoints", if (normalize) "/cell" else "", if (per.base) "/bp" else "/bin")
}

atac_region_plot_data <- function(x, normalize, per.base, groups) {
  if (!inherits(x, "scop_atac_regions") || !is.list(x$matrices) || !length(x$matrices)) stop("x must be a RunATACRegionMatrix result.", call. = FALSE)
  for (flag in list(normalize, per.base)) if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("normalize and per.base must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(x$positions) || !length(x$positions) || any(!is.finite(x$positions)) ||
    any(diff(x$positions) <= 0) || !is.numeric(x$bin.widths)) {
    stop("Malformed region positions or bin widths.", call. = FALSE)
  }
  reference.ids <- rownames(x$matrices[[1L]])
  if (is.null(groups)) groups <- names(x$matrices)
  if (!atac_unique_ids(groups) || !all(groups %in% names(x$matrices))) stop("groups must be unique existing group names.", call. = FALSE)
  result <- lapply(groups, function(group) {
    m <- x$matrices[[group]]
    if (!is.matrix(m) || !is.numeric(m) || any(!is.finite(m)) || any(m < 0) ||
      ncol(m) != length(x$positions) || ncol(m) != length(x$bin.widths) ||
      !atac_unique_ids(rownames(m)) || !identical(rownames(m), reference.ids) || any(!is.finite(x$bin.widths) | x$bin.widths <= 0)) {
      stop("Malformed region matrix or bin metadata.", call. = FALSE)
    }
    if (normalize) {
      n <- x$cells.per.group[[group]]
      if (length(n) != 1L || !is.finite(n) || n <= 0) stop("Invalid group cell denominator.", call. = FALSE)
      m <- m / n
    }
    if (per.base) m <- sweep(m, 2L, x$bin.widths, "/")
    data.frame(
      region = rep(rownames(m), times = ncol(m)),
      position = rep(x$positions, each = nrow(m)),
      bin.width = rep(x$bin.widths, each = nrow(m)), value = as.vector(m),
      group = group, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, result)
  result$group <- factor(result$group, levels = groups)
  result
}
