#' Count fragments in fixed-width genome tiles
#'
#' @md
#' @description Construct explicit non-overlapping genome tiles and delegate
#'   validated sparse fragment counting to [RecountATAC()]. No peaks are needed
#'   for defining the tile universe; the source ChromatinAssay supplies verified
#'   fragment mappings and reference provenance.
#' @param object A Seurat object with a validated ChromatinAssay and fragments.
#' @param assay Source chromatin assay; NULL uses the active assay.
#' @param seqlengths Optional named positive integer chromosome lengths. NULL
#'   uses declared assay lengths, which must all be known for selected contigs.
#'   Supplied lengths must agree with known assay lengths. Names are not used
#'   to guess assembly or chromosome style.
#' @param chromosomes Optional unique contig names in desired output order.
#' @param bin.size Positive integer tile width. The last tile per chromosome is
#'   truncated at that chromosome's end, never dropped or extended.
#' @param cells Selected source-assay cells; see [RecountATAC()].
#' @param new.assay Optional unused assay name; NULL returns counts/ranges/receipt.
#' @param process_n FeatureMatrix chunk size; see [RecountATAC()].
#' @param max.tiles Positive integer upper bound on tile count before allocation.
#' @param verbose Whether to print progress messages.
#' @return The [RecountATAC()] return value with tile parameters in its receipt.
#'   Counts are overlapping fragment records, not insertion endpoints or unique
#'   fragments across bins. A fragment crossing a boundary can count in two bins.
#' @seealso [RecountATAC()]
#' @export
RunGenomeTiles <- function(object, assay = NULL, seqlengths = NULL,
  chromosomes = NULL, bin.size = 5000L, cells = NULL, new.assay = NULL,
  process_n = 2000L, max.tiles = 2000000L, verbose = TRUE
) {
  report <- CheckATACData(object, assay = assay)
  if (!report$valid) stop(paste(report$errors, collapse = " "), call. = FALSE)
  assay <- report$assay
  atac_tiles_integer(bin.size, "bin.size", 1)
  atac_tiles_integer(max.tiles, "max.tiles", 1)
  source.lengths <- atac_backend_fun("GenomeInfoDb", "seqlengths")(object[[assay]])
  if (is.null(seqlengths)) seqlengths <- source.lengths
  if (!is.numeric(seqlengths) || !atac_unique_ids(names(seqlengths))) stop("seqlengths must be a named numeric vector with unique contigs.", call. = FALSE)
  if (!is.null(chromosomes)) {
    if (!atac_unique_ids(chromosomes) || !all(chromosomes %in% names(seqlengths))) stop("chromosomes must select unique known contigs.", call. = FALSE)
    seqlengths <- seqlengths[chromosomes]
  }
  if (any(!is.finite(seqlengths) | seqlengths < 1 | seqlengths %% 1 != 0 |
      seqlengths > .Machine$integer.max)) stop("Selected chromosome lengths must be known positive R-range integers.", call. = FALSE)
  if (any(source.lengths[names(seqlengths)] != seqlengths, na.rm = TRUE)) stop("Supplied chromosome lengths conflict with assay reference.", call. = FALSE)
  ntiles <- sum(ceiling(seqlengths / bin.size))
  if (ntiles > max.tiles) stop("Tile count exceeds max.tiles; increase bin.size or explicitly select contigs.", call. = FALSE)
  ranges <- GenomicRanges::tileGenome(seqlengths, tilewidth = bin.size, cut.last.tile.in.chrom = TRUE)
  if (length(report$genome) != 1L) stop("Genome tiles require one verified assay genome label.", call. = FALSE)
  GenomeInfoDb::genome(ranges) <- report$genome
  result <- RecountATAC(object, ranges, assay = assay, cells = cells,
    new.assay = new.assay, process_n = process_n, verbose = verbose)
  tile.receipt <- list(bin.size = bin.size, seqlengths = seqlengths,
    chromosomes = names(seqlengths), n_tiles = length(ranges),
    last_tile = "truncated at chromosome end")
  if (is.null(new.assay)) {
    result$receipt$tiles <- tile.receipt
  } else {
    receipt <- SeuratObject::Misc(result[[new.assay]], "scop_atac_recount")
    receipt$tiles <- tile.receipt
    # Enrich only the receipt of the newly created assay. Avoid Misc's
    # overwrite warning for this expected addition of tile provenance.
    result[[new.assay]]@misc$scop_atac_recount <- receipt
  }
  result
}

atac_tiles_integer <- function(x, label, minimum = 0) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
      x %% 1 != 0 || x > .Machine$integer.max) stop(paste(label, "must be an integer >=", minimum), call. = FALSE)
  invisible(x)
}
