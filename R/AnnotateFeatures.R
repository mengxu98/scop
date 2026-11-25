#' @title Annotate Features
#'
#' @description
#' Annotate features in a Seurat object with additional metadata from databases or a GTF file.
#'
#' @md
#' @inheritParams PrepareDB
#' @param srt Deprecated alias for `object`; supply exactly one of the two. It
#' will be removed in scop 1.0.0.
#' @param object Seurat object to be annotated.
#' @param IDtype Type of identifier to use for annotation.
#' Options are `"symbol"`, `"ensembl_id"`, or `"entrez_id"`.
#' @param mirror URL of the mirror to use for Ensembl database.
#' @param gtf Path to the GTF file to be used for annotation.
#' @param merge_gtf_by Column name to merge the GTF file by.
#' @param columns Vector of column names to be used from the GTF file.
#' Default is `"seqname"`, `"feature"`, `"start"`, `"end"`,
#' `"strand"`, `"gene_id"`, `"gene_name"`, `"gene_type"`.
#' @param assays Character vector of assay names to be annotated.
#' @param overwrite Whether to overwrite existing metadata.
#'
#' @seealso
#' [PrepareDB], [ListDB]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- AnnotateFeatures(
#'   pancreas_sub,
#'   species = "Mus_musculus",
#'   db = "TF"
#' )
#' head(
#'   GetFeaturesData(
#'     pancreas_sub
#'   )
#' )
#'
#' # Annotate features using a GTF file
#' gtf_file <- "path/to/genes.gtf"
#' \dontrun{
#' pancreas_sub <- AnnotateFeatures(
#'   pancreas_sub,
#'   gtf = gtf_file
#' )
#' head(
#'   GetFeaturesData(
#'     pancreas_sub
#'   )
#' )
#' }
AnnotateFeatures <- function(
  object,
  species = "Homo_sapiens",
  IDtype = c("symbol", "ensembl_id", "entrez_id"),
  db = NULL,
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  gtf = NULL,
  merge_gtf_by = "gene_name",
  columns = c("seqname", "feature", "start", "end", "strand", "gene_id", "gene_name", "gene_type"),
  assays = "RNA",
  overwrite = FALSE,
  verbose = TRUE,
  srt = NULL,
  ...
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  species <- normalize_species_name(species)
  IDtype <- match.arg(IDtype)
  if (is.null(db) && is.null(gtf)) {
    log_message(
      "One of {.arg db} and {.arg gtf} must be specified",
      message_type = "error"
    )
  }

  if (!is.null(db)) {
    db_list <- PrepareDB(
      species = species,
      db = db,
      db_update = db_update,
      db_version = db_version,
      convert_species = convert_species,
      db_IDtypes = IDtype,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      ...
    )
    db_notfound <- setdiff(db, names(db_list[[species]]))
    if (length(db_notfound) > 0) {
      log_message(
        "The following databases are not found: {.val {db_notfound}}",
        message_type = "warning",
        verbose = verbose
      )
    }
    for (single_db in names(db_list[[species]])) {
      TERM2GENE <- unique(db_list[[species]][[single_db]][["TERM2GENE"]])
      TERM2NAME <- unique(db_list[[species]][[single_db]][["TERM2NAME"]])
      rownames(TERM2NAME) <- TERM2NAME[, 1]
      TERM2GENE[, single_db] <- TERM2NAME[TERM2GENE[, 1], 2]
      db_df <- stats::aggregate(
        x = TERM2GENE[,
          !colnames(TERM2GENE) %in%
            c("Term", "entrez_id", "symbol", "ensembl_id"),
          drop = FALSE
        ],
        by = list(rowid = TERM2GENE[[IDtype]]),
        FUN = function(x) {
          paste0(unique(x), collapse = ";")
        }
      )
      rownames(db_df) <- db_df[["rowid"]]
      db_df[["rowid"]] <- NULL
      for (assay in assays) {
        meta_features <- GetFeaturesData(srt, assay = assay)
        if (
          any(colnames(db_df) %in% colnames(meta_features)) && isTRUE(overwrite)
        ) {
          meta_features <- meta_features[, setdiff(
            colnames(meta_features),
            colnames(db_df)
          )]
        }
        db_sub <- db_df[
          rownames(db_df) %in% rownames(meta_features), ,
          drop = FALSE
        ]
        if (nrow(db_sub) == 0) {
          log_message(
            "No data to append was found in {.cls Seurat}. Please check if the species name is correct. The expected feature names are {.val {utils::head(rownames(db_df), 10)}}",
            message_type = "error"
          )
        }
        meta_features <- cbind(
          meta_features,
          db_sub[
            rownames(meta_features),
            setdiff(colnames(db_sub), colnames(meta_features)),
            drop = FALSE
          ]
        )
        srt <- AddFeaturesData(
          srt,
          features = meta_features,
          assay = assay
        )
      }
    }
  }

  if (!is.null(gtf)) {
    gtf_columns_collapse <- gtf_feature_annotations(gtf, columns, merge_gtf_by)
    for (assay in assays) {
      meta_features <- GetFeaturesData(srt, assay = assay)
      if (
        length(intersect(
          colnames(meta_features),
          colnames(gtf_columns_collapse)
        )) >
          0 &&
          isTRUE(overwrite)
      ) {
        meta_features <- meta_features[, setdiff(
          colnames(meta_features),
          colnames(gtf_columns_collapse)
        ), drop = FALSE]
      }
      annotations <- gtf_columns_collapse[
        match(rownames(meta_features), rownames(gtf_columns_collapse)),
        setdiff(colnames(gtf_columns_collapse), colnames(meta_features)),
        drop = FALSE
      ]
      rownames(annotations) <- rownames(meta_features)
      meta_features <- cbind(meta_features, annotations)
      srt <- AddFeaturesData(
        srt,
        features = meta_features,
        assay = assay
      )
    }
  }
  return(srt)
}

gtf_feature_annotations <- function(gtf, columns, merge_gtf_by) {
  if (!is.character(columns) || !length(columns) || anyNA(columns) || any(!nzchar(columns))) {
    log_message("{.arg columns} must contain GTF column or attribute names.", message_type = "error")
  }
  if (!is.character(merge_gtf_by) || length(merge_gtf_by) != 1L ||
    is.na(merge_gtf_by) || !nzchar(merge_gtf_by)) {
    log_message("{.arg merge_gtf_by} must name one GTF column or attribute.", message_type = "error")
  }
  columns <- unique(columns)
  records <- read_gtf_feature_records(gtf)
  attrs <- setdiff(c(columns, merge_gtf_by), names(records))
  if (length(attrs)) {
    # Match semicolon-delimited keys without splitting quoted values.
    pattern <- r"{(?:^|;)[[:space:]]*([^;[:space:]]+)[[:space:]]+("(?:[^"\\]|\\.)*"|[^";]*?)[[:space:]]*(?=;|$)}"
    attributes <- regmatches(records$attribute, gregexec(pattern, records$attribute, perl = TRUE))
    for (j in seq_along(attrs)) {
      records[[attrs[j]]] <- vapply(attributes, function(fields) {
        if (!length(fields)) {
          return(NA_character_)
        }
        value <- fields[3L, match(attrs[j], fields[2L, ])]
        sub('^"(.*)"$', '\\1', value)
      }, character(1))
    }
  }
  keys <- records[[merge_gtf_by]]
  found <- !is.na(keys) & nzchar(keys) & keys != "."
  if (!any(found)) {
    log_message(
      "No GTF records contain values for {.arg merge_gtf_by} ({.val {merge_gtf_by}}).",
      message_type = "error"
    )
  }
  out <- stats::aggregate(
    records[found, columns, drop = FALSE],
    by = list(rowid = keys[found]),
    FUN = function(x) {
      values <- unique(x[!is.na(x)])
      if (length(values)) paste(values, collapse = ";") else NA_character_
    }
  )
  # aggregate() sorts grouping values; keep the actual key attached to each row.
  rownames(out) <- out[[1L]]
  out[, -1L, drop = FALSE]
}

read_gtf_feature_records <- function(gtf) {
  fields <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  precedence <- c("gene", "transcript", "exon", "CDS")
  retained <- list()
  best <- length(precedence) + 1L
  con <- file(gtf, open = "rt")
  on.exit(close(con))
  repeat {
    lines <- readLines(con, n = 4096L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!grepl("^[[:space:]]*(#|$)", lines)]
    if (!length(lines)) next
    parts <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(parts) != length(fields))) {
      log_message("GTF records must contain exactly nine tab-separated fields.", message_type = "error")
    }
    block <- do.call(rbind, parts)
    rank <- match(block[, 3L], precedence, nomatch = length(precedence) + 1L)
    if (min(rank) < best) {
      retained <- list()
      best <- min(rank)
    }
    if (any(rank == best)) {
      retained[[length(retained) + 1L]] <- block[rank == best, , drop = FALSE]
    }
  }
  if (!length(retained)) {
    log_message("No annotation records were found in the GTF file.", message_type = "error")
  }
  result <- as.data.frame(do.call(rbind, retained), stringsAsFactors = FALSE)
  names(result) <- fields
  result
}
