#' @title List cached gene annotation databases
#'
#' @description
#' Lists the gene annotation databases cached by [PrepareDB] and
#' [PrepareCCCDB] in the R.cache root, optionally filtered by species and
#' database name. Each row describes one cached database with its version and
#' creation date.
#'
#' @md
#' @param species Species.
#' Can be `"Homo_sapiens"` or `"Mus_musculus"`.
#' @param db Database names (for example `"GO_BP"`,
#' `"KEGG"`, `"CellChat"`), or a regular expression. If `NULL`, all databases
#' are listed.
#'
#' @return A data frame with columns `Database`, `Species`, `Version`, and
#' `Date`.
#' @export
#'
#' @examples
#' ListDB(species = "Homo_sapiens")
#' ListDB(species = c("Homo_sapiens", "Mus_musculus"))
#' ListDB(species = "Mus_musculus", db = "GO_BP")
ListDB <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = NULL
) {
  dbinfo <- list_db_cache_entries(species = species, db = db)
  if (is.null(dbinfo) || nrow(dbinfo) == 0L) {
    return(data.frame(
      Database = character(),
      Species = character(),
      Version = character(),
      Date = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(
    Database = dbinfo[["DB"]],
    Species = dbinfo[["Species"]],
    Version = dbinfo[["db_version"]],
    Date = dbinfo[["date"]],
    stringsAsFactors = FALSE
  )
  out <- out[order(out[["Species"]], out[["Database"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Full R.cache metadata for the databases (file, timestamp, version, comment)
# used by PrepareDB/PrepareCCCDB cache loading and by ListDB.
list_db_cache_entries <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = NULL
) {
  pathnames <- dir(
    path = R.cache::getCacheRootPath(),
    pattern = "[.]Rcache$",
    full.names = TRUE
  )
  if (length(pathnames) == 0) {
    return(NULL)
  }
  dbinfo <- lapply(
    pathnames, function(x) {
      info <- R.cache::readCacheHeader(x)
      info[["date"]] <- as.character(info[["timestamp"]])
      info[["db_version"]] <- strsplit(info[["comment"]], "\\|")[[1]][1]
      info[["db_name"]] <- strsplit(info[["comment"]], "\\|")[[1]][2]
      info
    }
  )
  dbinfo <- do.call(rbind.data.frame, dbinfo)
  dbinfo[["file"]] <- pathnames

  db_name_parts <- strsplit(as.character(dbinfo[["db_name"]]), "-")
  dbinfo[["Species"]] <- vapply(db_name_parts, function(x) x[1], character(1))
  dbinfo[["DB"]] <- vapply(db_name_parts, function(x) {
    paste(x[-1], collapse = "-")
  }, character(1))

  if (is.null(db)) {
    db <- ".*"
  }
  patterns <- as.vector(outer(species, db, function(s, d) {
    paste0("^", s, "-", d, "$")
  }))
  matched_rows <- unique(unlist(lapply(patterns, function(pat) {
    grep(pat, dbinfo[["db_name"]])
  })))
  if (length(matched_rows) == 0) {
    return(dbinfo[0, , drop = FALSE])
  }
  dbinfo <- dbinfo[matched_rows, , drop = FALSE]
  dbinfo <- dbinfo[
    order(dbinfo[["Species"]], -as.numeric(dbinfo[["timestamp"]])), ,
    drop = FALSE
  ]
  rownames(dbinfo) <- NULL
  dbinfo
}
