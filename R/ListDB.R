#' @title List cached databases
#'
#' @description
#' Retrieves information about databases based on a given species and database name.
#'
#' @md
#' @param species A character vector of species for which to retrieve database information.
#' Default is `c("Homo_sapiens", "Mus_musculus")`.
#' @param db The pattern to match against the database names.
#' Default is `NULL`, which matches all databases.
#'
#' @return A data frame containing information about the databases,
#' including a `Species` column and a `DB` column.
#'
#' @seealso [PrepareDB]
#'
#' @export
#' @examples
#' ListDB(species = "Homo_sapiens")
#' ListDB(species = c("Homo_sapiens", "Mus_musculus"))
#' ListDB(species = "Mus_musculus", db = "GO_BP")
ListDB <- function(
    species = c("Homo_sapiens", "Mus_musculus"),
    db = NULL) {
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
  return(dbinfo)
}
