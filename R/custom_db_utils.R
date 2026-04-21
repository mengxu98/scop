normalize_custom_db_input <- function(
  TERM2GENE,
  TERM2NAME = NULL,
  IDtype,
  remove_na = FALSE
) {
  colnames(TERM2GENE) <- c("Term", IDtype)
  TERM2GENE <- unique(TERM2GENE)

  if (is.null(TERM2NAME)) {
    TERM2NAME <- TERM2GENE[, c(1, 1), drop = FALSE]
  }
  colnames(TERM2NAME)[seq_len(2)] <- c("Term", "Name")
  TERM2NAME <- unique(TERM2NAME)

  if (isTRUE(remove_na)) {
    TERM2GENE <- stats::na.omit(TERM2GENE)
    TERM2NAME <- stats::na.omit(TERM2NAME)
  }

  list(TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME)
}

create_custom_db_list <- function(
  species,
  db = "custom",
  TERM2GENE,
  TERM2NAME = NULL,
  IDtype,
  version = NULL,
  remove_na = FALSE
) {
  custom_db <- normalize_custom_db_input(
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    IDtype = IDtype,
    remove_na = remove_na
  )

  db_list <- list()
  db_list[[species]] <- list()
  db_list[[species]][[db]] <- custom_db
  if (!is.null(version)) {
    db_list[[species]][[db]][["version"]] <- version
  }

  list(
    db_list = db_list,
    TERM2GENE = custom_db[["TERM2GENE"]],
    TERM2NAME = custom_db[["TERM2NAME"]]
  )
}
