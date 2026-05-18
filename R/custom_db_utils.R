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

normalize_feature_list_input <- function(features) {
  if (!is.list(features) || length(features) == 0L) {
    log_message(
      "{.arg features} must be a non-empty named list",
      message_type = "error"
    )
  }
  feature_names <- names(features)
  if (is.null(feature_names) || any(is.na(feature_names)) || any(!nzchar(trimws(feature_names)))) {
    log_message(
      "{.arg features} must be a named list",
      message_type = "error"
    )
  }

  features <- lapply(features, function(x) {
    genes <- unique(as.character(x))
    genes <- genes[!is.na(genes) & nzchar(trimws(genes))]
    genes
  })
  names(features) <- as.character(feature_names)
  features
}

feature_list_to_term_tables <- function(features, IDtype) {
  features <- normalize_feature_list_input(features)
  term_ids <- paste0("features_", seq_along(features))
  TERM2NAME <- data.frame(
    Term = term_ids,
    Name = names(features),
    stringsAsFactors = FALSE
  )
  TERM2GENE_list <- Map(
    function(term_id, genes) {
      if (length(genes) == 0L) {
        return(NULL)
      }
      out <- data.frame(
        Term = rep(term_id, length(genes)),
        value = genes,
        stringsAsFactors = FALSE
      )
      colnames(out)[2] <- IDtype
      out
    },
    term_ids,
    features
  )
  TERM2GENE_list <- TERM2GENE_list[!vapply(TERM2GENE_list, is.null, logical(1))]
  if (length(TERM2GENE_list) == 0L) {
    TERM2GENE <- data.frame(
      Term = character(0),
      stringsAsFactors = FALSE
    )
    TERM2GENE[[IDtype]] <- character(0)
  } else {
    TERM2GENE <- do.call(rbind, TERM2GENE_list)
  }

  list(
    features = features,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME
  )
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

create_custom_db_list_from_features <- function(
  species,
  features,
  IDtype,
  db = "custom",
  version = NULL,
  remove_na = FALSE
) {
  custom_tables <- feature_list_to_term_tables(
    features = features,
    IDtype = IDtype
  )
  custom_db <- create_custom_db_list(
    species = species,
    db = db,
    TERM2GENE = custom_tables[["TERM2GENE"]],
    TERM2NAME = custom_tables[["TERM2NAME"]],
    IDtype = IDtype,
    version = version,
    remove_na = remove_na
  )
  custom_db[["features"]] <- custom_tables[["features"]]
  custom_db
}
