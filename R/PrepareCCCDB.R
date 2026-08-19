#' @title Prepare cell-cell communication databases
#'
#' @description
#' Prepares the ligand-receptor interaction databases used by the
#' cell-cell communication wrappers: the CellTalkDB ligand-receptor pairs
#' (`"CellTalk"`) and the CellChat curated ligand-receptor database
#' (`"CellChat"`). Each database is returned as a `TERM2GENE`/`TERM2NAME`
#' mapping (with a `"ligand_*"`/`"receptor_*"` term convention) and cached with
#' [R.cache] so that [ListDB] can list it like any other annotation database.
#'
#' @md
#' @inheritParams PrepareDB
#' @param db Cell-cell communication
#' databases to prepare. Can be `"CellTalk"` and/or `"CellChat"`.
#' @param db_update Whether the databases should be forcefully updated.
#' If `FALSE`, cached databases are reused when available.
#'
#' @return A list with the same structure as [PrepareDB]: for each species a
#' named list of databases, each with `TERM2GENE`, `TERM2NAME`, and `version`
#' entries.
#' @export
#'
#' @seealso [PrepareDB], [ListCCCDB], [ListDB]
#'
#' @examples
#' \dontrun{
#' ccc_db <- PrepareCCCDB(
#'   species = "Homo_sapiens",
#'   db = "CellChat"
#' )
#' head(ccc_db[["Homo_sapiens"]][["CellChat"]][["TERM2GENE"]])
#' }
PrepareCCCDB <- function(
  species = c("Homo_sapiens", "Mus_musculus", "Danio_rerio"),
  db = c("CellChat", "CellTalk"),
  convert_species = TRUE,
  data_dir = NULL,
  db_version = "latest",
  db_update = FALSE,
  verbose = TRUE,
  ...
) {
  species <- match.arg(species, several.ok = TRUE)
  db <- match.arg(db, several.ok = TRUE)
  ccc_validate_flag(db_update, "db_update")
  ccc_validate_flag(convert_species, "convert_species")
  db_species <- stats::setNames(rep(species[1], length(db)), nm = db)
  db_list <- stats::setNames(vector("list", length(species)), species)
  for (sps in species) {
    db_list[[sps]] <- list()
    if (isFALSE(db_update)) {
      for (term in db) {
        dbinfo <- list_db_cache_entries(species = sps, db = term)
        if (nrow(dbinfo) > 0 && !is.null(dbinfo)) {
          pathname <- if (identical(db_version, "latest")) {
            dbinfo[order(dbinfo[["timestamp"]], decreasing = TRUE)[1], "file"]
          } else {
            hit <- dbinfo[grep(db_version, dbinfo[["db_version"]], fixed = TRUE)[1], "file"]
            if (is.na(hit)) {
              log_message(
                "There is no {.val {db_version}} version of the database. Use the latest version",
                message_type = "warning",
                verbose = verbose
              )
              dbinfo[order(dbinfo[["timestamp"]], decreasing = TRUE)[1], "file"]
            } else {
              hit
            }
          }
          if (!is.na(pathname)) {
            header <- R.cache::readCacheHeader(pathname)
            cached_version <- strsplit(header[["comment"]], "\\|")[[1]][1]
            timestamp <- format(header[["timestamp"]], "%Y-%m-%d %H:%M:%S")
            log_message(
              "Loading cached: {.pkg {term}} version: {.pkg {cached_version}} created: {.pkg {timestamp}}",
              verbose = verbose
            )
            db_loaded <- R.cache::loadCache(pathname = pathname)
            Sys.sleep(0.5)
            db_list[[sps]][[term]] <- db_loaded
          }
        }
      }
    }
    for (term in db) {
      if (term %in% names(db_list[[sps]])) {
        next
      }
      if (!sps %in% c("Homo_sapiens", "Mus_musculus", "Danio_rerio")) {
        if (isTRUE(convert_species)) {
          log_message(
            "Use the human annotation to create the {.pkg {term}} database for {.val {sps}}",
            message_type = "warning",
            verbose = verbose
          )
          db_species[term] <- "Homo_sapiens"
        } else {
          log_message(
            "{.pkg {term}} database only support {.val {c('Homo_sapiens', 'Mus_musculus', 'Danio_rerio')}}. Consider using {.arg convert_species=TRUE}",
            message_type = "error"
          )
        }
      }
      prepared <- if (identical(term, "CellTalk")) {
        ccc_prepare_celltalk(
          species = db_species[term],
          data_dir = data_dir,
          verbose = verbose
        )
      } else {
        ccc_prepare_cellchat(
          species = db_species[term],
          data_dir = data_dir,
          db_version = db_version,
          verbose = verbose
        )
      }
      db_list[[db_species[term]]][[term]] <- prepared
      if (sps == db_species[term]) {
        R.cache::saveCache(
          prepared,
          key = list(
            prepared$version,
            as.character(db_species[term]),
            term
          ),
          comment = paste0(
            prepared$version,
            " nterm:",
            length(prepared$TERM2NAME[[1]]),
            "|",
            db_species[term],
            "-",
            term
          )
        )
      }
    }
  }
  db_list
}

ccc_validate_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg}} must be a single logical value",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

ccc_prepare_celltalk <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  data_dir = NULL,
  verbose = TRUE
) {
  species <- match.arg(species)
  log_message("Preparing {.pkg CellTalk} database", verbose = verbose)
  url <- switch(species,
    "Homo_sapiens" = "https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/human_lr_pair.rds",
    "Mus_musculus" = "https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/mouse_lr_pair.rds"
  )
  celltalk_file <- switch(species,
    "Homo_sapiens" = "human_lr_pair.rds",
    "Mus_musculus" = "mouse_lr_pair.rds"
  )
  source_file <- preparedb_local_source_file(
    data_dir = data_dir,
    db = "CellTalk",
    pattern = switch(species,
      "Homo_sapiens" = "^human_lr_pair\\.rds$",
      "Mus_musculus" = "^mouse_lr_pair\\.rds$"
    ),
    verbose = verbose
  )
  source_is_temp <- FALSE
  if (is.null(source_file)) {
    source_file <- tempfile()
    source_is_temp <- TRUE
    download(url = url, destfile = source_file)
  }
  lr <- readRDS(source_file)
  if (isTRUE(source_is_temp)) {
    unlink(source_file)
  }
  version <- "v1.0"

  lr[["ligand_gene_symbol2"]] <- paste0(
    "ligand_",
    lr[["ligand_gene_symbol"]]
  )
  lr[["receptor_gene_symbol2"]] <- paste0(
    "receptor_",
    lr[["receptor_gene_symbol"]]
  )
  TERM2GENE <- rbind(
    data.frame(
      "Term" = lr[["ligand_gene_symbol2"]],
      "symbol" = lr[["receptor_gene_symbol"]]
    ),
    data.frame(
      "Term" = lr[["receptor_gene_symbol2"]],
      "symbol" = lr[["ligand_gene_symbol"]]
    )
  )
  TERM2NAME <- TERM2GENE[, c(1, 1)]
  colnames(TERM2GENE) <- c("Term", "symbol")
  colnames(TERM2NAME) <- c("Term", "Name")
  TERM2GENE <- stats::na.omit(unique(TERM2GENE))
  TERM2NAME <- stats::na.omit(unique(TERM2NAME))
  list(
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    version = version
  )
}

ccc_prepare_cellchat <- function(
  species = c("Homo_sapiens", "Mus_musculus", "Danio_rerio"),
  data_dir = NULL,
  db_version = "latest",
  verbose = TRUE
) {
  species <- match.arg(species)
  log_message("Preparing {.pkg CellChat} database", verbose = verbose)
  url <- paste0(
    "https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.",
    switch(species,
      "Homo_sapiens" = "human.rda",
      "Mus_musculus" = "mouse.rda",
      "Danio_rerio" = "zebrafish.rda"
    )
  )
  cellchat_file <- paste0(
    "CellChatDB.",
    switch(species,
      "Homo_sapiens" = "human.rda",
      "Mus_musculus" = "mouse.rda",
      "Danio_rerio" = "zebrafish.rda"
    )
  )
  source_file <- preparedb_local_source_file(
    data_dir = data_dir,
    db = "CellChat",
    pattern = switch(species,
      "Homo_sapiens" = "^CellChatDB\\.human\\.rda$",
      "Mus_musculus" = "^CellChatDB\\.mouse\\.rda$",
      "Danio_rerio" = "^CellChatDB\\.zebrafish\\.rda$"
    ),
    verbose = verbose
  )
  source_is_temp <- FALSE
  if (is.null(source_file)) {
    source_file <- tempfile()
    source_is_temp <- TRUE
    download(url = url, destfile = source_file)
  }
  load(source_file)
  lr <- get(paste0(
    "CellChatDB.",
    switch(species,
      "Homo_sapiens" = "human",
      "Mus_musculus" = "mouse",
      "Danio_rerio" = "zebrafish"
    )
  ))[["interaction"]]
  if (isTRUE(source_is_temp)) {
    temp <- source_file
    download(
      url = "https://raw.githubusercontent.com/sqjin/CellChat/master/DESCRIPTION",
      destfile = temp
    )
    version <- grep(
      pattern = "Version",
      x = readLines(temp),
      value = TRUE
    )
    version <- gsub(
      pattern = "(.*Version: )|(</td>)",
      replacement = "",
      x = version
    )
    unlink(temp)
  } else {
    version <- if (identical(db_version, "latest")) "local" else db_version
  }

  lr_list <- strsplit(lr$interaction_name, split = "_")
  lr[["ligand_gene_symbol"]] <- paste0(
    "ligand_",
    sapply(lr_list, function(x) x[[1]])
  )
  lr[["receptor_list"]] <- lapply(
    lr_list,
    function(x) paste0("receptor_", x[2:length(x)])
  )
  lr <- unnest_fun(data = lr, cols = "receptor_list", keep_empty = FALSE)
  TERM2GENE <- rbind(
    data.frame(
      "Term" = lr[["ligand_gene_symbol"]],
      "symbol" = gsub(
        pattern = "receptor_",
        replacement = "",
        lr[["receptor_list"]]
      )
    ),
    data.frame(
      "Term" = lr[["receptor_list"]],
      "symbol" = gsub(
        pattern = "ligand_",
        replacement = "",
        lr[["ligand_gene_symbol"]]
      )
    )
  )

  if (species == "Homo_sapiens") {
    TERM2GENE[["symbol"]] <- toupper(TERM2GENE[["symbol"]])
  } else if (species == "Mus_musculus") {
    TERM2GENE[["symbol"]] <- capitalize(
      TERM2GENE[["symbol"]],
      force_tolower = TRUE
    )
  } else if (species == "Danio_rerio") {
    TERM2GENE[["symbol"]] <- tolower(TERM2GENE[["symbol"]])
  }
  TERM2NAME <- TERM2GENE[, c(1, 1)]
  colnames(TERM2GENE) <- c("Term", "symbol")
  colnames(TERM2NAME) <- c("Term", "Name")
  TERM2GENE <- stats::na.omit(unique(TERM2GENE))
  TERM2NAME <- stats::na.omit(unique(TERM2NAME))
  list(
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    version = version
  )
}
