#' @title Prepare KEGG databases
#' @description Prepare KEGG annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"KEGG"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareKEGG(species = "Homo_sapiens", db = "KEGG")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareKEGG <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "KEGG",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "KEGG")) {
    log_message("Unsupported {.arg db} selector for PrepareKEGG", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "KEGG") && (!"KEGG" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg KEGG} database", verbose = verbose)
        check_r("httr", verbose = FALSE)
        orgs <- kegg_get("https://rest.kegg.jp/list/genome")
        orgs_parsed <- strsplit(orgs[, 2], "; ")
        orgs_code <- vapply(orgs_parsed, `[`, "", 1)
        orgs_species <- vapply(orgs_parsed, `[`, "", 2)
        kegg_sp <- orgs_code[
          grep(gsub(pattern = "_", replacement = " ", x = sps), orgs_species)
        ]
        if (length(kegg_sp) == 0) {
          db_species_name <- db_species["KEGG"]
          log_message(
            "Failed to prepare the KEGG database for {.val {db_species_name}}",
            message_type = "warning",
            verbose = verbose
          )
          if (isTRUE(convert_species) && db_species_name != "Homo_sapiens") {
            log_message(
              "Use the human annotation to create the KEGG database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["KEGG"] <- "Homo_sapiens"
            kegg_sp <- "hsa"
            return(NULL)
          } else {
            log_message(
              "Stop the preparation",
              message_type = "error"
            )
          }
        }
        kegg_db <- "pathway"

        kegg_pathwaygene_url <- paste0(
          "https://rest.kegg.jp/link/",
          kegg_sp,
          "/",
          kegg_db,
          collapse = ""
        )
        TERM2GENE <- kegg_get(kegg_pathwaygene_url)
        colnames(TERM2GENE) <- c("Pathway", "KEGG_ID")
        kegg_geneconversion_url <- paste0(
          "https://rest.kegg.jp/conv/ncbi-geneid/",
          kegg_sp
        )
        GENECONV <- kegg_get(kegg_geneconversion_url)
        colnames(GENECONV) <- c("KEGG_ID", "ENTREZID")
        TERM2GENE <- merge(
          x = TERM2GENE,
          y = GENECONV,
          by = "KEGG_ID",
          all.x = TRUE
        )
        TERM2GENE[, "Pathway"] <- gsub(
          pattern = "[^:]+:",
          replacement = "",
          x = TERM2GENE[, "Pathway"]
        )
        TERM2GENE[, "ENTREZID"] <- gsub(
          pattern = "[^:]+:",
          replacement = "",
          x = TERM2GENE[, "ENTREZID"]
        )
        TERM2GENE <- TERM2GENE[, c("Pathway", "ENTREZID")]

        kegg_pathwayname_url <- paste0(
          "https://rest.kegg.jp/list/",
          kegg_db,
          "/",
          kegg_sp,
          collapse = ""
        )
        TERM2NAME <- kegg_get(kegg_pathwayname_url)
        colnames(TERM2NAME) <- c("Pathway", "Name")
        TERM2NAME[, "Pathway"] <- gsub(
          pattern = "[^:]+:",
          replacement = "",
          x = TERM2NAME[, "Pathway"]
        )
        TERM2NAME[, "Name"] <- gsub(
          pattern = paste0(
            " - ",
            paste0(
              unlist(strsplit(db_species["KEGG"], split = "_")),
              collapse = " "
            ),
            ".*$"
          ),
          replacement = "",
          x = TERM2NAME[, "Name"]
        )
        TERM2NAME <- TERM2NAME[
          TERM2NAME[, "Pathway"] %in% TERM2GENE[, "Pathway"], ,
          drop = FALSE
        ]

        colnames(TERM2GENE) <- c("Term", default_id_types[["KEGG"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        kegg_info <- strsplit(
          httr::content(httr::GET(paste0(
            "https://rest.kegg.jp/info/",
            kegg_sp
          ))),
          split = "\n"
        )[[1]]
        version <- kegg_release_version(kegg_info)
        db_list[[db_species["KEGG"]]][["KEGG"]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species["KEGG"]]][["KEGG"]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species["KEGG"]]][["KEGG"]][["version"]] <- version
        if (sps == db_species["KEGG"]) {
          preparedb_cache_annotation(
            db_list[[db_species["KEGG"]]][["KEGG"]],
            species = as.character(db_species["KEGG"]), db = "KEGG"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare WikiPathway databases
#' @description Prepare WikiPathway annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"WikiPathway"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareWikiPathway(species = "Homo_sapiens", db = "WikiPathway")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareWikiPathway <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "WikiPathway",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "WikiPathway")) {
    log_message("Unsupported {.arg db} selector for PrepareWikiPathway", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (
        any(db == "WikiPathway") &&
          (!"WikiPathway" %in% names(db_list[[sps]]))
      ) {
        log_message("Preparing {.pkg WikiPathway} database", verbose = verbose)
        tempdir <- tempdir()
        gmt_files <- list.files(tempdir)[grep(
          ".gmt",
          x = list.files(tempdir)
        )]
        if (length(gmt_files) > 0) {
          file.remove(paste0(tempdir, "/", gmt_files))
        }
        temp <- tempfile()
        wiki_source_url <- if (is.null(mirror)) {
          "https://data.wikipathways.org/current/gmt"
        } else {
          mirror
        }
        wiki_source_url <- sub("/+$", "", wiki_source_url)
        wiki_file_url <- NULL
        download(
          quiet = TRUE, url = wiki_source_url,
          destfile = temp
        )
        lines <- paste0(readLines(temp, warn = FALSE), collapse = " ")
        gmtfiles <- unlist(regmatches(
          lines,
          m = gregexpr(
            "wikipathways-[^\"'<>[:space:]]+\\.gmt\\b",
            lines,
            perl = TRUE
          )
        ))
        gmtfiles <- unique(gmtfiles)
        if (
          length(gmtfiles) == 0 &&
            identical(
              wiki_source_url,
              "https://wikipathways-data.wmcloud.org/current/gmt"
            )
        ) {
          wiki_source_url <- "https://data.wikipathways.org/current/gmt"
          download(
            quiet = TRUE, url = wiki_source_url,
            destfile = temp
          )
          lines <- paste0(readLines(temp, warn = FALSE), collapse = " ")
          gmtfiles <- unlist(regmatches(
            lines,
            m = gregexpr(
              "wikipathways-[^\"'<>[:space:]]+\\.gmt\\b",
              lines,
              perl = TRUE
            )
          ))
          gmtfiles <- unique(gmtfiles)
        }
        if (
          length(gmtfiles) == 0 &&
            grepl("\\.gmt([?#].*)?$", wiki_source_url, ignore.case = TRUE)
        ) {
          wiki_file_url <- sub("[?#].*$", "", wiki_source_url)
          gmtfiles <- basename(wiki_file_url)
        }
        wiki_sp <- sps
        gmtfile <- gmtfiles[grep(wiki_sp, gmtfiles, fixed = TRUE)]
        if (length(gmtfile) == 0) {
          db_species_name <- db_species["WikiPathway"]
          log_message(
            "Failed to prepare the WikiPathway database for {.val {db_species_name}}",
            message_type = "warning",
            verbose = verbose
          )
          if (isTRUE(convert_species) && db_species_name != "Homo_sapiens") {
            log_message(
              "Use the human annotation to create the WikiPathway database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["WikiPathway"] <- "Homo_sapiens"
            wiki_sp <- "Homo_sapiens"
            gmtfile <- gmtfiles[grep(wiki_sp, gmtfiles, fixed = TRUE)]
          } else {
            log_message(
              "Stop the preparation",
              message_type = "error"
            )
          }
        }
        if (length(gmtfile) == 0) {
          log_message(
            c(
              "No {.pkg WikiPathway} GMT file is available for {.val {wiki_sp}}",
              "Check whether {.arg mirror} points to a GMT file or a directory index"
            ),
            message_type = "error"
          )
        }
        gmtfile <- gmtfile[[1]]
        version_parts <- strsplit(gmtfile, split = "-", fixed = TRUE)[[1]]
        version <- if (length(version_parts) >= 2) {
          version_parts[[2]]
        } else {
          tools::file_path_sans_ext(gmtfile)
        }
        if (is.null(wiki_file_url)) {
          wiki_file_url <- paste0(wiki_source_url, "/", gmtfile)
        }
        download(
          quiet = TRUE, url = wiki_file_url,
          destfile = temp
        )
        wiki_gmt <- clusterProfiler::read.gmt(temp)
        unlink(temp)
        wiki_gmt <- apply(wiki_gmt, 1, function(x) {
          wikiid <- strsplit(x[["term"]], split = "%")[[1]][3]
          wikiterm <- strsplit(x[["term"]], split = "%")[[1]][1]
          gmt <- x[["gene"]]
          data.frame(
            v0 = wikiid,
            v1 = gmt,
            v2 = wikiterm,
            stringsAsFactors = FALSE
          )
        })
        bg <- do.call(rbind.data.frame, wiki_gmt)
        TERM2GENE <- bg[, c(1, 2)]
        TERM2NAME <- bg[, c(1, 3)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["WikiPathway"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        db_list[[db_species["WikiPathway"]]][["WikiPathway"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["WikiPathway"]]][["WikiPathway"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["WikiPathway"]]][["WikiPathway"]][[
          "version"
        ]] <- version
        if (sps == db_species["WikiPathway"]) {
          preparedb_cache_annotation(
            db_list[[db_species["WikiPathway"]]][["WikiPathway"]],
            species = as.character(db_species["WikiPathway"]), db = "WikiPathway"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare Reactome databases
#' @description Prepare Reactome annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character selector; must be `"Reactome"`.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareReactome(species = "Homo_sapiens", db = "Reactome")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareReactome <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "Reactome",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(db %in% "Reactome")) {
    log_message("Unsupported {.arg db} selector for PrepareReactome", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      if (any(db == "Reactome") && (!"Reactome" %in% names(db_list[[sps]]))) {
        log_message("Preparing {.pkg Reactome} database", verbose = verbose)
        reactome_sp <- gsub(pattern = "_", replacement = " ", x = sps)
        reactome_db <- get_namespace_fun("reactome.db", "reactome.db")
        df_all <- suppressMessages(
          AnnotationDbi::select(
            reactome_db,
            keys = AnnotationDbi::keys(reactome_db),
            columns = c("PATHID", "PATHNAME")
          )
        )
        df <- df_all[
          grepl(
            pattern = paste0("^", reactome_sp, ": "),
            x = df_all$PATHNAME
          ), ,
          drop = FALSE
        ]
        if (nrow(df) == 0) {
          if (
            isTRUE(convert_species) &&
              db_species["Reactome"] != "Homo_sapiens"
          ) {
            log_message(
              "Use the human annotation to create the Reactome database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["Reactome"] <- "Homo_sapiens"
            reactome_sp <- gsub(
              pattern = "_",
              replacement = " ",
              x = "Homo_sapiens"
            )
            df <- df_all[
              grepl(
                pattern = paste0("^", reactome_sp, ": "),
                x = df_all$PATHNAME
              ), ,
              drop = FALSE
            ]
          } else {
            log_message(
              "Stop the preparation",
              message_type = "error"
            )
          }
        }
        df <- stats::na.omit(df)
        df$PATHNAME <- gsub(
          x = df$PATHNAME,
          pattern = paste0("^", reactome_sp, ": "),
          replacement = "",
          perl = TRUE
        )
        TERM2GENE <- df[, c(2, 1)]
        TERM2NAME <- df[, c(2, 3)]
        colnames(TERM2GENE) <- c("Term", default_id_types[["Reactome"]])
        colnames(TERM2NAME) <- c("Term", "Name")
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        version <- utils::packageVersion("reactome.db")
        db_list[[db_species["Reactome"]]][["Reactome"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["Reactome"]]][["Reactome"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["Reactome"]]][["Reactome"]][[
          "version"
        ]] <- version
        if (sps == db_species["Reactome"]) {
          preparedb_cache_annotation(
            db_list[[db_species["Reactome"]]][["Reactome"]],
            species = as.character(db_species["Reactome"]), db = "Reactome"
          )
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}

#' @title Prepare MSigDB databases
#' @description Prepare MSigDB annotations, reusing the common cache and
#' species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @param db Character vector selecting `"MSigDB"` or `"MSigDB_<collection>"`.
#' @details Collections can be selected with `MSigDB_<collection>`, including
#' nested prefixes such as `MSigDB_M2_CP_BIOCARTA`. Colon forms are accepted.
#' @return A named species list containing the selected databases. Annotation
#' entries contain `TERM2GENE`, `TERM2NAME` and `version`.
#' @seealso [PrepareDB], [ListDB]
#' @md
#' @export
#' @examples
#' \dontrun{
#' databases <- PrepareMSigDB(species = "Homo_sapiens", db = "MSigDB")
#' names(databases[["Homo_sapiens"]])
#' }
PrepareMSigDB <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = "MSigDB",
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest", db_update = FALSE, data_dir = NULL,
  convert_species = TRUE, Ensembl_version = NULL, mirror = NULL,
  biomart = NULL, max_tries = 5, verbose = TRUE, ...
) {
  if (!is.character(db) || !length(db) || anyNA(db) || !all(grepl("^MSigDB($|_)", db))) {
    log_message("Unsupported {.arg db} selector for PrepareMSigDB", message_type = "error")
  }
  preparedb_annotation(
    species = species,
    db = db,
    db_IDtypes = db_IDtypes,
    db_version = db_version,
    db_update = db_update,
    data_dir = data_dir,
    convert_species = convert_species,
    Ensembl_version = Ensembl_version,
    mirror = mirror,
    biomart = biomart,
    max_tries = max_tries,
    verbose = verbose, ...,
    prepare = function(db_list, db_species, default_id_types, db, sps, org_sp, org_key, orgdb, biomart) {
      msigdb_requested <- db[grepl("^MSigDB($|_)", db)]
      if (
        length(msigdb_requested) > 0L &&
          any(!msigdb_requested %in% names(db_list[[sps]]))
      ) {
        if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the human annotation to create the MSigDB database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
            db_species["MSigDB"] <- "Homo_sapiens"
          } else {
            log_message(
              "{.pkg MSigDB} database only support {.val {c('Homo_sapiens', 'Mus_musculus')}}. Consider setting {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        log_message("Preparing {.pkg MSigDB} database", verbose = verbose)

        msigdb_release_species <- switch(db_species[["MSigDB"]],
          "Homo_sapiens" = "Hs",
          "Mus_musculus" = "Mm"
        )
        if (is.null(msigdb_release_species)) {
          log_message(
            "{.pkg MSigDB} database only support {.val {c('Homo_sapiens', 'Mus_musculus')}}. Consider setting {.arg convert_species=TRUE}",
            message_type = "error"
          )
        }

        msigdb_pattern <- if (identical(db_version, "latest")) {
          paste0("^msigdb\\.v.+\\.", msigdb_release_species, "\\.json$")
        } else {
          utils::glob2rx(paste0("msigdb.v", db_version, ".json"))
        }
        msigdb_local_file <- preparedb_local_source_file(
          data_dir = data_dir,
          db = "MSigDB",
          pattern = msigdb_pattern,
          verbose = verbose
        )

        if (!is.null(msigdb_local_file)) {
          version <- sub("^msigdb\\.v(.+)\\.json$", "\\1", basename(msigdb_local_file))
          log_message(
            "Using local {.pkg MSigDB} JSON file: {.path {msigdb_local_file}}",
            verbose = verbose
          )
          lines <- paste0(readLines(msigdb_local_file, warn = FALSE), collapse = "")
        } else {
          temp <- tempfile()
          download(
            quiet = TRUE, url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",
            destfile = temp
          )
          version <- readLines(temp)
          version <- version[grep("alt=\"\\[DIR\\]\"", version)]
          version <- version[grep(
            msigdb_release_species,
            version
          )]
          version <- version[length(version)]
          version <- regmatches(
            version,
            m = regexpr("(?<=href\\=\")\\S+(?=/\"\\>)", version, perl = TRUE)
          )

          url <- paste0(
            "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",
            version,
            "/msigdb.v",
            version,
            ".json"
          )
          download(quiet = TRUE, url = url, destfile = temp)
          lines <- paste0(readLines(temp, warn = FALSE), collapse = "")
        }
        lines <- gsub("\"", "", lines)
        lines <- gsub("^\\{|\\}$", "", lines)
        terms <- strsplit(lines, "\\},")[[1]]
        term_list <- list()
        for (idx in seq_along(terms)) {
          term_content <- terms[[idx]]
          term_name <- trimws(gsub(":|_|\\{.*", " ", term_content))
          term_id <- trimws(regmatches(
            term_content,
            m = regexpr(
              "(?<=systematicName:)\\S+?(?=,)",
              term_content,
              perl = TRUE
            )
          ))
          term_gene <- trimws(regmatches(
            term_content,
            m = regexpr(
              pattern = "(?<=geneSymbols:\\[)\\S+?(?=\\],)",
              term_content,
              perl = TRUE
            )
          ))
          term_collection <- trimws(regmatches(
            term_content,
            m = regexpr(
              pattern = "(?<=collection:)\\S+?(?=\\,)",
              term_content,
              perl = TRUE
            )
          ))
          term_list[[idx]] <- c(
            id = term_id,
            "name" = term_name,
            "gene" = term_gene,
            "collection" = term_collection
          )
        }
        df <- as.data.frame(do.call(rbind, term_list))
        df$gene <- strsplit(df$gene, split = ",")
        df <- unnest_fun(df, cols = "gene")

        TERM2NAME <- df[, c(1, 2, 4)]
        TERM2GENE <- df[, c(1, 3)]
        colnames(TERM2NAME) <- c("Term", "Name", "Collection")
        colnames(TERM2GENE) <- c("Term", default_id_types[["MSigDB"]])
        TERM2NAME <- stats::na.omit(unique(TERM2NAME))
        TERM2GENE <- stats::na.omit(unique(TERM2GENE))

        db_list[[db_species["MSigDB"]]][["MSigDB"]][[
          "TERM2GENE"
        ]] <- TERM2GENE
        db_list[[db_species["MSigDB"]]][["MSigDB"]][[
          "TERM2NAME"
        ]] <- TERM2NAME
        db_list[[db_species["MSigDB"]]][["MSigDB"]][["version"]] <- version
        if (sps == db_species["MSigDB"]) {
          preparedb_cache_annotation(
            db_list[[db_species["MSigDB"]]][["MSigDB"]],
            species = as.character(db_species["MSigDB"]), db = "MSigDB"
          )
        }

        msigdb_subsets <- preparedb_msigdb_collection_subsets(
          TERM2GENE = TERM2GENE,
          TERM2NAME = TERM2NAME
        )
        for (collection_db in names(msigdb_subsets)) {
          db_species[collection_db] <- db_species["MSigDB"]
          default_id_types[[collection_db]] <- default_id_types[["MSigDB"]]
          TERM2NAME_sub <- msigdb_subsets[[collection_db]][["TERM2NAME"]]
          TERM2GENE_sub <- msigdb_subsets[[collection_db]][["TERM2GENE"]]
          db_list[[db_species["MSigDB"]]][[collection_db]][[
            "TERM2GENE"
          ]] <- TERM2GENE_sub
          db_list[[db_species["MSigDB"]]][[collection_db]][[
            "TERM2NAME"
          ]] <- TERM2NAME_sub
          db_list[[db_species["MSigDB"]]][[collection_db]][[
            "version"
          ]] <- version
          if (sps == db_species["MSigDB"]) {
            preparedb_cache_annotation(
              db_list[[db_species["MSigDB"]]][[collection_db]],
              species = as.character(db_species["MSigDB"]), db = collection_db
            )
          }
        }
      }
      list(db_list = db_list, db_species = db_species, default_id_types = default_id_types)
    }
  )
}
