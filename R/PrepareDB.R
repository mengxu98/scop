#' @title Prepare databases and reference resources
#'
#' @description
#' Prepare species-specific databases and reference resources from annotation
#' packages, cached downloads, or supplied files.
#'
#' @md
#' @inheritParams GeneConvert
#' @inheritParams thisutils::log_message
#' @param species `"Homo_sapiens"` or `"Mus_musculus"`.
#' @param db Character vector of database or resource selectors: `"GO"`,
#' `"GO_BP"`, `"GO_CC"`, `"GO_MF"`, `"KEGG"`, `"WikiPathway"`, `"Reactome"`,
#' `"CORUM"`, `"MP"`, `"DO"`, `"HPO"`, `"PFAM"`, `"Chromosome"`, `"GeneType"`,
#' `"Enzyme"`, `"TF"`, `"CSPA"`, `"Surfaceome"`, `"SPRomeDB"`, `"VerSeDa"`,
#' `"TFLink"`, `"hTFtarget"`, `"TRRUST"`, `"JASPAR"`, `"ENCODE"`, `"MSigDB"`,
#' `"CellTalk"`, `"CellChat"`, `"CytoTRACE2"` or
#' `"MSigDB_<collection>"`. A vector may combine sources. See the corresponding
#' preparation function for source-specific settings.
#' @param db_IDtypes Gene ID types to include.
#' @param db_version Database version to retrieve.
#' @param db_update Force a refresh. `FALSE` loads the cache when available.
#' @param data_dir Directory or named list of local source files. Searches
#' `data_dir/<db>/` then `data_dir`. Named lists override a path, e.g.
#' `list(MSigDB = "~/db/msigdb")`.
#' @param convert_species Use a species-converted database when the annotation is
#' missing for `species`.
#' @param Ensembl_version Ensembl version. `NULL` uses the latest.
#' @param custom_TERM2GENE,custom_TERM2NAME Custom mappings for `custom_species`.
#' @param custom_species,custom_IDtype,custom_version Metadata for a custom database.
#' @param ... Passed to helper functions.
#'
#' @return A named list of prepared resources. Gene annotation databases are
#' nested by species and database, with `TERM2GENE` (gene-to-term mappings),
#' `TERM2NAME` (term names) and `version`. Additional resource fields follow the
#' return contract of the corresponding preparation function.
#'
#' @seealso [ListDB], [PrepareGO], [PrepareKEGG], [PrepareMSigDB]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' databases <- PrepareDB(species = "Homo_sapiens", db = c("GO_BP", "KEGG"))
#' names(databases[["Homo_sapiens"]])
#' ListDB(species = "Homo_sapiens", db = c("GO_BP", "KEGG"))
#' }
PrepareDB <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db = c(
    "GO",
    "GO_BP",
    "GO_CC",
    "GO_MF",
    "KEGG",
    "WikiPathway",
    "Reactome",
    "CORUM",
    "MP",
    "DO",
    "HPO",
    "PFAM",
    "CSPA",
    "Surfaceome",
    "SPRomeDB",
    "VerSeDa",
    "TFLink",
    "hTFtarget",
    "TRRUST",
    "JASPAR",
    "ENCODE",
    "MSigDB",
    "CellTalk",
    "CellChat",
    "Chromosome",
    "GeneType",
    "Enzyme",
    "TF",
    "CytoTRACE2"
  ),
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest",
  db_update = FALSE,
  data_dir = NULL,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  biomart = NULL,
  max_tries = 5,
  custom_TERM2GENE = NULL,
  custom_TERM2NAME = NULL,
  custom_species = NULL,
  custom_IDtype = NULL,
  custom_version = NULL,
  verbose = TRUE,
  ...
) {
  species <- normalize_species_name(species)
  if (!is.null(db)) {
    db <- as.character(db)
    if (anyNA(db)) log_message("{.arg db} cannot contain missing values", message_type = "error")
  }
  if (!is.null(custom_TERM2GENE)) {
    return(PrepareCustomDB(
      species, db, db_IDtypes, db_version, db_update, data_dir,
      convert_species, Ensembl_version, mirror, biomart, max_tries, custom_TERM2GENE,
      custom_TERM2NAME, custom_species, custom_IDtype, custom_version, verbose, ...
    ))
  }
  db_requested <- db
  db <- normalize_msigdb_db_names(db)
  db_list <- list()
  if ("CytoTRACE2" %in% db) {
    db_list <- PrepareCytoTRACE2(db_update = db_update, verbose = verbose)
    db <- setdiff(db, "CytoTRACE2")
  }
  annotation_db <- db
  cached_order <- list()
  if (length(annotation_db) && isFALSE(db_update)) {
    check_r("R.cache", verbose = FALSE)
    info <- list_db_cache_entries(species = species, db = annotation_db)
    if (!is.null(info) && nrow(info)) {
      for (sps in species) {
        cached <- info$DB[info$Species == sps]
        cached_order[[sps]] <- annotation_db[vapply(
          annotation_db,
          function(term) any(grepl(paste0("^", term, "$"), cached)), logical(1)
        )]
      }
    }
  }
  sources <- list(
    GO = PrepareGO,
    KEGG = PrepareKEGG,
    WikiPathway = PrepareWikiPathway,
    Reactome = PrepareReactome,
    CORUM = PrepareCORUM,
    MP = PrepareMP,
    DO = PrepareDO,
    HPO = PrepareHPO,
    PFAM = PreparePFAM,
    Chromosome = PrepareChromosome,
    GeneType = PrepareGeneType,
    Enzyme = PrepareEnzyme,
    TF = PrepareTF,
    CSPA = PrepareCSPA,
    Surfaceome = PrepareSurfaceome,
    SPRomeDB = PrepareSPRomeDB,
    VerSeDa = PrepareVerSeDa,
    TFLink = PrepareTFLink,
    hTFtarget = PrepareHTFtarget,
    TRRUST = PrepareTRRUST,
    JASPAR = PrepareJASPAR,
    ENCODE = PrepareENCODE,
    MSigDB = PrepareMSigDB,
    CellTalk = PrepareCellTalk,
    CellChat = PrepareCellChat
  )
  source_names <- ifelse(grepl("^MSigDB_", annotation_db), "MSigDB", annotation_db)
  source_names[source_names %in% c("GO", "GO_BP", "GO_CC", "GO_MF")] <- "GO"
  prepared_sources <- unique(source_names)
  prepared_sources <- c(intersect(names(sources), prepared_sources), setdiff(prepared_sources, names(sources)))
  for (source in prepared_sources) {
    selected <- annotation_db[source_names == source]
    prepare <- sources[[source]]
    if (is.null(prepare)) prepare <- preparedb_annotation
    prepared <- prepare(
      species = species, db = selected, db_IDtypes = db_IDtypes,
      db_version = db_version, db_update = db_update, data_dir = data_dir,
      convert_species = convert_species, Ensembl_version = Ensembl_version,
      mirror = mirror, biomart = biomart, max_tries = max_tries, verbose = verbose, ...
    )
    for (sps in names(prepared)) {
      db_list[[sps]] <- c(db_list[[sps]], prepared[[sps]])
    }
  }
  db_requested <- setdiff(db_requested, "CytoTRACE2")
  for (sps in intersect(species, names(db_list))) {
    ordered <- unique(c(cached_order[[sps]], names(db_list[[sps]])))
    db_list[[sps]] <- db_list[[sps]][intersect(ordered, names(db_list[[sps]]))]
    db_list <- preparedb_alias_requested_msigdb(db_list, sps, db_requested, db)
  }
  db_list
}

#' @title Prepare custom gene annotation databases
#' @description Build custom mappings using the common cache and species/identifier conversion pipeline.
#' @inheritParams PrepareDB
#' @return A named species list containing the selected annotation database.
#' @seealso [PrepareDB], [ListDB]
#' @export
#' @examples
#' mappings <- data.frame(Term = c("Response", "Response"), symbol = c("Isg15", "Ifit3"))
#' databases <- PrepareCustomDB(
#'   species = "Mus_musculus", db = "Response", db_IDtypes = "symbol",
#'   custom_TERM2GENE = mappings, custom_species = "Mus_musculus",
#'   custom_IDtype = "symbol", custom_version = "v1", verbose = FALSE
#' )
#' databases[["Mus_musculus"]][["Response"]]$TERM2GENE
PrepareCustomDB <- function(
  species = c("Homo_sapiens", "Mus_musculus"), db,
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"), db_version = "latest",
  db_update = FALSE, data_dir = NULL, convert_species = TRUE, Ensembl_version = NULL,
  mirror = NULL, biomart = NULL, max_tries = 5, custom_TERM2GENE = NULL,
  custom_TERM2NAME = NULL, custom_species = NULL, custom_IDtype = NULL,
  custom_version = NULL, verbose = TRUE, ...
) {
  preparedb_annotation(
    species = species, db = db, db_IDtypes = db_IDtypes,
    db_version = db_version, db_update = db_update, data_dir = data_dir,
    convert_species = convert_species, Ensembl_version = Ensembl_version,
    mirror = mirror, biomart = biomart, max_tries = max_tries,
    custom_TERM2GENE = custom_TERM2GENE, custom_TERM2NAME = custom_TERM2NAME,
    custom_species = custom_species, custom_IDtype = custom_IDtype,
    custom_version = custom_version, verbose = verbose, ...
  )
}

preparedb_annotation <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  db,
  db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
  db_version = "latest",
  db_update = FALSE,
  data_dir = NULL,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  biomart = NULL,
  max_tries = 5,
  custom_TERM2GENE = NULL,
  custom_TERM2NAME = NULL,
  custom_species = NULL,
  custom_IDtype = NULL,
  custom_version = NULL,
  verbose = TRUE,
  prepare = NULL,
  ...
) {
  check_r("R.cache", verbose = FALSE)
  species <- normalize_species_name(species)
  db_list <- list()

  if (!is.null(db)) {
    db <- as.character(db)
    if (anyNA(db)) log_message("{.arg db} cannot contain missing values", message_type = "error")
  }
  db_requested <- db
  db <- normalize_msigdb_db_names(db)

  for (sps in species) {
    log_message(
      "Species: {.val {sps}}",
      verbose = verbose
    )
    default_id_types <- list(
      "GO" = "entrez_id",
      "GO_BP" = "entrez_id",
      "GO_CC" = "entrez_id",
      "GO_MF" = "entrez_id",
      "KEGG" = "entrez_id",
      "WikiPathway" = "entrez_id",
      "Reactome" = "entrez_id",
      "CORUM" = "symbol",
      "MP" = "symbol",
      "DO" = "symbol",
      "HPO" = "symbol",
      "PFAM" = "entrez_id",
      "Chromosome" = "entrez_id",
      "GeneType" = "entrez_id",
      "Enzyme" = "entrez_id",
      "TF" = "symbol",
      "CSPA" = "symbol",
      "Surfaceome" = "symbol",
      "SPRomeDB" = "entrez_id",
      "VerSeDa" = "symbol",
      "TFLink" = "symbol",
      "hTFtarget" = "symbol",
      "TRRUST" = "symbol",
      "JASPAR" = "symbol",
      "ENCODE" = "symbol",
      "MSigDB" = "symbol",
      "CellTalk" = "symbol",
      "CellChat" = "symbol"
    )
    if (!is.null(custom_TERM2GENE)) {
      if (length(db) > 1) {
        log_message(
          "When building a custom database, the length of {.arg db} must be 1",
          message_type = "error"
        )
      }
      if (is.null(custom_IDtype) || is.null(custom_species) || is.null(custom_version)) {
        log_message(
          "When building a custom database, {.arg custom_IDtype}, {.arg custom_species} and {.arg custom_version} must be provided",
          message_type = "error"
        )
      }
      custom_IDtype <- match.arg(
        custom_IDtype,
        choices = c("symbol", "entrez_id", "ensembl_id")
      )
      default_id_types[[db]] <- custom_IDtype
    }

    if (isFALSE(db_update) && is.null(custom_TERM2GENE)) {
      for (term in db) {
        cached <- preparedb_load_cache_entry(sps, term, db_version, verbose)
        if (!is.null(cached)) db_list[[sps]][[term]] <- cached
      }
    }

    db_species <- stats::setNames(object = rep(sps, length(db)), nm = db)
    if (any(grepl("^MSigDB($|_)", db)) && !"MSigDB" %in% names(db_species)) {
      db_species["MSigDB"] <- sps
    }

    sp <- unlist(strsplit(sps, split = "_"))
    org_sp <- paste0(
      "org.",
      paste0(substring(sp, 1, 1), collapse = ""),
      ".eg.db"
    )
    org_key <- "ENTREZID"
    if (sps == "Arabidopsis_thaliana") {
      biomart <- "plants_mart"
      org_sp <- "org.At.tair.db"
      org_key <- "TAIR"
      default_id_types[c(
        "GO",
        "GO_BP",
        "GO_CC",
        "GO_MF",
        "PFAM",
        "Chromosome",
        "GeneType",
        "Enzyme"
      )] <- "tair_locus"
    }
    if (sps == "Saccharomyces_cerevisiae") {
      org_sp <- "org.Sc.sgd.db"
    }

    if (any(!sps %in% names(db_list)) || any(!db %in% names(db_list[[sps]]))) {
      orgdb_dependent <- c(
        "GO",
        "GO_BP",
        "GO_CC",
        "GO_MF",
        "PFAM",
        "Chromosome",
        "GeneType",
        "Enzyme"
      )
      if (any(orgdb_dependent %in% db)) {
        check_r(c(org_sp, "GO.db", "GOSemSim"), verbose = FALSE)
        if (!isTRUE(all(unlist(check_r(org_sp, install = FALSE, verbose = FALSE), use.names = FALSE)))) {
          log_message(
            "Annotation package {.pkg {org_sp}} is not installed. Install it with {.code BiocManager::install('{org_sp}')}.",
            message_type = "warning"
          )
          if (isTRUE(convert_species)) {
            db_to_convert <- intersect(db, orgdb_dependent)
            log_message(
              "Use the human annotation to create the {.pkg {db_to_convert}} database for {.val {sps}}",
              message_type = "warning"
            )
            org_sp <- "org.Hs.eg.db"
            db_species[db_to_convert] <- "Homo_sapiens"
          } else {
            log_message(
              "Stop the preparation",
              message_type = "error"
            )
            log_message(
              "Required annotation package is not available: ",
              org_sp,
              message_type = "error"
            )
          }
        }
        check_r(org_sp, verbose = FALSE)
        orgdb <- get_namespace_fun(org_sp, org_sp)
      }
      if ("PFAM" %in% db) {
        check_r("PFAM.db", verbose = FALSE)
      }
      if ("Reactome" %in% db) {
        check_r("reactome.db", verbose = FALSE)
      }

      if (is.null(custom_TERM2GENE)) {
        if (!is.null(prepare)) {
          prepared <- prepare(
            db_list, db_species, default_id_types, db, sps,
            org_sp, org_key, if (exists("orgdb", inherits = FALSE)) orgdb else NULL, biomart
          )
          if (is.null(prepared)) next
          db_list <- prepared$db_list
          db_species <- prepared$db_species
          default_id_types <- prepared$default_id_types
        }
      } else {
        db_species[db] <- custom_species
        if (sps != custom_species) {
          if (isTRUE(convert_species)) {
            log_message(
              "Use the {.val {custom_species}} annotation to create the {.val {db}} database for {.val {sps}}",
              message_type = "warning",
              verbose = verbose
            )
          } else {
            log_message(
              "{.pkg {db}} database only support {.val {custom_species}}. Consider using {.arg convert_species=TRUE}",
              message_type = "error"
            )
          }
        }
        custom_db <- normalize_custom_db_input(
          TERM2GENE = custom_TERM2GENE,
          TERM2NAME = custom_TERM2NAME,
          IDtype = custom_IDtype,
          remove_na = TRUE
        )
        TERM2GENE <- custom_db[["TERM2GENE"]]
        TERM2NAME <- custom_db[["TERM2NAME"]]
        db_list[[db_species[db]]][[db]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species[db]]][[db]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species[db]]][[db]][["version"]] <- custom_version
        if (sps == db_species[db]) {
          preparedb_cache_annotation(
            db_list[[db_species[db]]][[db]],
            species = as.character(db_species[db]), db = db
          )
        }
      }
    }

    db_list <- preparedb_ensure_msigdb_parent(
      db_list = db_list,
      species = sps,
      db_version = db_version,
      verbose = verbose
    )

    if (!all(db_species == sps)) {
      for (term in names(db_species[db_species != sps])) {
        log_message(
          "Convert species for the {.pkg {term}} database",
          verbose = verbose
        )
        sp_from <- db_species[term]
        db_info <- db_list[[sp_from]][[names(sp_from)]]
        TERM2GENE <- preparedb_normalize_term2gene_id_columns(
          db_info[["TERM2GENE"]]
        )
        db_info[["TERM2GENE"]] <- TERM2GENE
        TERM2NAME <- db_info[["TERM2NAME"]]
        IDtype <- preparedb_source_idtype(
          term = term,
          TERM2GENE = TERM2GENE,
          default_id_types = default_id_types
        )
        TERM2GENE_map <- NULL
        if (grepl("^MSigDB_", term)) {
          TERM2GENE_map <- db_list[[sps]][["MSigDB"]][["TERM2GENE"]]
        }
        if (!is.null(TERM2GENE_map)) {
          TERM2GENE <- TERM2GENE_map[
            TERM2GENE_map[["Term"]] %in% TERM2GENE[["Term"]], ,
            drop = FALSE
          ]
          TERM2NAME <- TERM2NAME[
            TERM2NAME[["Term"]] %in% TERM2GENE[["Term"]], ,
            drop = FALSE
          ]
        } else {
          res <- GeneConvert(
            geneID = as.character(unique(TERM2GENE[, 2])),
            geneID_from_IDtype = IDtype,
            geneID_to_IDtype = "ensembl_id",
            species_from = sp_from,
            species_to = sps,
            Ensembl_version = Ensembl_version,
            mirror = mirror,
            biomart = biomart,
            max_tries = max_tries
          )
          if (is.null(res$geneID_res)) {
            log_message(
              "Failed to convert species for the database: {.val {term}}",
              message_type = "warning",
              verbose = verbose
            )
            next
          }
          map <- res$geneID_collapse
          TERM2GENE[["ensembl_id-converted"]] <- map[
            as.character(TERM2GENE[, 2]),
            "ensembl_id"
          ]
          TERM2GENE <- unnest_fun(
            TERM2GENE,
            cols = "ensembl_id-converted",
            keep_empty = FALSE
          )
          TERM2GENE <- TERM2GENE[, c("Term", "ensembl_id-converted")]
          colnames(TERM2GENE) <- c("Term", "ensembl_id")
          TERM2NAME <- TERM2NAME[
            TERM2NAME[["Term"]] %in% TERM2GENE[["Term"]], ,
            drop = FALSE
          ]
        }

        db_info[["TERM2GENE"]] <- unique(TERM2GENE)
        db_info[["TERM2NAME"]] <- unique(TERM2NAME)
        version <- paste0(
          db_info[["version"]],
          "(converted from ",
          sp_from,
          ")"
        )
        db_info[["version"]] <- version
        db_list[[sps]][[term]] <- db_info
        default_id_types[[term]] <- "ensembl_id"
        preparedb_cache_annotation(
          db_list[[sps]][[term]],
          species = sps, db = term
        )
      }
    }

    for (term in names(db_list[[sps]])) {
      if (is.null(db_list[[sps]][[term]][["TERM2GENE"]])) next
      db_list[[sps]][[term]][["TERM2GENE"]] <-
        preparedb_normalize_term2gene_id_columns(
          db_list[[sps]][[term]][["TERM2GENE"]]
        )
      IDtypes <- db_IDtypes[
        !db_IDtypes %in% colnames(db_list[[sps]][[term]][["TERM2GENE"]])
      ]
      if (length(IDtypes) > 0) {
        log_message(
          "Convert ID types for the {.pkg {term}} database",
          verbose = verbose
        )
        TERM2GENE <- db_list[[sps]][[term]][["TERM2GENE"]]
        TERM2NAME <- db_list[[sps]][[term]][["TERM2NAME"]]
        IDtype <- preparedb_source_idtype(
          term = term,
          TERM2GENE = TERM2GENE,
          default_id_types = default_id_types
        )
        parent_term2gene <- preparedb_msigdb_parent_term2gene(
          parent = db_list[[sps]][["MSigDB"]][["TERM2GENE"]],
          id_types = IDtypes
        )
        if (grepl("^MSigDB_", term) && !is.null(parent_term2gene)) {
          map <- parent_term2gene[, -1, drop = FALSE]
          map <- stats::aggregate(
            map,
            by = list(map[[1]]),
            FUN = function(x) list(unique(x))
          )
          rownames(map) <- map[, 1]
        } else {
          map <- preparedb_local_orgdb_id_map(
            geneID = as.character(unique(TERM2GENE[, 2])),
            geneID_from_IDtype = IDtype,
            geneID_to_IDtype = IDtypes,
            org_sp = org_sp,
            org_key = org_key,
            verbose = verbose
          )
          if (is.null(map)) {
            res <- GeneConvert(
              geneID = as.character(unique(TERM2GENE[, 2])),
              geneID_from_IDtype = IDtype,
              geneID_to_IDtype = IDtypes,
              species_from = sps,
              species_to = sps,
              Ensembl_version = Ensembl_version,
              mirror = mirror,
              biomart = biomart,
              max_tries = max_tries
            )
            if (is.null(res$geneID_res)) {
              log_message(
                "Failed to convert ID types for the database: {.val {term}}",
                message_type = "warning",
                verbose = verbose
              )
              next
            }
            map <- res$geneID_collapse
          }
          if (is.null(map) || nrow(map) == 0) {
            log_message(
              "Failed to convert ID types for the database: {.val {term}}",
              message_type = "warning",
              verbose = verbose
            )
            next
          }
        }
        for (type in IDtypes) {
          TERM2GENE[[type]] <- map[as.character(TERM2GENE[, 2]), type]
          TERM2GENE <- unnest_fun(TERM2GENE, cols = type, keep_empty = TRUE)
        }
        db_list[[sps]][[term]][["TERM2GENE"]] <- TERM2GENE
        version <- db_list[[sps]][[term]][["version"]]
        preparedb_cache_annotation(
          db_list[[sps]][[term]],
          species = sps, db = term
        )
      }
    }
    preparedb_require_msigdb_names(
      db_list = db_list,
      species = sps,
      db = db
    )
    db_list <- preparedb_alias_requested_msigdb(
      db_list = db_list,
      species = sps,
      db_requested = db_requested,
      db_normalized = db
    )
  }
  db_list <- preparedb_keep_requested_dbs(
    db_list = db_list,
    species = species,
    db_names = unique(c(db, db_requested))
  )
  return(db_list)
}

normalize_species_name <- function(species) {
  if (
    !is.character(species) ||
      length(species) == 0L ||
      any(is.na(species)) ||
      any(!nzchar(trimws(species)))
  ) {
    log_message(
      "{.arg species} must be a non-empty character vector",
      message_type = "error"
    )
  }

  vapply(species, function(x) {
    x <- trimws(x)
    x <- gsub("[[:space:].-]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0L) {
      return(x)
    }
    parts <- tolower(parts)
    parts[1] <- paste0(
      toupper(substr(parts[1], 1, 1)),
      substr(parts[1], 2, nchar(parts[1]))
    )
    paste(parts, collapse = "_")
  }, character(1), USE.NAMES = FALSE)
}

preparedb_normalize_term2gene_id_columns <- function(TERM2GENE) {
  legacy_msigdb_col <- "symbol.ensembl_id"
  if (
    legacy_msigdb_col %in% colnames(TERM2GENE) &&
      !"symbol" %in% colnames(TERM2GENE)
  ) {
    colnames(TERM2GENE)[colnames(TERM2GENE) == legacy_msigdb_col] <- "symbol"
  }
  TERM2GENE
}

normalize_msigdb_db_names <- function(db) {
  if (is.null(db) || length(db) == 0L) {
    return(db)
  }
  db <- as.character(db)
  nested <- grepl(
    "^MSigDB_[A-Za-z0-9_]+(?::[A-Za-z0-9_]+)+$",
    db,
    perl = TRUE
  )
  db[nested] <- gsub(":", "_", db[nested], fixed = TRUE)
  db
}

msigdb_collection_db_names <- function(collection) {
  collections <- as.character(collection)
  unlist(lapply(collections, function(x) {
    if (is.na(x) || !nzchar(x)) {
      return(character(0))
    }
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0L) {
      return(character(0))
    }
    vapply(seq_along(parts), function(i) {
      paste0("MSigDB_", paste(parts[seq_len(i)], collapse = "_"))
    }, character(1))
  }), use.names = FALSE)
}

preparedb_msigdb_collection_subsets <- function(TERM2GENE, TERM2NAME) {
  collections <- unique(as.character(TERM2NAME[["Collection"]]))
  collections <- collections[!is.na(collections) & nzchar(collections)]
  db_names <- unique(msigdb_collection_db_names(collections))
  normalized <- gsub(":", "_", as.character(TERM2NAME[["Collection"]]), fixed = TRUE)
  stats::setNames(lapply(db_names, function(db_name) {
    suffix <- sub("^MSigDB_", "", db_name)
    keep <- !is.na(normalized) & (
      normalized == suffix | startsWith(normalized, paste0(suffix, "_"))
    )
    term2name <- TERM2NAME[keep, , drop = FALSE]
    term2gene <- TERM2GENE[
      TERM2GENE[["Term"]] %in% term2name[["Term"]], ,
      drop = FALSE
    ]
    list(TERM2GENE = term2gene, TERM2NAME = term2name)
  }), db_names)
}

preparedb_require_msigdb_names <- function(db_list, species, db) {
  requested <- as.character(db)
  requested <- requested[grepl("^MSigDB($|_)", requested)]
  if (length(requested) == 0L) {
    return(invisible(NULL))
  }
  available_names <- names(db_list[[species]])
  missing <- requested[!requested %in% available_names]
  if (length(missing) == 0L) {
    return(invisible(NULL))
  }
  available <- grep("^MSigDB($|_)", available_names, value = TRUE)
  if (length(available) == 0L) {
    available <- "(none)"
  }
  log_message(
    paste0(
      "MSigDB database {.val {missing}} is not available. ",
      "Nested collections use underscores, for example {.val MSigDB_M2_CGP}, ",
      "{.val MSigDB_M2_CP} and {.val MSigDB_M2_CP_BIOCARTA}. ",
      "Available MSigDB databases: {.val {available}}."
    ),
    message_type = "error"
  )
}

preparedb_cache_annotation <- function(object, species, db) {
  R.cache::saveCache(object,
    key = list(object$version, as.character(species), db),
    comment = paste0(
      object$version, " nterm:", length(object$TERM2NAME[[1]]),
      "|", species, "-", db
    )
  )
}

preparedb_load_cache_entry <- function(
  species,
  db,
  db_version = "latest",
  verbose = TRUE
) {
  dbinfo <- list_db_cache_entries(species = species, db = db)
  if (is.null(dbinfo) || nrow(dbinfo) == 0L) {
    return(NULL)
  }
  if (identical(as.character(db_version), "latest")) {
    pathname <- dbinfo[
      order(dbinfo[["timestamp"]], decreasing = TRUE)[1],
      "file"
    ]
  } else {
    pathname <- dbinfo[
      grep(db_version, dbinfo[["db_version"]], fixed = TRUE)[1],
      "file"
    ]
    if (is.na(pathname)) {
      log_message(
        "There is no {.val {db_version}} version of the database. Use the latest version",
        message_type = "warning",
        verbose = verbose
      )
      pathname <- dbinfo[
        order(dbinfo[["timestamp"]], decreasing = TRUE)[1],
        "file"
      ]
    }
  }
  if (length(pathname) == 0L || is.na(pathname)) {
    return(NULL)
  }
  header <- R.cache::readCacheHeader(pathname)
  cached_version <- strsplit(header[["comment"]], "\\|")[[1]][1]
  timestamp <- format(header[["timestamp"]], "%Y-%m-%d %H:%M:%S")
  log_message(
    "Loading cached: {.pkg {db}} version: {.pkg {cached_version}} created: {.pkg {timestamp}}",
    verbose = verbose
  )
  R.cache::loadCache(pathname = pathname)
}

preparedb_ensure_msigdb_parent <- function(
  db_list,
  species,
  db_version = "latest",
  verbose = TRUE
) {
  present <- names(db_list[[species]])
  needs_parent <- any(grepl("^MSigDB_", present)) && !"MSigDB" %in% present
  if (!isTRUE(needs_parent)) {
    return(db_list)
  }
  parent <- preparedb_load_cache_entry(
    species = species,
    db = "MSigDB",
    db_version = db_version,
    verbose = verbose
  )
  if (!is.null(parent)) {
    db_list[[species]][["MSigDB"]] <- parent
  }
  db_list
}

preparedb_msigdb_parent_term2gene <- function(parent, id_types) {
  if (
    is.null(parent) ||
      !is.data.frame(parent) ||
      length(id_types) == 0L ||
      !all(id_types %in% colnames(parent))
  ) {
    return(NULL)
  }
  parent
}

preparedb_keep_requested_dbs <- function(db_list, species, db_names) {
  db_names <- unique(as.character(db_names))
  db_names <- db_names[!is.na(db_names) & nzchar(db_names)]
  if (length(db_names) == 0L) {
    return(db_list)
  }
  for (sps in intersect(as.character(species), names(db_list))) {
    keep <- intersect(names(db_list[[sps]]), db_names)
    db_list[[sps]] <- db_list[[sps]][keep]
  }
  drop <- setdiff(names(db_list), c(as.character(species), "CytoTRACE2"))
  if (length(drop) > 0L) {
    db_list[drop] <- NULL
  }
  db_list
}

preparedb_alias_requested_msigdb <- function(
  db_list,
  species,
  db_requested,
  db_normalized
) {
  if (length(db_requested) == 0L || is.null(db_list[[species]])) {
    return(db_list)
  }
  for (i in seq_along(db_requested)) {
    requested <- db_requested[[i]]
    normalized <- db_normalized[[i]]
    if (identical(requested, normalized)) {
      next
    }
    entry <- db_list[[species]][[normalized]]
    if (!is.null(entry)) {
      db_list[[species]][[requested]] <- entry
    }
  }
  db_list
}

preparedb_local_source_file <- function(
  data_dir,
  db,
  pattern,
  verbose = TRUE
) {
  if (is.null(data_dir)) {
    return(NULL)
  }

  data_source <- data_dir
  if (is.list(data_dir)) {
    if (is.null(names(data_dir)) || !db %in% names(data_dir)) {
      return(NULL)
    }
    data_source <- data_dir[[db]]
    if (is.list(data_source)) {
      if (!is.null(data_source[["file"]])) {
        data_source <- data_source[["file"]]
      } else if (!is.null(data_source[["path"]])) {
        data_source <- data_source[["path"]]
      } else {
        return(NULL)
      }
    }
  }

  if (!is.character(data_source) || length(data_source) != 1 || is.na(data_source) || !nzchar(data_source)) {
    log_message(
      "{.arg data_dir} must be one directory path, one file path, or a named list of paths",
      message_type = "error"
    )
  }
  data_source <- path.expand(data_source)
  if (file.exists(data_source) && !dir.exists(data_source)) {
    if (any(grepl(pattern, basename(data_source), perl = TRUE))) {
      return(normalizePath(data_source, mustWork = TRUE))
    }
    log_message(
      "Local source file {.path {data_source}} does not match the expected {.pkg {db}} file name; downloading instead",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
  if (!dir.exists(data_source)) {
    log_message(
      "{.arg data_dir} must point to an existing local directory or file",
      message_type = "error"
    )
  }

  candidate_dirs <- unique(c(file.path(data_source, db), data_source))
  candidate_dirs <- candidate_dirs[dir.exists(candidate_dirs)]
  local_files <- unlist(lapply(candidate_dirs, function(dir) {
    list.files(
      dir,
      pattern = pattern,
      full.names = TRUE
    )
  }), use.names = FALSE)
  local_files <- unique(local_files)
  if (length(local_files) == 0) {
    log_message(
      "No local source file for {.pkg {db}} was found in {.arg data_dir}; downloading instead",
      message_type = "warning",
      verbose = verbose
    )
    return(NULL)
  }
  local_files[order(file.info(local_files)[["mtime"]], decreasing = TRUE)][[1]]
}

preparedb_read_gmt_source <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    temp <- tempfile(fileext = ".gz")
    file.copy(path, temp, overwrite = TRUE)
    R.utils::gunzip(temp)
    unzipped <- sub("\\.gz$", "", temp)
    on.exit(unlink(unzipped), add = TRUE)
    return(clusterProfiler::read.gmt(unzipped))
  }
  clusterProfiler::read.gmt(path)
}

preparedb_source_idtype <- function(term, TERM2GENE, default_id_types) {
  default_idtype <- default_id_types[[term]]
  if (length(default_idtype) > 0 && !all(is.na(default_idtype))) {
    default_idtype <- default_idtype[default_idtype %in% colnames(TERM2GENE)]
    if (length(default_idtype) > 0) {
      return(default_idtype[[1]])
    }
  }
  colnames(TERM2GENE)[[2]]
}

preparedb_local_orgdb_id_map <- function(
  geneID,
  geneID_from_IDtype,
  geneID_to_IDtype,
  org_sp,
  org_key,
  verbose = TRUE
) {
  if (is.null(org_sp) || !isTRUE(all(unlist(check_r(org_sp, install = FALSE, verbose = FALSE), use.names = FALSE)))) {
    return(NULL)
  }
  idtype_to_orgdb_column <- function(idtype) {
    if (length(idtype) != 1) {
      return(NA_character_)
    }
    switch(tolower(idtype),
      "symbol" = "SYMBOL",
      "ensembl_id" = "ENSEMBL",
      "entrez_id" = org_key,
      "tair_locus" = "TAIR",
      "sgd_gene" = "SGD",
      NA_character_
    )
  }
  from_column <- idtype_to_orgdb_column(geneID_from_IDtype)
  to_columns <- stats::setNames(
    vapply(geneID_to_IDtype, idtype_to_orgdb_column, character(1)),
    geneID_to_IDtype
  )
  if (is.na(from_column) || any(is.na(to_columns))) {
    return(NULL)
  }
  orgdb <- get_namespace_fun(org_sp, org_sp)
  columns_available <- AnnotationDbi::columns(orgdb)
  columns_needed <- unique(c(from_column, to_columns))
  if (any(!columns_needed %in% columns_available)) {
    return(NULL)
  }
  geneID <- unique(stats::na.omit(as.character(geneID)))
  geneID <- geneID[nzchar(geneID)]
  if (length(geneID) == 0) {
    return(NULL)
  }
  map <- tryCatch(
    {
      suppressMessages(
        AnnotationDbi::select(
          orgdb,
          keys = geneID,
          keytype = from_column,
          columns = unique(to_columns)
        )
      )
    },
    error = function(e) NULL
  )
  if (is.null(map) || nrow(map) == 0) {
    return(NULL)
  }
  map <- map[
    !is.na(map[[from_column]]) & map[[from_column]] %in% geneID, ,
    drop = FALSE
  ]
  if (nrow(map) == 0) {
    return(NULL)
  }
  for (type in names(to_columns)) {
    if (!identical(type, to_columns[[type]])) {
      map[[type]] <- map[[to_columns[[type]]]]
    }
  }
  map <- map[, unique(c(from_column, names(to_columns))), drop = FALSE]
  map <- stats::aggregate(
    map[, names(to_columns), drop = FALSE],
    by = list(map[[from_column]]),
    FUN = function(x) {
      list(unique(x[!is.na(x) & nzchar(as.character(x))]))
    }
  )
  rownames(map) <- map[, 1]
  map <- map[, -1, drop = FALSE]
  log_message(
    "Converted ID types using local annotation package {.pkg {org_sp}}",
    verbose = verbose
  )
  map
}

kegg_get <- function(url) {
  temp <- tempfile()
  on.exit(unlink(temp))
  download(quiet = TRUE, url = url, destfile = temp)
  content <- as.data.frame(
    do.call(
      rbind,
      strsplit(readLines(temp), split = "\t")
    )
  )
  content
}

kegg_release_version <- function(info_lines) {
  kegg_release_from_info(info_lines) %||%
    kegg_release_from_relnote() %||%
    paste0("Retrieved ", Sys.Date())
}

kegg_release_from_info <- function(info_lines) {
  release_line <- info_lines[grepl("Release", x = info_lines)]
  if (length(release_line) == 0) {
    return(NULL)
  }
  gsub(".*(?=Release)", replacement = "", x = release_line, perl = TRUE)
}

kegg_release_from_relnote <- function(url = "https://www.kegg.jp/kegg/docs/relnote.html") {
  page <- tryCatch(
    httr::content(httr::GET(url, httr::timeout(10))),
    error = function(e) NULL
  )
  if (is.null(page)) {
    return(NULL)
  }
  page <- paste0(page, collapse = "")
  release <- regmatches(page, regexpr("Release [0-9]+\\.[0-9]+", page))
  if (length(release) == 0) NULL else release
}

#' @title Prepare CytoTRACE2 model resources
#' @description Download or load the versioned model assets used by [RunCytoTRACE()].
#' @inheritParams PrepareDB
#' @return A list with a `CytoTRACE2` entry containing `data_dir`, `files` and `version`.
#' @details Model assets are species-independent and cached in the user data
#' directory. This entry retains the structure used by [RunCytoTRACE()].
#' @seealso [PrepareDB], [RunCytoTRACE]
#' @export
#' @examples
#' \dontrun{
#' model <- PrepareCytoTRACE2()
#' model[["CytoTRACE2"]]$version
#' }
PrepareCytoTRACE2 <- function(db_update = FALSE, verbose = TRUE) {
  check_r("R.cache", verbose = FALSE)
  db_list <- list()
  cyto_version <- "1.1.0"
  cyto_cache_key <- list(cyto_version, "CytoTRACE2", "CytoTRACE2")
  cyto_data_dir <- file.path(
    tools::R_user_dir("scop", "data"),
    "CytoTRACE2"
  )
  cyto_files <- c(
    "model_parameters.rds",
    "features_model_training_17.csv",
    "mt_dict_human_to_mouse.csv",
    "mt_human_alias.csv",
    "mt_mouse_alias.csv"
  )
  cyto_url <- "https://raw.githubusercontent.com/mengxu98/datasets/main/CytoTRACE2"

  if (isFALSE(db_update)) {
    cyto_cached <- R.cache::loadCache(key = cyto_cache_key)
    cyto_cached_dir <- if (!is.null(cyto_cached$data_dir)) {
      normalizePath(cyto_cached$data_dir, mustWork = FALSE)
    } else {
      NULL
    }
    cyto_data_dir_norm <- normalizePath(cyto_data_dir, mustWork = FALSE)
    if (
      !is.null(cyto_cached) &&
        identical(cyto_cached_dir, cyto_data_dir_norm) &&
        dir.exists(cyto_cached_dir) &&
        all(file.exists(file.path(cyto_cached_dir, cyto_files)))
    ) {
      log_message(
        "Loading cached: {.pkg CytoTRACE2} version: {.pkg {cyto_version}}",
        verbose = verbose
      )
      db_list[["CytoTRACE2"]] <- cyto_cached
    } else if (!is.null(cyto_cached_dir) && dir.exists(cyto_cached_dir)) {
      log_message(
        "Ignoring legacy CytoTRACE2 cache outside the datasets cache: {.path {cyto_cached_dir}}",
        verbose = verbose
      )
    }
  }

  if (is.null(db_list[["CytoTRACE2"]])) {
    log_message(
      "Preparing {.pkg CytoTRACE2} database",
      verbose = verbose
    )

    if (!dir.exists(cyto_data_dir) ||
      !all(file.exists(file.path(cyto_data_dir, cyto_files))) ||
      isTRUE(db_update)) {
      log_message(
        "Downloading CytoTRACE2 model data from datasets GitHub repository...",
        verbose = verbose
      )
      dir.create(cyto_data_dir, showWarnings = FALSE, recursive = TRUE)
      old_timeout <- getOption("timeout")
      options(timeout = max(600, old_timeout))
      on.exit(options(timeout = old_timeout), add = TRUE)
      for (fname in cyto_files) {
        url <- paste0(cyto_url, "/", fname)
        dest <- file.path(cyto_data_dir, fname)
        log_message(
          "  Downloading {.path {fname}} ...",
          expr = utils::download.file(
            url = url,
            destfile = dest,
            mode = "wb",
            quiet = TRUE
          ),
          verbose = verbose
        )
      }
      log_message(
        "CytoTRACE2 data cached at {.path {cyto_data_dir}}",
        message_type = "success",
        verbose = verbose
      )
    } else {
      log_message(
        "Using cached CytoTRACE2 data from {.path {cyto_data_dir}}",
        verbose = verbose
      )
    }

    cyto_cache <- list(
      data_dir = cyto_data_dir,
      files = cyto_files,
      version = cyto_version
    )
    R.cache::saveCache(
      cyto_cache,
      key = cyto_cache_key,
      comment = paste0(
        cyto_version,
        " nterm:",
        length(cyto_files),
        "|CytoTRACE2-CytoTRACE2"
      )
    )
    db_list[["CytoTRACE2"]] <- cyto_cache
  }
  db_list
}
